#!/usr/bin/env python3

"""Render synchronized real-boundary and digital-twin evolution.

Command-line usage:
    python3 render_boundary_twin_evolution.py [--max-frames N] [--fps N]
        [--episode N] [--output FILE] [--rebuild-mesh-cache] [--mesh-only]

Requires `output/policy_rollout/data_history.csv` and `output/replay.jsonl`.
The default outputs are `output/policy_rollout/boundary_twin_evolution.gif`
and its final PNG; `--mesh-only` selects the boundary-mesh filenames, and
`--output` overrides the GIF path.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from dt_airfoil.geometry import (
    Action,
    Airfoil,
    Point,
    read_airfoil,
    write_airfoil,
)
from dt_airfoil.learning import TwinPredictor, train_twin
from dt_airfoil.replay import TransitionRecord, load_replay
from dt_airfoil.trajectory import TrajectoryFrame, read_trajectory


@dataclass(frozen=True)
class TwinSnapshot:
    trajectory: TrajectoryFrame
    samples: int
    trial_count: int
    upper_score: np.ndarray
    lower_score: np.ndarray
    next_action: Action | None
    next_score: float | None
    next_predicted_reward: float | None
    mesh_nodes: np.ndarray
    mesh_elements: np.ndarray


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct twin checkpoints from replay.jsonl and render a "
            "synchronized boundary/action-ranking GIF."
        )
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=41,
        help="maximum evenly spaced accepted states to render",
    )
    parser.add_argument("--fps", type=int, default=2)
    parser.add_argument(
        "--episode",
        type=int,
        default=None,
        help="twin-controller episode; default is the latest one",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="output GIF path",
    )
    parser.add_argument(
        "--rebuild-mesh-cache",
        action="store_true",
        help="replay all accepted mesh motions even when a cache is available",
    )
    parser.add_argument(
        "--mesh-only",
        action="store_true",
        help="render only the real local boundary mesh panel",
    )
    arguments = parser.parse_args()
    if arguments.max_frames < 2:
        raise SystemExit("--max-frames must be at least two")
    if arguments.fps < 1:
        raise SystemExit("--fps must be positive")
    return arguments


def latest_controller_episode(
    records: list[TransitionRecord],
) -> int:
    episodes = [
        record.episode
        for record in records
        if record.phase == "twin_controller"
    ]
    if not episodes:
        raise ValueError("replay contains no twin-controller episode")
    return max(episodes)


def select_frame_indices(
    frame_count: int,
    maximum: int,
    mandatory: list[int] | None = None,
) -> list[int]:
    required = {
        int(value)
        for value in (mandatory or [])
        if 0 <= value < frame_count
    }
    required.update((0, frame_count - 1))
    remaining_count = max(0, min(frame_count, maximum) - len(required))
    candidates = [
        value for value in range(frame_count) if value not in required
    ]
    if remaining_count and candidates:
        positions = np.rint(
            np.linspace(
                0,
                len(candidates) - 1,
                min(remaining_count, len(candidates)),
            )
        ).astype(int)
        required.update(candidates[int(position)] for position in positions)
    return sorted(required)


def frame_airfoil(frame: TrajectoryFrame) -> Airfoil:
    return Airfoil(
        title=f"accepted state {frame.accepted_step}",
        upper=[
            Point(float(x_value), float(y_value))
            for x_value, y_value in zip(
                frame.upper_x,
                frame.upper_y,
            )
        ],
        lower=[
            Point(float(x_value), float(y_value))
            for x_value, y_value in zip(
                frame.lower_x,
                frame.lower_y,
            )
        ],
    )


def load_mesh_csv(
    nodes_file: Path,
    elements_file: Path,
) -> tuple[np.ndarray, np.ndarray]:
    nodes = np.loadtxt(
        nodes_file,
        delimiter=",",
        skiprows=1,
        ndmin=2,
    )
    elements = np.loadtxt(
        elements_file,
        delimiter=",",
        skiprows=1,
        dtype=np.int64,
        ndmin=2,
    )
    return (
        np.asarray(nodes[:, 1:3], dtype=np.float64),
        np.asarray(elements[:, 1:4], dtype=np.int64),
    )


def run_checked(command: list[str]) -> None:
    result = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): "
            f"{' '.join(command)}\n{result.stdout}"
        )


def reconstruct_mesh_states(
    *,
    root: Path,
    trajectories: list[TrajectoryFrame],
    controller_records: list[TransitionRecord],
    selected_frames: list[int],
    controller_episode: int,
    force: bool,
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Replay the accepted fixed-topology meshes and exact remesh checkpoints."""

    accepted_records = [
        record for record in controller_records if record.accepted
    ]
    if len(accepted_records) + 1 != len(trajectories):
        raise ValueError(
            "mesh replay requires one trajectory frame per accepted action"
        )

    cache = (
        root
        / "output"
        / "policy_rollout"
        / ".plot_cache"
        / "mesh_evolution"
    )
    signature_file = cache / "signature.json"
    replay_file = root / "output" / "replay.jsonl"
    signature = {
        "episode": controller_episode,
        "accepted_records": len(accepted_records),
        "replay_size": replay_file.stat().st_size,
        "selected_frames": selected_frames,
    }
    mesh_files = {
        frame_index: cache / f"mesh_{frame_index:04d}.npz"
        for frame_index in selected_frames
    }
    cache_matches = False
    if not force and signature_file.exists():
        try:
            cache_matches = (
                json.loads(signature_file.read_text()) == signature
                and all(path.exists() for path in mesh_files.values())
            )
        except (OSError, json.JSONDecodeError):
            cache_matches = False
    if cache_matches:
        return {
            frame_index: (
                np.load(filename)["nodes"],
                np.load(filename)["elements"],
            )
            for frame_index, filename in mesh_files.items()
        }

    if cache.exists():
        shutil.rmtree(cache)
    cache.mkdir(parents=True)
    work = cache / "work"
    work.mkdir()

    initial = root / "output" / "reference" / "initial"
    current_mesh = work / "mesh_current.mesh"
    shutil.copy2(initial / "mesh_smoothed.mesh", current_mesh)
    initial_nodes, initial_elements = load_mesh_csv(
        initial / "mesh_smoothed_nodes.csv",
        initial / "mesh_smoothed_elements.csv",
    )
    if 0 in mesh_files:
        np.savez_compressed(
            mesh_files[0],
            nodes=initial_nodes,
            elements=initial_elements,
        )

    accepted_archives = sorted(
        (root / "output" / "history").glob("accepted_*.dat")
    )
    if not accepted_archives:
        raise ValueError("accepted mesh history was not found")
    maximum_global_step = max(
        int(path.stem.rsplit("_", 1)[1])
        for path in accepted_archives
    )
    global_offset = maximum_global_step - len(accepted_records)

    generate = root / "backend" / "generate_airfoil_geometry"
    move = root / "backend" / "move_and_smooth"
    history = root / "output" / "history"
    selected = set(selected_frames)

    for accepted_index, record in enumerate(
        accepted_records,
        start=1,
    ):
        step = work / "step"
        if step.exists():
            shutil.rmtree(step)
        step.mkdir()

        if record.remeshed:
            global_step = global_offset + accepted_index
            stem = history / f"remeshed_{global_step:04d}"
            archived_mesh = stem.with_suffix(".mesh")
            archived_nodes = Path(f"{stem}_nodes.csv")
            archived_elements = Path(f"{stem}_elements.csv")
            if not (
                archived_mesh.exists()
                and archived_nodes.exists()
                and archived_elements.exists()
            ):
                raise FileNotFoundError(
                    f"remesh checkpoint {stem.name} is incomplete"
                )
            shutil.copy2(archived_mesh, current_mesh)
            nodes, elements = load_mesh_csv(
                archived_nodes,
                archived_elements,
            )
        else:
            current_dat = step / "current.dat"
            next_dat = step / "next.dat"
            write_airfoil(
                current_dat,
                frame_airfoil(trajectories[accepted_index - 1]),
            )
            write_airfoil(
                next_dat,
                frame_airfoil(trajectories[accepted_index]),
            )
            run_checked(
                [
                    str(generate),
                    str(current_dat),
                    str(step),
                    "96",
                    "0.35",
                    "0.12",
                    "0.0",
                    "0.0",
                    str(next_dat),
                ]
            )
            run_checked(
                [
                    str(move),
                    str(current_mesh),
                    str(step / "boundary_initial.dat"),
                    str(step / "boundary_moved.dat"),
                    str(step),
                    "200",
                    "0.45",
                ]
            )
            shutil.copy2(step / "mesh_smoothed.mesh", current_mesh)
            nodes, elements = load_mesh_csv(
                step / "mesh_smoothed_nodes.csv",
                step / "mesh_smoothed_elements.csv",
            )

        if accepted_index in selected:
            np.savez_compressed(
                mesh_files[accepted_index],
                nodes=nodes,
                elements=elements,
            )

    signature_file.write_text(json.dumps(signature, indent=2))
    shutil.rmtree(work)
    return {
        frame_index: (
            np.load(filename)["nodes"],
            np.load(filename)["elements"],
        )
        for frame_index, filename in mesh_files.items()
    }


def controller_partition(
    records: list[TransitionRecord],
    episode: int,
) -> tuple[list[TransitionRecord], list[TransitionRecord]]:
    indices = [
        index
        for index, record in enumerate(records)
        if record.phase == "twin_controller"
        and record.episode == episode
    ]
    if not indices:
        raise ValueError(f"controller episode {episode} was not found")
    if indices != list(range(indices[0], indices[-1] + 1)):
        raise ValueError("controller records are not contiguous in replay")
    return records[: indices[0]], records[indices[0] : indices[-1] + 1]


def frame_trial_positions(
    controller_records: list[TransitionRecord],
) -> list[int]:
    positions = [0]
    for trial_count, record in enumerate(controller_records, start=1):
        if record.accepted:
            positions.append(trial_count)
    return positions


def action_grid(
    center_count: int = 21,
    magnitudes: tuple[float, ...] = (0.005, 0.01, 0.02, 0.03, 0.04),
) -> tuple[np.ndarray, list[Action], list[Action]]:
    centers = np.linspace(0.0, 1.0, center_count)
    upper = [
        Action("U", float(center), -magnitude)
        for center in centers
        for magnitude in magnitudes
    ]
    lower = [
        Action("L", float(center), magnitude)
        for center in centers
        for magnitude in magnitudes
    ]
    return centers, upper, lower


def state_at_frame(
    frame_index: int,
    controller_records: list[TransitionRecord],
) -> np.ndarray:
    if frame_index == 0:
        return np.asarray(
            controller_records[0].state_before,
            dtype=np.float32,
        )
    accepted = [
        record
        for record in controller_records
        if record.accepted
    ]
    return np.asarray(
        accepted[frame_index - 1].state_after_trial,
        dtype=np.float32,
    )


def normalized_scores(
    upper: np.ndarray,
    lower: np.ndarray,
    selected: float | None,
) -> tuple[np.ndarray, np.ndarray, float | None]:
    combined = np.concatenate((upper, lower))
    low = float(np.min(combined))
    high = float(np.max(combined))
    span = high - low
    if span < 1.0e-14:
        return (
            np.full_like(upper, 0.5),
            np.full_like(lower, 0.5),
            0.5 if selected is not None else None,
        )
    selected_score = (
        None
        if selected is None
        else float(np.clip((selected - low) / span, 0.0, 1.0))
    )
    return (
        (upper - low) / span,
        (lower - low) / span,
        selected_score,
    )


def local_mesh_segments(
    nodes: np.ndarray,
    elements: np.ndarray,
) -> np.ndarray:
    triangles = np.asarray(elements, dtype=np.int64)
    edges = np.concatenate(
        (
            triangles[:, [0, 1]],
            triangles[:, [1, 2]],
            triangles[:, [2, 0]],
        ),
        axis=0,
    )
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)
    segments = nodes[edges]
    midpoints = np.mean(segments, axis=1)
    local = (
        (midpoints[:, 0] >= -0.10)
        & (midpoints[:, 0] <= 1.10)
        & (midpoints[:, 1] >= -0.58)
        & (midpoints[:, 1] <= 0.58)
    )
    return segments[local]


def reconstruct_snapshots(
    *,
    root: Path,
    records: list[TransitionRecord],
    controller_episode: int,
    trajectories: list[TrajectoryFrame],
    selected_frames: list[int],
    mesh_states: dict[int, tuple[np.ndarray, np.ndarray]],
    configuration: dict,
) -> tuple[np.ndarray, list[TwinSnapshot]]:
    prefix, controller = controller_partition(
        records,
        controller_episode,
    )
    positions = frame_trial_positions(controller)
    if len(positions) != len(trajectories):
        raise ValueError(
            "trajectory/replay mismatch: "
            f"{len(trajectories)} states but "
            f"{len(positions)} accepted-state positions"
        )

    centers, upper_actions, lower_actions = action_grid()
    magnitudes_per_center = len(upper_actions) // len(centers)
    config_arguments = configuration.get("arguments", {})
    stage_seeds = configuration.get("stage_seeds", {})
    controller_seed = int(stage_seeds.get("controller", 2026))
    controller_epochs = int(
        config_arguments.get("controller_twin_epochs", 150)
    )
    optimization_seed = int(
        stage_seeds.get("twin_optimization", 23)
    )
    optimization_episodes = int(
        config_arguments.get("optimization_episodes", 4)
    )

    cache = (
        root
        / "output"
        / "policy_rollout"
        / ".plot_cache"
        / "twin_evolution"
    )
    cache.mkdir(parents=True, exist_ok=True)
    checkpoint = cache / "snapshot.pt"

    snapshots: list[TwinSnapshot] = []
    for frame_index in selected_frames:
        trial_count = positions[frame_index]
        training_records = prefix + controller[:trial_count]
        if trial_count == 0:
            epochs = 250
            seed = optimization_seed + optimization_episodes
        else:
            epochs = controller_epochs
            seed = controller_seed + trial_count
        train_twin(
            training_records,
            checkpoint,
            epochs=epochs,
            seed=seed,
        )
        predictor = TwinPredictor(checkpoint)
        state = state_at_frame(frame_index, controller)
        upper_prediction = predictor.predict(state, upper_actions).reshape(
            len(centers),
            magnitudes_per_center,
        )
        lower_prediction = predictor.predict(state, lower_actions).reshape(
            len(centers),
            magnitudes_per_center,
        )
        upper_best = np.max(upper_prediction, axis=1)
        lower_best = np.max(lower_prediction, axis=1)

        next_action = None
        next_prediction = None
        if trial_count < len(controller):
            record = controller[trial_count]
            next_action = record.action
            next_prediction = float(
                predictor.predict(state, [next_action])[0]
            )
        upper_score, lower_score, next_score = normalized_scores(
            upper_best,
            lower_best,
            next_prediction,
        )
        snapshots.append(
            TwinSnapshot(
                trajectory=trajectories[frame_index],
                samples=len(training_records),
                trial_count=trial_count,
                upper_score=upper_score,
                lower_score=lower_score,
                next_action=next_action,
                next_score=next_score,
                next_predicted_reward=next_prediction,
                mesh_nodes=mesh_states[frame_index][0],
                mesh_elements=mesh_states[frame_index][1],
            )
        )
    return centers, snapshots


def render_mesh_only(
    *,
    root: Path,
    trajectories: list[TrajectoryFrame],
    selected_frames: list[int],
    mesh_states: dict[int, tuple[np.ndarray, np.ndarray]],
    target_file: Path,
    animation_file: Path,
    final_figure: Path,
    fps: int,
) -> None:
    cache_root = (
        root / "output" / "policy_rollout" / ".plot_cache"
    )
    os.environ.setdefault(
        "MPLCONFIGDIR",
        str(cache_root / "matplotlib"),
    )
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.collections import LineCollection

    target = read_airfoil(target_file)
    animation_file.parent.mkdir(parents=True, exist_ok=True)
    final_figure.parent.mkdir(parents=True, exist_ok=True)

    navy = "#17324D"
    blue = "#2F6B9A"
    teal = "#2A9D8F"
    red = "#D95D4F"

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "legend.fontsize": 11,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )
    figure, axis = plt.subplots(figsize=(7.5, 6.2), dpi=120)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")
    axis.set_xlim(-0.05, 1.05)
    axis.set_ylim(-0.55, 0.55)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"$x/c$")
    axis.set_ylabel(r"$y/c$")
    axis.grid(False)

    mesh_collection = LineCollection(
        [],
        colors=blue,
        linewidths=0.55,
        alpha=0.78,
        zorder=1,
    )
    axis.add_collection(mesh_collection)
    target_upper, = axis.plot(
        [point.x for point in target.upper],
        [point.y for point in target.upper],
        color=teal,
        linewidth=2.3,
        linestyle="--",
        label="target",
        zorder=3,
    )
    axis.plot(
        [point.x for point in target.lower],
        [point.y for point in target.lower],
        color=teal,
        linewidth=2.3,
        linestyle="--",
        zorder=3,
    )
    current_upper, = axis.plot(
        [],
        [],
        color=navy,
        linewidth=2.9,
        label="current",
        zorder=4,
    )
    current_lower, = axis.plot(
        [],
        [],
        color=navy,
        linewidth=2.9,
        zorder=4,
    )
    accepted_marker, = axis.plot(
        [],
        [],
        marker="D",
        markersize=7,
        linestyle="None",
        color=red,
        label="last accepted update",
        zorder=5,
    )
    status = axis.text(
        0.02,
        0.03,
        "",
        transform=axis.transAxes,
        color=navy,
        va="bottom",
    )
    axis.legend(
        handles=[current_upper, target_upper, accepted_marker],
        loc="upper right",
        frameon=False,
    )
    title = axis.set_title("")

    def draw(frame_number: int):
        frame_index = selected_frames[frame_number]
        frame = trajectories[frame_index]
        nodes, elements = mesh_states[frame_index]
        mesh_collection.set_segments(
            local_mesh_segments(nodes, elements)
        )
        current_upper.set_data(frame.upper_x, frame.upper_y)
        current_lower.set_data(frame.lower_x, frame.lower_y)
        if frame.action is None:
            accepted_marker.set_data([], [])
        else:
            x_values = np.asarray(
                frame.upper_x
                if frame.action.surface == "U"
                else frame.lower_x
            )
            y_values = np.asarray(
                frame.upper_y
                if frame.action.surface == "U"
                else frame.lower_y
            )
            accepted_marker.set_data(
                [frame.action.center],
                [
                    float(
                        np.interp(
                            frame.action.center,
                            x_values,
                            y_values,
                        )
                    )
                ],
            )
        status.set_text(frame.source)
        title.set_text(
            f"Accepted step {frame.accepted_step}  |  "
            f"data MSE = {frame.loss:.3e}"
        )
        return (
            mesh_collection,
            current_upper,
            current_lower,
            accepted_marker,
            status,
            title,
        )

    figure.tight_layout()
    animation = FuncAnimation(
        figure,
        draw,
        frames=len(selected_frames),
        interval=1000 / fps,
        blit=False,
        repeat=True,
    )
    animation.save(
        animation_file,
        writer=PillowWriter(fps=fps),
        dpi=120,
    )
    draw(len(selected_frames) - 1)
    figure.savefig(
        final_figure,
        dpi=180,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)


def render(
    *,
    root: Path,
    centers: np.ndarray,
    snapshots: list[TwinSnapshot],
    target_file: Path,
    animation_file: Path,
    final_figure: Path,
    fps: int,
) -> None:
    cache_root = (
        root / "output" / "policy_rollout" / ".plot_cache"
    )
    os.environ.setdefault(
        "MPLCONFIGDIR",
        str(cache_root / "matplotlib"),
    )
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.collections import LineCollection

    target = read_airfoil(target_file)
    animation_file.parent.mkdir(parents=True, exist_ok=True)
    final_figure.parent.mkdir(parents=True, exist_ok=True)

    navy = "#17324D"
    blue = "#2F6B9A"
    teal = "#2A9D8F"
    gold = "#E9A23B"
    red = "#D95D4F"
    grid = "#D8DEE6"

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "legend.fontsize": 11,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )
    figure, (shape_axis, twin_axis) = plt.subplots(
        1,
        2,
        figsize=(12.0, 5.2),
        dpi=110,
        gridspec_kw={"width_ratios": (1.05, 0.95)},
    )
    figure.patch.set_facecolor("white")
    for axis in (shape_axis, twin_axis):
        axis.set_facecolor("white")
        axis.grid(True, color=grid, linewidth=0.8)

    shape_axis.set_title("REAL: local boundary mesh")
    shape_axis.set_xlim(-0.05, 1.05)
    shape_axis.set_ylim(-0.55, 0.55)
    shape_axis.set_aspect("equal", adjustable="box")
    shape_axis.set_xlabel(r"$x/c$")
    shape_axis.set_ylabel(r"$y/c$")
    shape_axis.grid(False)
    mesh_collection = LineCollection(
        [],
        colors=blue,
        linewidths=0.48,
        alpha=0.72,
        zorder=1,
    )
    shape_axis.add_collection(mesh_collection)
    target_upper, = shape_axis.plot(
        [point.x for point in target.upper],
        [point.y for point in target.upper],
        color=teal,
        linewidth=2.3,
        linestyle="--",
        label="target",
        zorder=3,
    )
    shape_axis.plot(
        [point.x for point in target.lower],
        [point.y for point in target.lower],
        color=teal,
        linewidth=2.3,
        linestyle="--",
        zorder=3,
    )
    current_upper, = shape_axis.plot(
        [],
        [],
        color=navy,
        linewidth=2.8,
        label="current",
        zorder=4,
    )
    current_lower, = shape_axis.plot(
        [],
        [],
        color=navy,
        linewidth=2.8,
        zorder=4,
    )
    accepted_marker, = shape_axis.plot(
        [],
        [],
        marker="D",
        markersize=7,
        linestyle="None",
        color=red,
        label="last accepted update",
        zorder=5,
    )
    shape_status = shape_axis.text(
        0.02,
        0.03,
        "",
        transform=shape_axis.transAxes,
        color=navy,
        va="bottom",
    )
    shape_axis.legend(
        handles=[current_upper, target_upper, accepted_marker],
        loc="upper right",
        frameon=False,
    )

    twin_axis.set_title("DIGITAL: action ranking after update")
    twin_axis.set_xlim(-0.02, 1.02)
    twin_axis.set_ylim(-0.05, 1.05)
    twin_axis.set_xlabel(r"action center $x/c$")
    twin_axis.set_ylabel("relative predicted reward")
    upper_line, = twin_axis.plot(
        [],
        [],
        color=blue,
        linewidth=2.5,
        marker="o",
        markersize=3.5,
        label="upper surface",
    )
    lower_line, = twin_axis.plot(
        [],
        [],
        color=gold,
        linewidth=2.5,
        marker="s",
        markersize=3.5,
        label="lower surface",
    )
    proposal_marker, = twin_axis.plot(
        [],
        [],
        marker="D",
        markersize=8,
        linestyle="None",
        color=red,
        label="next real action",
    )
    twin_status = twin_axis.text(
        0.02,
        0.03,
        "",
        transform=twin_axis.transAxes,
        color=navy,
        va="bottom",
    )
    twin_axis.legend(loc="upper right", frameon=False)
    main_title = figure.suptitle("", color=navy, fontweight="bold")

    def draw(frame_number: int):
        snapshot = snapshots[frame_number]
        frame = snapshot.trajectory
        mesh_collection.set_segments(
            local_mesh_segments(
                snapshot.mesh_nodes,
                snapshot.mesh_elements,
            )
        )
        current_upper.set_data(frame.upper_x, frame.upper_y)
        current_lower.set_data(frame.lower_x, frame.lower_y)
        if frame.action is None:
            accepted_marker.set_data([], [])
        else:
            x_values = np.asarray(
                frame.upper_x
                if frame.action.surface == "U"
                else frame.lower_x
            )
            y_values = np.asarray(
                frame.upper_y
                if frame.action.surface == "U"
                else frame.lower_y
            )
            marker_y = float(
                np.interp(frame.action.center, x_values, y_values)
            )
            accepted_marker.set_data(
                [frame.action.center],
                [marker_y],
            )
        shape_status.set_text(frame.source)

        upper_line.set_data(centers, snapshot.upper_score)
        lower_line.set_data(centers, snapshot.lower_score)
        if snapshot.next_action is None or snapshot.next_score is None:
            proposal_marker.set_data([], [])
            twin_status.set_text("controller converged")
        else:
            proposal_marker.set_data(
                [snapshot.next_action.center],
                [snapshot.next_score],
            )
            twin_status.set_text(
                f"next: {snapshot.next_action.surface}  "
                f"x={snapshot.next_action.center:.3f}  "
                f"dy={snapshot.next_action.shift:+.4f}\n"
                f"predicted reward="
                f"{snapshot.next_predicted_reward:+.2e}"
            )

        main_title.set_text(
            f"Accepted step {frame.accepted_step}  |  "
            f"data MSE = {frame.loss:.3e}  |  "
            f"twin samples = {snapshot.samples}"
        )
        return (
            mesh_collection,
            current_upper,
            current_lower,
            accepted_marker,
            upper_line,
            lower_line,
            proposal_marker,
            shape_status,
            twin_status,
            main_title,
        )

    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    animation = FuncAnimation(
        figure,
        draw,
        frames=len(snapshots),
        interval=1000 / fps,
        blit=False,
        repeat=True,
    )
    animation.save(
        animation_file,
        writer=PillowWriter(fps=fps),
        dpi=110,
    )
    draw(len(snapshots) - 1)
    figure.savefig(
        final_figure,
        dpi=180,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)


def main() -> None:
    arguments = parse_arguments()
    root = Path(__file__).resolve().parent
    rollout = root / "output" / "policy_rollout"
    history_file = rollout / "data_history.csv"
    replay_file = root / "output" / "replay.jsonl"
    configuration_file = root / "output" / "run_configuration.json"
    target_file = root / "data" / "target_naca0012.dat"
    if not history_file.exists() or not replay_file.exists():
        raise SystemExit(
            "run run_experiment.py before rendering twin evolution"
        )
    configuration = (
        json.loads(configuration_file.read_text())
        if configuration_file.exists()
        else {}
    )
    records = load_replay(replay_file)
    episode = (
        latest_controller_episode(records)
        if arguments.episode is None
        else arguments.episode
    )
    trajectories = read_trajectory(history_file)
    _, controller_records = controller_partition(records, episode)
    accepted_records = [
        record for record in controller_records if record.accepted
    ]
    remesh_frames = [
        frame_index
        for frame_index, record in enumerate(
            accepted_records,
            start=1,
        )
        if record.remeshed
    ]
    selected_frames = select_frame_indices(
        len(trajectories),
        arguments.max_frames,
        mandatory=remesh_frames
        + [frame_index - 1 for frame_index in remesh_frames],
    )
    mesh_states = reconstruct_mesh_states(
        root=root,
        trajectories=trajectories,
        controller_records=controller_records,
        selected_frames=selected_frames,
        controller_episode=episode,
        force=arguments.rebuild_mesh_cache,
    )
    animation_file = (
        arguments.output
        if arguments.output is not None
        else rollout
        / (
            "boundary_mesh_evolution.gif"
            if arguments.mesh_only
            else "boundary_twin_evolution.gif"
        )
    )
    final_figure = animation_file.with_name(
        f"{animation_file.stem}_final.png"
    )
    if arguments.mesh_only:
        render_mesh_only(
            root=root,
            trajectories=trajectories,
            selected_frames=selected_frames,
            mesh_states=mesh_states,
            target_file=target_file,
            animation_file=animation_file,
            final_figure=final_figure,
            fps=arguments.fps,
        )
    else:
        centers, snapshots = reconstruct_snapshots(
            root=root,
            records=records,
            controller_episode=episode,
            trajectories=trajectories,
            selected_frames=selected_frames,
            mesh_states=mesh_states,
            configuration=configuration,
        )
        render(
            root=root,
            centers=centers,
            snapshots=snapshots,
            target_file=target_file,
            animation_file=animation_file,
            final_figure=final_figure,
            fps=arguments.fps,
        )
    print(f"controller episode: {episode}")
    print(f"rendered snapshots: {len(selected_frames)}")
    print(f"animation: {animation_file}")
    print(f"final PNG: {final_figure}")


if __name__ == "__main__":
    main()
