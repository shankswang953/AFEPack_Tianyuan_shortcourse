"""Record accepted airfoil states and render rollout figures.

This import-only module writes caller-supplied CSV, GIF, and PNG paths,
normally below `output/policy_rollout/`. Rendering requires Matplotlib and
ImageIO/Pillow support; use `render_trajectory.py` as the command-line entry.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

from .geometry import Action, Airfoil, read_airfoil


FIELDNAMES = [
    "frame",
    "accepted_step",
    "loss",
    "reward",
    "source",
    "action_surface",
    "action_center",
    "action_shift",
    "point_surface",
    "point_index",
    "x",
    "y",
]


class TrajectoryRecorder:
    """Store every accepted airfoil state in one long-format CSV file."""

    def __init__(self, filename: Path) -> None:
        self.filename = filename
        self.frame = 0

    def _write_state(
        self,
        stream,
        *,
        airfoil: Airfoil,
        accepted_step: int,
        loss: float,
        reward: float,
        source: str,
        action: Action | None,
    ) -> None:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        for surface, points in (
            ("U", airfoil.upper),
            ("L", airfoil.lower),
        ):
            for index, point in enumerate(points):
                writer.writerow(
                    {
                        "frame": self.frame,
                        "accepted_step": accepted_step,
                        "loss": f"{loss:.16g}",
                        "reward": f"{reward:.16g}",
                        "source": source,
                        "action_surface": (
                            "" if action is None else action.surface
                        ),
                        "action_center": (
                            "" if action is None else f"{action.center:.16g}"
                        ),
                        "action_shift": (
                            "" if action is None else f"{action.shift:.16g}"
                        ),
                        "point_surface": surface,
                        "point_index": index,
                        "x": f"{point.x:.16g}",
                        "y": f"{point.y:.16g}",
                    }
                )

    def start(self, data_file: Path, *, loss: float) -> None:
        self.filename.parent.mkdir(parents=True, exist_ok=True)
        self.frame = 0
        with self.filename.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
            writer.writeheader()
            self._write_state(
                stream,
                airfoil=read_airfoil(data_file),
                accepted_step=0,
                loss=loss,
                reward=0.0,
                source="initial",
                action=None,
            )

    def append(
        self,
        data_file: Path,
        *,
        accepted_step: int,
        loss: float,
        reward: float,
        source: str,
        action: Action,
    ) -> None:
        self.frame += 1
        with self.filename.open("a", newline="") as stream:
            self._write_state(
                stream,
                airfoil=read_airfoil(data_file),
                accepted_step=accepted_step,
                loss=loss,
                reward=reward,
                source=source,
                action=action,
            )


@dataclass(frozen=True)
class TrajectoryFrame:
    frame: int
    accepted_step: int
    loss: float
    reward: float
    source: str
    action: Action | None
    upper_x: list[float]
    upper_y: list[float]
    lower_x: list[float]
    lower_y: list[float]


def read_trajectory(filename: Path) -> list[TrajectoryFrame]:
    grouped: dict[int, list[dict[str, str]]] = {}
    with filename.open(newline="") as stream:
        for row in csv.DictReader(stream):
            grouped.setdefault(int(row["frame"]), []).append(row)

    frames: list[TrajectoryFrame] = []
    for frame_index in sorted(grouped):
        rows = grouped[frame_index]
        first = rows[0]
        upper = sorted(
            (
                row
                for row in rows
                if row["point_surface"] == "U"
            ),
            key=lambda row: int(row["point_index"]),
        )
        lower = sorted(
            (
                row
                for row in rows
                if row["point_surface"] == "L"
            ),
            key=lambda row: int(row["point_index"]),
        )
        action = None
        if first["action_surface"]:
            action = Action(
                surface=first["action_surface"],
                center=float(first["action_center"]),
                shift=float(first["action_shift"]),
            )
        frames.append(
            TrajectoryFrame(
                frame=frame_index,
                accepted_step=int(first["accepted_step"]),
                loss=float(first["loss"]),
                reward=float(first["reward"]),
                source=first["source"],
                action=action,
                upper_x=[float(row["x"]) for row in upper],
                upper_y=[float(row["y"]) for row in upper],
                lower_x=[float(row["x"]) for row in lower],
                lower_y=[float(row["y"]) for row in lower],
            )
        )
    if not frames:
        raise ValueError(f"trajectory contains no frames: {filename}")
    return frames


def render_trajectory(
    history_file: Path,
    target_file: Path,
    animation_file: Path,
    final_figure: Path,
    *,
    frames_per_second: int = 2,
    animation_dpi: int = 200,
) -> None:
    """Create a GIF and a final PNG from the recorded CSV."""

    cache_root = animation_file.parent / ".plot_cache"
    matplotlib_cache = cache_root / "matplotlib"
    font_cache = cache_root / "fontconfig"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    font_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation, PillowWriter

    frames = read_trajectory(history_file)
    target = read_airfoil(target_file)
    animation_file.parent.mkdir(parents=True, exist_ok=True)
    final_figure.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "legend.fontsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
        }
    )
    figure, axis = plt.subplots(figsize=(8.0, 6.0), dpi=150)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")
    axis.set_xlim(-0.05, 1.05)
    axis.set_ylim(-0.55, 0.55)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"$x/c$")
    axis.set_ylabel(r"$y/c$")
    axis.grid(True, color="#d9e0e6", linewidth=0.8)

    target_upper, = axis.plot(
        [point.x for point in target.upper],
        [point.y for point in target.upper],
        color="#1b9e77",
        linewidth=2.4,
        linestyle="--",
        label="target",
    )
    axis.plot(
        [point.x for point in target.lower],
        [point.y for point in target.lower],
        color="#1b9e77",
        linewidth=2.4,
        linestyle="--",
    )
    current_upper, = axis.plot(
        [],
        [],
        color="#2166ac",
        linewidth=2.8,
        label="current",
    )
    current_lower, = axis.plot(
        [],
        [],
        color="#2166ac",
        linewidth=2.8,
    )
    fixed_points = axis.scatter(
        [0.0, 1.0],
        [0.0, 0.0],
        s=42,
        color="#111111",
        zorder=5,
        label="fixed endpoints",
    )
    action_marker, = axis.plot(
        [],
        [],
        marker="D",
        markersize=7,
        linestyle="None",
        color="#d95f02",
        label="accepted update",
    )
    title = axis.set_title("")
    axis.legend(
        handles=[
            current_upper,
            target_upper,
            fixed_points,
            action_marker,
        ],
        loc="upper right",
        frameon=False,
    )

    def draw(frame_number: int):
        frame = frames[frame_number]
        current_upper.set_data(frame.upper_x, frame.upper_y)
        current_lower.set_data(frame.lower_x, frame.lower_y)
        if frame.action is None:
            action_marker.set_data([], [])
        else:
            x_values = np.asarray(
                (
                    frame.upper_x
                    if frame.action.surface == "U"
                    else frame.lower_x
                )
            )
            y_values = np.asarray(
                (
                    frame.upper_y
                    if frame.action.surface == "U"
                    else frame.lower_y
                )
            )
            marker_y = float(
                np.interp(frame.action.center, x_values, y_values)
            )
            action_marker.set_data([frame.action.center], [marker_y])
        title.set_text(
            f"Accepted step {frame.accepted_step}  |  "
            f"data MSE = {frame.loss:.6f}"
        )
        return (
            current_upper,
            current_lower,
            action_marker,
            title,
        )

    animation = FuncAnimation(
        figure,
        draw,
        frames=len(frames),
        interval=1000 / max(1, frames_per_second),
        blit=False,
        repeat=True,
    )
    animation.save(
        animation_file,
        writer=PillowWriter(fps=max(1, frames_per_second)),
        dpi=animation_dpi,
    )
    draw(len(frames) - 1)
    figure.savefig(
        final_figure,
        dpi=180,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)
