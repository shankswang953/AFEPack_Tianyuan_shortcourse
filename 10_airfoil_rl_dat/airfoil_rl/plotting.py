"""Write training plots, rollout tables, `dat` histories, and animations.

This import-only module writes only caller-supplied paths, normally below
`output/` and `output/evaluation/`. Matplotlib and ImageIO/Pillow support
are needed for PNG/GIF rendering.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np

from .geometry import Action, Airfoil


def render_training_history(
    rows: list[dict[str, float | int]],
    filename: Path,
) -> None:
    cache = filename.parent / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    filename.parent.mkdir(parents=True, exist_ok=True)
    episodes = [int(row["episode"]) for row in rows]
    returns = [float(row["return"]) for row in rows]
    final_losses = [float(row["final_loss"]) for row in rows]
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), dpi=140)
    axes[0].plot(episodes, returns, color="#2F6B9A", linewidth=1.2)
    axes[0].set_xlabel("training episode")
    axes[0].set_ylabel("episode return")
    axes[0].grid(True, color="#D8DEE6", linewidth=0.7)
    axes[1].semilogy(
        episodes,
        np.maximum(final_losses, 1.0e-12),
        color="#D95D4F",
        linewidth=1.2,
    )
    axes[1].set_xlabel("training episode")
    axes[1].set_ylabel("terminal data MSE")
    axes[1].grid(True, color="#D8DEE6", linewidth=0.7)
    figure.tight_layout()
    figure.savefig(filename, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def write_rollout_csv(
    filename: Path,
    actions: list[Action],
    losses: list[float],
    rewards: list[float],
) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "step",
                "loss",
                "reward",
                "surface",
                "center",
                "shift",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "step": 0,
                "loss": f"{losses[0]:.16g}",
                "reward": "0",
                "surface": "",
                "center": "",
                "shift": "",
            }
        )
        for step, (action, loss, reward) in enumerate(
            zip(actions, losses[1:], rewards),
            start=1,
        ):
            writer.writerow(
                {
                    "step": step,
                    "loss": f"{loss:.16g}",
                    "reward": f"{reward:.16g}",
                    "surface": action.surface,
                    "center": f"{action.center:.16g}",
                    "shift": f"{action.shift:.16g}",
                }
            )


def render_rollout(
    *,
    states: list[Airfoil],
    actions: list[Action],
    losses: list[float],
    target: Airfoil,
    animation_file: Path,
    final_figure: Path,
    loss_figure: Path,
    fps: int = 3,
    maximum_frames: int = 61,
    animation_dpi: int = 180,
) -> None:
    cache = animation_file.parent / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    animation_file.parent.mkdir(parents=True, exist_ok=True)
    frame_count = min(maximum_frames, len(states))
    selected = sorted(
        set(
            np.rint(
                np.linspace(0, len(states) - 1, frame_count)
            ).astype(int)
        )
    )
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "legend.fontsize": 11,
        }
    )
    figure, axis = plt.subplots(figsize=(8.0, 5.4), dpi=130)
    axis.set_xlim(-0.05, 1.05)
    axis.set_ylim(-0.55, 0.55)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"$x/c$")
    axis.set_ylabel(r"$y/c$")
    axis.grid(True, color="#D8DEE6", linewidth=0.75)
    target_upper, = axis.plot(
        target.x,
        target.upper_y,
        color="#2A9D8F",
        linewidth=2.4,
        linestyle="--",
        label="target",
    )
    axis.plot(
        target.x,
        target.lower_y,
        color="#2A9D8F",
        linewidth=2.4,
        linestyle="--",
    )
    current_upper, = axis.plot(
        [],
        [],
        color="#17324D",
        linewidth=2.9,
        label="RL state",
    )
    current_lower, = axis.plot(
        [],
        [],
        color="#17324D",
        linewidth=2.9,
    )
    action_marker, = axis.plot(
        [],
        [],
        marker="D",
        markersize=7,
        linestyle="None",
        color="#D95D4F",
        label="selected action",
    )
    axis.legend(
        handles=[current_upper, target_upper, action_marker],
        loc="upper right",
        frameon=False,
    )
    status = axis.text(
        0.02,
        0.03,
        "",
        transform=axis.transAxes,
        color="#17324D",
        va="bottom",
    )
    title = axis.set_title("")

    def draw(animation_index: int):
        state_index = selected[animation_index]
        state = states[state_index]
        current_upper.set_data(state.x, state.upper_y)
        current_lower.set_data(state.x, state.lower_y)
        if state_index == 0:
            action_marker.set_data([], [])
            status.set_text("reset: circular initial state")
        else:
            action = actions[state_index - 1]
            values = (
                state.upper_y
                if action.surface == "U"
                else state.lower_y
            )
            marker_y = float(np.interp(action.center, state.x, values))
            action_marker.set_data([action.center], [marker_y])
            status.set_text(
                f"{action.surface}  x={action.center:.2f}  "
                f"dy={action.shift:+.3f}"
            )
        title.set_text(
            f"Greedy RL rollout: step {state_index}  |  "
            f"data MSE = {losses[state_index]:.3e}"
        )
        return (
            current_upper,
            current_lower,
            action_marker,
            status,
            title,
        )

    figure.tight_layout()
    animation = FuncAnimation(
        figure,
        draw,
        frames=len(selected),
        interval=1000 / fps,
        blit=False,
        repeat=True,
    )
    animation.save(
        animation_file,
        writer=PillowWriter(fps=fps),
        dpi=animation_dpi,
    )
    draw(len(selected) - 1)
    figure.savefig(
        final_figure,
        dpi=180,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7.2, 4.2), dpi=150)
    axis.semilogy(
        range(len(losses)),
        np.maximum(losses, 1.0e-12),
        color="#2F6B9A",
        linewidth=2.2,
    )
    axis.axhline(
        2.0e-6,
        color="#D95D4F",
        linestyle="--",
        linewidth=1.5,
        label="terminal threshold",
    )
    axis.set_xlabel("RL step")
    axis.set_ylabel("data MSE")
    axis.grid(True, color="#D8DEE6", linewidth=0.7)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(
        loss_figure,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)
