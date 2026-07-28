#!/usr/bin/env python3

"""Plot and animate the barycentric fixed-topology mesh-motion results.

Command-line usage:
    python3 plot_results.py [OUTPUT_DIR] [--gif-dpi N]

OUTPUT_DIR defaults to `output`. The script reads `continuation.csv` and
`snapshots/`, then writes PNG and GIF teaching artifacts below
OUTPUT_DIR/figures/. It does not run AFEPack or EasyMesh.
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

cache = Path(tempfile.gettempdir()) / "afepack-barycentric-plot-cache"
cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(cache))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.collections import LineCollection


def rows(filename: Path) -> list[dict[str, str]]:
    with filename.open(newline="") as stream:
        return list(csv.DictReader(stream))


def read_nodes(filename: Path) -> list[tuple[float, float, int]]:
    return [
        (float(row["x"]), float(row["y"]), int(row["boundary_mark"]))
        for row in rows(filename)
    ]


def read_elements(filename: Path) -> list[tuple[int, int, int]]:
    return [
        (int(row["v0"]), int(row["v1"]), int(row["v2"]))
        for row in rows(filename)
    ]


def read_airfoil(
    filename: Path,
) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    with filename.open() as stream:
        stream.readline()
        counts = stream.readline().split()
        upper_count = int(round(float(counts[0].rstrip("."))))
        lower_count = int(round(float(counts[1].rstrip("."))))
        points = [
            (float(fields[0]), float(fields[1]))
            for line in stream
            if len(fields := line.split()) >= 2
        ]
    if len(points) != upper_count + lower_count:
        raise ValueError(f"invalid airfoil data: {filename}")
    return points[:upper_count], points[upper_count:]


def unique_edges(
    elements: list[tuple[int, int, int]],
) -> list[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for first, second, third in elements:
        for left, right in (
            (first, second),
            (second, third),
            (third, first),
        ):
            edges.add((min(left, right), max(left, right)))
    return sorted(edges)


def mesh_segments(
    nodes: list[tuple[float, float, int]],
    edges: list[tuple[int, int]],
) -> list[list[tuple[float, float]]]:
    return [
        [
            (nodes[first][0], nodes[first][1]),
            (nodes[second][0], nodes[second][1]),
        ]
        for first, second in edges
    ]


def airfoil_segments(
    nodes: list[tuple[float, float, int]],
    elements: list[tuple[int, int, int]],
) -> list[list[tuple[float, float]]]:
    result = []
    for first, second in unique_edges(elements):
        if nodes[first][2] == 3 and nodes[second][2] == 3:
            result.append(
                [
                    (nodes[first][0], nodes[first][1]),
                    (nodes[second][0], nodes[second][1]),
                ]
            )
    return result


def configure_mesh_axis(axis: plt.Axes) -> None:
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(-0.08, 1.08)
    # Keep one physical scale for every panel and include the complete
    # circular starting boundary.  A tighter limit makes the circle look
    # disconnected because its top and bottom are clipped.
    axis.set_ylim(-0.53, 0.53)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)


def selected_rows(
    history: list[dict[str, str]],
) -> list[dict[str, str]]:
    selected = []
    for target in (0.0, 0.25, 0.5, 0.75, 1.0):
        selected.append(
            min(history, key=lambda row: abs(float(row["theta"]) - target))
        )
    return selected


def snapshot_prefix(snapshot_dir: Path, row: dict[str, str]) -> Path:
    return snapshot_dir / f"stage_{int(row['stage']):03d}"


def plot_data_path(
    figure_dir: Path,
    snapshot_dir: Path,
    selected: list[dict[str, str]],
) -> None:
    figure, axis = plt.subplots(figsize=(10.5, 4.5), constrained_layout=True)
    colors = ("#1d4ed8", "#0f766e", "#65a30d", "#d97706", "#be123c")
    for color, row in zip(colors, selected):
        prefix = snapshot_prefix(snapshot_dir, row)
        upper, lower = read_airfoil(Path(str(prefix) + ".dat"))
        theta = float(row["theta"])
        for points in (upper, lower):
            axis.plot(
                [point[0] for point in points],
                [point[1] for point in points],
                color=color,
                linewidth=2.2,
                label=rf"$\theta={theta:.2f}$" if points is upper else None,
            )
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(-0.03, 1.03)
    axis.set_ylim(-0.53, 0.53)
    axis.set_xlabel(r"$x/c$")
    axis.set_ylabel(r"$y/c$")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.legend(loc="upper center", ncol=5, frameon=False)
    figure.savefig(
        figure_dir / "01_barycentric_data_path.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_mesh_path(
    figure_dir: Path,
    snapshot_dir: Path,
    selected: list[dict[str, str]],
) -> None:
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.5, 7.8),
        constrained_layout=True,
    )
    flat_axes = list(axes.flat)
    for axis, row in zip(flat_axes, selected):
        prefix = snapshot_prefix(snapshot_dir, row)
        nodes = read_nodes(Path(str(prefix) + "_nodes.csv"))
        elements = read_elements(Path(str(prefix) + "_elements.csv"))
        edges = unique_edges(elements)
        axis.add_collection(
            LineCollection(
                mesh_segments(nodes, edges),
                colors="#315f86",
                linewidths=0.48,
            )
        )
        axis.add_collection(
            LineCollection(
                airfoil_segments(nodes, elements),
                colors="#111827",
                linewidths=1.8,
            )
        )
        axis.set_title(
            rf"$\theta={float(row['theta']):.2f}$"
            + "\n"
            + rf"$q_{{\min}}={float(row['min_quality']):.3f}$"
        )
        configure_mesh_axis(axis)
    flat_axes[-1].axis("off")
    flat_axes[-1].text(
        0.5,
        0.57,
        "1,781 nodes\n3,402 triangles",
        ha="center",
        va="center",
        transform=flat_axes[-1].transAxes,
        fontsize=18,
    )
    flat_axes[-1].text(
        0.5,
        0.34,
        "identical connectivity\nat every stage",
        ha="center",
        va="center",
        transform=flat_axes[-1].transAxes,
        color="#315f86",
    )
    figure.savefig(
        figure_dir / "02_fixed_topology_mesh_path.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_quality(
    figure_dir: Path,
    history: list[dict[str, str]],
) -> None:
    theta = [float(row["theta"]) for row in history]
    quality = [float(row["min_quality"]) for row in history]
    figure, axis = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    axis.plot(theta, quality, color="#0f766e", linewidth=3.0)
    axis.scatter(theta, quality, color="#0f766e", s=16, zorder=3)
    axis.axhline(
        0.40,
        color="#be123c",
        linewidth=1.8,
        linestyle="--",
        label="conservative remesh threshold",
    )
    axis.axhline(
        0.0,
        color="#111827",
        linewidth=1.2,
        label="inversion limit",
    )
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 0.95)
    axis.set_xlabel(r"barycentric coordinate $\theta$")
    axis.set_ylabel(r"minimum triangle quality $q_{\min}$")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(axis="y", color="#d1d5db", linewidth=0.7)
    axis.legend(frameon=False)
    figure.savefig(
        figure_dir / "03_mesh_quality_along_path.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def equilateral_targets(
    nodes: list[tuple[float, float, int]],
    elements: list[tuple[int, int, int]],
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    targets = []
    for element in elements:
        boundary = [vertex for vertex in element if nodes[vertex][2] == 3]
        interior = [vertex for vertex in element if nodes[vertex][2] == 0]
        if len(boundary) != 2 or len(interior) != 1:
            continue
        first = nodes[boundary[0]]
        second = nodes[boundary[1]]
        current = nodes[interior[0]]
        edge_x = second[0] - first[0]
        edge_y = second[1] - first[1]
        edge_length = (edge_x * edge_x + edge_y * edge_y) ** 0.5
        midpoint = (
            0.5 * (first[0] + second[0]),
            0.5 * (first[1] + second[1]),
        )
        normal = (-edge_y / edge_length, edge_x / edge_length)
        side = (
            1.0
            if (
                (current[0] - midpoint[0]) * normal[0]
                + (current[1] - midpoint[1]) * normal[1]
            )
            >= 0.0
            else -1.0
        )
        height = 0.5 * 3.0**0.5 * edge_length
        target = (
            midpoint[0] + side * height * normal[0],
            midpoint[1] + side * height * normal[1],
        )
        targets.append(((current[0], current[1]), target))
    return targets


def draw_close_mesh(
    axis: plt.Axes,
    nodes: list[tuple[float, float, int]],
    elements: list[tuple[int, int, int]],
    y_limit: float = 0.20,
) -> None:
    edges = unique_edges(elements)
    axis.add_collection(
        LineCollection(
            mesh_segments(nodes, edges),
            colors="#315f86",
            linewidths=0.55,
        )
    )
    axis.add_collection(
        LineCollection(
            airfoil_segments(nodes, elements),
            colors="#111827",
            linewidths=2.0,
        )
    )
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(-0.03, 1.03)
    axis.set_ylim(-y_limit, y_limit)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)


def plot_boundary_quality_comparison(
    figure_dir: Path,
    snapshot_dir: Path,
    final_row: dict[str, str],
) -> None:
    prefix = snapshot_prefix(snapshot_dir, final_row)
    unsmoothed_nodes = read_nodes(
        Path(str(prefix) + "_before_nodes.csv")
    )
    smoothed_nodes = read_nodes(Path(str(prefix) + "_nodes.csv"))
    elements = read_elements(Path(str(prefix) + "_elements.csv"))
    quality = {
        row["stage"]: float(row["min_shape_quality"])
        for row in rows(Path(str(prefix) + "_quality.csv"))
    }

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(15.5, 4.4),
        constrained_layout=True,
    )
    draw_close_mesh(axes[0], unsmoothed_nodes, elements)
    axes[0].set_title(
        "Boundary moved"
        + "\n"
        + rf"$q_{{\min}}={quality['moved_unsmoothed']:.3f}$"
    )

    draw_close_mesh(axes[1], unsmoothed_nodes, elements)
    targets = equilateral_targets(unsmoothed_nodes, elements)
    axes[1].scatter(
        [target[0] for _, target in targets],
        [target[1] for _, target in targets],
        color="#d97706",
        s=12,
        zorder=4,
    )
    for current, target in targets[::4]:
        axes[1].annotate(
            "",
            xy=target,
            xytext=current,
            arrowprops={
                "arrowstyle": "->",
                "color": "#d97706",
                "linewidth": 1.1,
            },
            zorder=4,
        )
    axes[1].set_title("Equilateral third-vertex targets")

    draw_close_mesh(axes[2], smoothed_nodes, elements)
    axes[2].set_title(
        "Quality-aware smoothing"
        + "\n"
        + rf"$q_{{\min}}={quality['smoothed']:.3f}$"
    )
    figure.savefig(
        figure_dir / "04_boundary_quality_smoothing.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def animate_smoothing_mechanism(
    figure_dir: Path,
    snapshot_dir: Path,
    history: list[dict[str, str]],
    gif_dpi: int,
) -> None:
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(15.5, 4.4),
        constrained_layout=True,
    )
    writer = PillowWriter(fps=5)
    with writer.saving(
        figure,
        figure_dir / "barycentric_smoothing_mechanism.gif",
        dpi=gif_dpi,
    ):
        for row in history:
            prefix = snapshot_prefix(snapshot_dir, row)
            before = read_nodes(Path(str(prefix) + "_before_nodes.csv"))
            after = read_nodes(Path(str(prefix) + "_nodes.csv"))
            elements = read_elements(Path(str(prefix) + "_elements.csv"))
            quality = {
                item["stage"]: float(item["min_shape_quality"])
                for item in rows(Path(str(prefix) + "_quality.csv"))
            }

            for axis in axes:
                axis.clear()
            draw_close_mesh(axes[0], before, elements, y_limit=0.53)
            axes[0].set_title(
                "Boundary moved"
                + "\n"
                + rf"$q_{{\min}}={quality['moved_unsmoothed']:.3f}$"
            )

            draw_close_mesh(axes[1], before, elements, y_limit=0.53)
            targets = equilateral_targets(before, elements)
            axes[1].scatter(
                [target[0] for _, target in targets],
                [target[1] for _, target in targets],
                color="#d97706",
                s=10,
                zorder=4,
            )
            for current, target in targets[::4]:
                axes[1].annotate(
                    "",
                    xy=target,
                    xytext=current,
                    arrowprops={
                        "arrowstyle": "->",
                        "color": "#d97706",
                        "linewidth": 1.0,
                    },
                    zorder=4,
                )
            axes[1].set_title(
                "Equilateral third-vertex targets"
                + "\n"
                + rf"$\theta={float(row['theta']):.3f}$"
            )

            draw_close_mesh(axes[2], after, elements, y_limit=0.53)
            axes[2].set_title(
                "Quality-aware smoothing"
                + "\n"
                + rf"$q_{{\min}}={quality['smoothed']:.3f}$"
            )
            writer.grab_frame()
    plt.close(figure)


def animate_mesh(
    figure_dir: Path,
    snapshot_dir: Path,
    history: list[dict[str, str]],
    gif_dpi: int,
) -> None:
    first_prefix = snapshot_prefix(snapshot_dir, history[0])
    elements = read_elements(Path(str(first_prefix) + "_elements.csv"))
    edges = unique_edges(elements)
    figure, axis = plt.subplots(figsize=(9.2, 5.0), constrained_layout=True)
    configure_mesh_axis(axis)
    mesh_lines = LineCollection([], colors="#315f86", linewidths=0.50)
    boundary_lines = LineCollection([], colors="#111827", linewidths=1.8)
    axis.add_collection(mesh_lines)
    axis.add_collection(boundary_lines)
    title = axis.set_title("")

    writer = PillowWriter(fps=6)
    with writer.saving(
        figure,
        figure_dir / "barycentric_fixed_topology.gif",
        dpi=gif_dpi,
    ):
        for row in history:
            prefix = snapshot_prefix(snapshot_dir, row)
            nodes = read_nodes(Path(str(prefix) + "_nodes.csv"))
            mesh_lines.set_segments(mesh_segments(nodes, edges))
            boundary_lines.set_segments(airfoil_segments(nodes, elements))
            title.set_text(
                rf"$\theta={float(row['theta']):.3f}$"
                + "   "
                + rf"$q_{{\min}}={float(row['min_quality']):.3f}$"
            )
            writer.grab_frame()
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path, nargs="?", default=Path("output"))
    parser.add_argument(
        "--gif-dpi",
        type=int,
        default=150,
        help="GIF resolution in dots per inch (default: 150)",
    )
    arguments = parser.parse_args()
    if arguments.gif_dpi < 1:
        parser.error("--gif-dpi must be positive")
    output = arguments.output.resolve()
    history = rows(output / "continuation.csv")
    snapshot_dir = output / "snapshots"
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 13,
            "legend.fontsize": 11,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    selected = selected_rows(history)
    plot_data_path(figure_dir, snapshot_dir, selected)
    plot_mesh_path(figure_dir, snapshot_dir, selected)
    plot_quality(figure_dir, history)
    plot_boundary_quality_comparison(figure_dir, snapshot_dir, history[-1])
    animate_smoothing_mechanism(
        figure_dir,
        snapshot_dir,
        history,
        arguments.gif_dpi,
    )
    animate_mesh(
        figure_dir,
        snapshot_dir,
        history,
        arguments.gif_dpi,
    )
    print(f"figures: {figure_dir}")


if __name__ == "__main__":
    main()
