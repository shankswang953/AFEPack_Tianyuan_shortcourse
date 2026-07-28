#!/usr/bin/env python3

"""Create teaching figures from the mesh-motion CSV files.

Command-line usage:
    python3 visualize_results.py OUTPUT_DIR [--figure-dir FIGURE_DIR]

Inputs:
    Fit, node, element, displacement, and quality CSV files in OUTPUT_DIR.
Outputs:
    PNG files in FIGURE_DIR, or OUTPUT_DIR/figures/ when not specified.
This plotting helper does not run AFEPack or modify the persistent mesh.
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

# Keep font and Matplotlib caches local even when the home directory is
# mounted read-only.
_cache_root = Path(tempfile.gettempdir()) / "afepack-airfoil-plot-cache"
_cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache_root / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_root))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon


def read_csv(filename: Path) -> list[dict[str, str]]:
    with filename.open(newline="") as stream:
        return list(csv.DictReader(stream))


def read_nodes(filename: Path) -> list[tuple[float, float, int]]:
    return [
        (float(row["x"]), float(row["y"]), int(row["boundary_mark"]))
        for row in read_csv(filename)
    ]


def read_elements(filename: Path) -> list[tuple[int, int, int]]:
    return [
        (int(row["v0"]), int(row["v1"]), int(row["v2"]))
        for row in read_csv(filename)
    ]


def signed_twice_area(
    nodes: list[tuple[float, float, int]],
    element: tuple[int, int, int],
) -> float:
    a, b, c = (nodes[index] for index in element)
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def configure_axes(axis: plt.Axes) -> None:
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"$x/c$")
    axis.set_ylabel(r"$y/c$")
    axis.tick_params(direction="out", length=3)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def fit_points(
    fit_rows: list[dict[str, str]],
    kind: str,
    surface: str | None = None,
) -> list[tuple[float, float]]:
    return [
        (float(row["x"]), float(row["y"]))
        for row in fit_rows
        if row["kind"] == kind
        and (surface is None or row["surface"] == surface)
    ]


def draw_fit(axis: plt.Axes, fit_rows: list[dict[str, str]]) -> None:
    for surface in ("upper", "lower"):
        raw = fit_points(fit_rows, "raw", surface)
        curve = fit_points(fit_rows, "curve", surface)
        controls = fit_points(fit_rows, "control", surface)
        axis.plot(
            [p[0] for p in raw],
            [p[1] for p in raw],
            linestyle="none",
            marker="o",
            markersize=2.7,
            color="#4b5563",
            label="input data" if surface == "upper" else None,
            zorder=3,
        )
        axis.plot(
            [p[0] for p in controls],
            [p[1] for p in controls],
            marker="s",
            markersize=3.0,
            linewidth=0.75,
            color="#0f766e",
            alpha=0.75,
            label="control polygon" if surface == "upper" else None,
            zorder=2,
        )
        axis.plot(
            [p[0] for p in curve],
            [p[1] for p in curve],
            linewidth=2.0,
            color="#ea580c",
            label="Bezier curve" if surface == "upper" else None,
            zorder=4,
        )

    mesh_points = fit_points(fit_rows, "mesh-point")
    axis.scatter(
        [p[0] for p in mesh_points],
        [p[1] for p in mesh_points],
        s=10,
        facecolors="none",
        edgecolors="#2563eb",
        linewidths=0.7,
        label="EasyMesh points",
        zorder=5,
    )
    configure_axes(axis)
    axis.set_xlim(-0.04, 1.04)
    axis.set_ylim(-0.13, 0.13)
    axis.legend(loc="upper center", ncol=4, frameon=False)


def draw_fit_update(
    axis: plt.Axes,
    initial_rows: list[dict[str, str]],
    moved_rows: list[dict[str, str]],
) -> None:
    for surface in ("upper", "lower"):
        initial_curve = fit_points(initial_rows, "curve", surface)
        moved_curve = fit_points(moved_rows, "curve", surface)
        moved_raw = fit_points(moved_rows, "raw", surface)
        axis.plot(
            [point[0] for point in initial_curve],
            [point[1] for point in initial_curve],
            linewidth=1.5,
            linestyle="--",
            color="#4b5563",
            label="before" if surface == "upper" else None,
            zorder=2,
        )
        axis.plot(
            [point[0] for point in moved_curve],
            [point[1] for point in moved_curve],
            linewidth=2.2,
            color="#ea580c",
            label="after" if surface == "upper" else None,
            zorder=4,
        )
        axis.scatter(
            [point[0] for point in moved_raw],
            [point[1] for point in moved_raw],
            s=8,
            color="#2563eb",
            linewidths=0,
            label="updated data" if surface == "upper" else None,
            zorder=5,
        )

    configure_axes(axis)
    axis.set_xlim(-0.04, 1.04)
    axis.set_ylim(-0.13, 0.13)
    axis.legend(loc="upper center", ncol=3, frameon=False)


def draw_mesh(
    axis: plt.Axes,
    nodes: list[tuple[float, float, int]],
    elements: list[tuple[int, int, int]],
    reference_nodes: list[tuple[float, float, int]],
    show_inverted: bool,
) -> None:
    edges: set[tuple[int, int]] = set()
    for element in elements:
        for first, second in (
            (element[0], element[1]),
            (element[1], element[2]),
            (element[2], element[0]),
        ):
            edges.add((min(first, second), max(first, second)))

    segments = [
        [(nodes[first][0], nodes[first][1]), (nodes[second][0], nodes[second][1])]
        for first, second in edges
    ]
    axis.add_collection(
        LineCollection(segments, colors="#527394", linewidths=0.42, alpha=0.72)
    )

    if show_inverted:
        for element in elements:
            initial_sign = signed_twice_area(reference_nodes, element)
            current_sign = signed_twice_area(nodes, element)
            if initial_sign * current_sign <= 0.0:
                axis.add_patch(
                    Polygon(
                        [(nodes[index][0], nodes[index][1]) for index in element],
                        closed=True,
                        facecolor="#ef4444",
                        edgecolor="#b91c1c",
                        linewidth=0.5,
                        alpha=0.65,
                        zorder=2,
                    )
                )

    airfoil = [(x, y) for x, y, mark in nodes if mark == 3]
    axis.scatter(
        [point[0] for point in airfoil],
        [point[1] for point in airfoil],
        s=5,
        color="#ea580c",
        linewidths=0,
        zorder=3,
    )
    configure_axes(axis)
    axis.set_xlim(-0.18, 1.18)
    axis.set_ylim(-0.32, 0.32)


def save_single(
    filename: Path,
    drawer,
    *,
    width: float = 10.0,
    height: float = 4.8,
) -> None:
    figure, axis = plt.subplots(figsize=(width, height), constrained_layout=True)
    drawer(axis)
    figure.savefig(filename, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--figure-dir",
        type=Path,
        help="default: OUTPUT_DIR/figures",
    )
    arguments = parser.parse_args()

    output_dir = arguments.output_dir.resolve()
    figure_dir = (
        arguments.figure_dir.resolve()
        if arguments.figure_dir
        else output_dir / "figures"
    )
    figure_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    fit_rows = read_csv(output_dir / "fit_initial.csv")
    moved_fit_rows = read_csv(output_dir / "fit_moved.csv")
    elements = read_elements(output_dir / "mesh_initial_elements.csv")
    initial_nodes = read_nodes(output_dir / "mesh_initial_nodes.csv")
    moved_nodes = read_nodes(output_dir / "mesh_moved_unsmoothed_nodes.csv")
    smoothed_nodes = read_nodes(output_dir / "mesh_smoothed_nodes.csv")

    save_single(
        figure_dir / "00_data_update.png",
        lambda axis: draw_fit_update(axis, fit_rows, moved_fit_rows),
    )
    save_single(
        figure_dir / "01_bezier_fit.png",
        lambda axis: draw_fit(axis, fit_rows),
    )
    save_single(
        figure_dir / "02_mesh_initial.png",
        lambda axis: draw_mesh(
            axis, initial_nodes, elements, initial_nodes, False
        ),
    )
    save_single(
        figure_dir / "03_mesh_moved_unsmoothed.png",
        lambda axis: draw_mesh(
            axis, moved_nodes, elements, initial_nodes, True
        ),
    )
    save_single(
        figure_dir / "04_mesh_smoothed.png",
        lambda axis: draw_mesh(
            axis, smoothed_nodes, elements, initial_nodes, False
        ),
    )

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(13.0, 7.2),
        constrained_layout=True,
    )
    inverted_count = sum(
        signed_twice_area(initial_nodes, element)
        * signed_twice_area(moved_nodes, element)
        <= 0.0
        for element in elements
    )

    draw_fit_update(axes[0, 0], fit_rows, moved_fit_rows)
    axes[0, 0].set_title("Persistent data update and Bezier refit")
    draw_mesh(axes[0, 1], initial_nodes, elements, initial_nodes, False)
    axes[0, 1].set_title("Initial mesh")
    draw_mesh(axes[1, 0], moved_nodes, elements, initial_nodes, True)
    axes[1, 0].set_title(
        f"Moved boundary: {inverted_count} inverted elements"
    )
    draw_mesh(axes[1, 1], smoothed_nodes, elements, initial_nodes, False)
    axes[1, 1].set_title("After smoothing: no inverted elements")
    figure.savefig(
        figure_dir / "mesh_motion_overview.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(figure)

    print(f"Figures are in {figure_dir}")


if __name__ == "__main__":
    main()
