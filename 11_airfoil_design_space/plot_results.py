#!/usr/bin/env python3

"""Plot the accepted design set, validator behavior, and selected meshes."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

_cache = Path(tempfile.gettempdir()) / "afepack-design-space-plot-cache"
_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from airfoil_parameterization import (
    ENGINEERING_BOUNDS,
    PARAMETER_NAMES,
    AirfoilGeometry,
    read_uiuc,
    validate_airfoil,
)


LABELS = {
    "m": r"maximum camber $m/c$",
    "x_c": r"camber station $x_c/c$",
    "t": r"maximum thickness $t/c$",
    "x_t": r"thickness station $x_t/c$",
    "t_te": r"trailing-edge thickness $t_{TE}/c$",
}


def read_rows(filename: Path) -> list[dict[str, str]]:
    with filename.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def row_values(row: dict[str, str]) -> np.ndarray:
    return np.array([float(row[name]) for name in PARAMETER_NAMES])


def draw_outline(axis: plt.Axes, geometry: AirfoilGeometry, **kwargs) -> None:
    axis.plot(geometry.upper[:, 0], geometry.upper[:, 1], **kwargs)
    axis.plot(geometry.lower[:, 0], geometry.lower[:, 1], **kwargs)
    axis.plot(
        [geometry.upper[-1, 0], geometry.lower[-1, 0]],
        [geometry.upper[-1, 1], geometry.lower[-1, 1]],
        **kwargs,
    )
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(-0.03, 1.03)
    axis.set_ylim(-0.30, 0.30)
    axis.axis("off")


def plot_design_space(output: Path, figure_directory: Path) -> None:
    accepted = read_rows(output / "accepted.csv")
    selected = read_rows(output / "selected.csv")
    accepted_values = np.vstack([row_values(row) for row in accepted])
    selected_values = np.vstack([row_values(row) for row in selected])
    colors = plt.cm.tab10(np.arange(len(selected)) % 10)

    figure, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))
    axes[0].scatter(
        accepted_values[:, 0],
        accepted_values[:, 2],
        s=12,
        color="#9ca3af",
        alpha=0.42,
        linewidths=0,
        label="valid candidates",
    )
    axes[0].scatter(
        selected_values[:, 0],
        selected_values[:, 2],
        s=72,
        c=colors,
        edgecolor="#111827",
        linewidth=0.7,
        label="maximin mesh set",
        zorder=3,
    )
    for index, point in enumerate(selected_values):
        axes[0].annotate(str(index + 1), (point[0], point[2]), xytext=(5, 4),
                         textcoords="offset points", fontsize=9)
    axes[0].set_xlabel(LABELS["m"])
    axes[0].set_ylabel(LABELS["t"])
    axes[0].spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, loc="best")
    axes[0].set_title("Wide proposals, deterministic geometric filter")

    normalized = (selected_values - ENGINEERING_BOUNDS[:, 0]) / (
        ENGINEERING_BOUNDS[:, 1] - ENGINEERING_BOUNDS[:, 0]
    )
    locations = np.arange(len(PARAMETER_NAMES))
    for index, values in enumerate(normalized):
        axes[1].plot(locations, values, marker="o", color=colors[index],
                     linewidth=2.0, label=str(index + 1))
    axes[1].set_xticks(locations, [r"$m$", r"$x_c$", r"$t$", r"$x_t$", r"$t_{TE}$"])
    axes[1].set_ylim(-0.04, 1.04)
    axes[1].set_ylabel("normalized engineering range")
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].set_title("Selected points remain separated in 5-D")
    axes[1].legend(title="sample", ncol=2, frameon=False, loc="best")
    figure.tight_layout()
    figure.savefig(figure_directory / "design_space_selection.png", dpi=190)
    plt.close(figure)


def plot_selected_shapes(output: Path, figure_directory: Path) -> None:
    selected_ids = (output / "selected_ids.txt").read_text().split()
    columns = 3
    rows = int(np.ceil(len(selected_ids) / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(13.2, 2.9 * rows))
    axes = np.atleast_1d(axes).ravel()
    for index, (axis, sample_id) in enumerate(zip(axes, selected_ids)):
        geometry = read_uiuc(output / "data" / f"{sample_id}.dat")
        draw_outline(axis, geometry, color="#17324d", linewidth=2.3)
        axis.set_title(f"{index + 1}. {sample_id}", fontsize=11)
    for axis in axes[len(selected_ids):]:
        axis.axis("off")
    figure.tight_layout()
    figure.savefig(figure_directory / "selected_airfoils.png", dpi=190)
    plt.close(figure)


def read_easymesh(case_directory: Path) -> tuple[np.ndarray, np.ndarray]:
    node_file = (case_directory / "airfoil.n").read_text().splitlines()
    element_file = (case_directory / "airfoil.e").read_text().splitlines()
    node_count = int(node_file[0].split()[0])
    element_count = int(element_file[0].split()[0])
    # EasyMesh appends a dashed separator and a column legend.  Read the
    # declared number of records instead of relying on end-of-file.
    node_lines = node_file[1:1 + node_count]
    element_lines = element_file[1:1 + element_count]
    nodes = np.array(
        [[float(line.split()[1]), float(line.split()[2])] for line in node_lines if line.strip()]
    )
    elements = np.array(
        [[int(line.split()[1]), int(line.split()[2]), int(line.split()[3])]
         for line in element_lines if line.strip()],
        dtype=int,
    )
    return nodes, elements


def triangle_quality(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    triangles = nodes[elements]
    first = triangles[:, 1] - triangles[:, 0]
    second = triangles[:, 2] - triangles[:, 0]
    doubled_area = np.abs(first[:, 0] * second[:, 1] - first[:, 1] * second[:, 0])
    edge_square_sum = np.sum(
        (triangles[:, 1] - triangles[:, 0]) ** 2
        + (triangles[:, 2] - triangles[:, 1]) ** 2
        + (triangles[:, 0] - triangles[:, 2]) ** 2,
        axis=1,
    )
    return 2.0 * np.sqrt(3.0) * doubled_area / edge_square_sum


def plot_selected_meshes(output: Path, figure_directory: Path) -> None:
    selected_ids = (output / "selected_ids.txt").read_text().split()
    columns = 3
    rows = int(np.ceil(len(selected_ids) / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(13.2, 3.25 * rows))
    axes = np.atleast_1d(axes).ravel()
    quality_rows: list[tuple[str, int, int, float]] = []
    for index, (axis, sample_id) in enumerate(zip(axes, selected_ids)):
        nodes, elements = read_easymesh(output / "cases" / sample_id)
        quality = triangle_quality(nodes, elements)
        quality_rows.append(
            (sample_id, len(nodes), len(elements), float(np.min(quality)))
        )
        edge_indices: set[tuple[int, int]] = set()
        for triangle in elements:
            for first, second in ((0, 1), (1, 2), (2, 0)):
                edge_indices.add(tuple(sorted((triangle[first], triangle[second]))))
        segments = [[nodes[first], nodes[second]] for first, second in edge_indices]
        axis.add_collection(
            LineCollection(segments, colors="#2f6b9a", linewidths=0.48)
        )
        axis.set_xlim(-0.16, 1.16)
        axis.set_ylim(-0.35, 0.35)
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        axis.set_title(
            f"{index + 1}. {sample_id}  |  {len(nodes)} nodes, "
            f"{len(elements)} cells, $q_{{min}}={np.min(quality):.2f}$",
            fontsize=9.5,
        )
    for axis in axes[len(selected_ids):]:
        axis.axis("off")
    figure.tight_layout()
    figure.savefig(figure_directory / "selected_meshes.png", dpi=210)
    plt.close(figure)
    with (output / "mesh_quality.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("sample_id", "nodes", "elements", "minimum_quality"))
        writer.writerows(quality_rows)


def plot_validator_examples(output: Path, figure_directory: Path) -> None:
    sample_id = (output / "selected_ids.txt").read_text().split()[0]
    valid = read_uiuc(output / "data" / f"{sample_id}.dat")
    reversed_data = AirfoilGeometry(valid.upper[::-1].copy(), valid.lower.copy())
    tangled = AirfoilGeometry(valid.upper.copy(), valid.lower.copy())
    middle = len(tangled.upper) // 2
    tangled.upper[middle - 2:middle + 3, 1] = (
        tangled.lower[middle - 2:middle + 3, 1] - 0.025
    )
    examples = [
        ("accepted", valid),
        ("reversed point order", reversed_data),
        ("crossed / tangled surfaces", tangled),
    ]
    figure, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
    for axis, (label, geometry) in zip(axes, examples):
        result = validate_airfoil(geometry)
        draw_outline(
            axis,
            geometry,
            color="#2a9d8f" if result.valid else "#d95d4f",
            linewidth=2.2,
        )
        status = "PASS" if result.valid else "REJECT"
        reason = "" if result.valid else "\n" + ", ".join(result.reasons[:2])
        axis.set_title(f"{status}: {label}{reason}", fontsize=10)
    figure.tight_layout()
    figure.savefig(figure_directory / "validator_examples.png", dpi=190)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    figure_directory = output / "figures"
    figure_directory.mkdir(parents=True, exist_ok=True)
    plot_design_space(output, figure_directory)
    plot_selected_shapes(output, figure_directory)
    plot_selected_meshes(output, figure_directory)
    plot_validator_examples(output, figure_directory)
    print(f"Figures: {figure_directory}")


if __name__ == "__main__":
    main()
