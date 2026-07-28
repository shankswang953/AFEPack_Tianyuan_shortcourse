#!/usr/bin/env python3

"""Render the Tian/Yuan peer meshes and their common AFEPack mesh.

Command-line usage:
    python3 plot_tianyuan_merge.py OUTPUT_DIR

Inputs:
    OUTPUT_DIR/T_tian.[ne], T_yuan.[ne], T_common.[ne], and
    glyph_distance_field.csv, normally produced by `run.sh`.
Outputs:
    PNG teaching figures below OUTPUT_DIR/figures/.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import ListedColormap


COLORS = {
    "tian": "#2f6b9a",
    "yuan": "#2a9d8f",
    "common": "#d08a24",
    "root": "#6f8193",
}


def read_nodes(filename: Path) -> list[tuple[float, float]]:
    with filename.open() as stream:
        node_count = int(next(stream).split()[0])
        nodes: list[tuple[float, float]] = []
        for _ in range(node_count):
            fields = next(stream).split()
            nodes.append((float(fields[1]), float(fields[2])))
    return nodes


def read_elements(filename: Path) -> list[tuple[int, int, int]]:
    with filename.open() as stream:
        element_count = int(next(stream).split()[0])
        elements: list[tuple[int, int, int]] = []
        for _ in range(element_count):
            fields = next(stream).split()
            elements.append((int(fields[1]), int(fields[2]), int(fields[3])))
    return elements


def read_mesh_segments(
    basename: Path,
) -> tuple[list[list[tuple[float, float]]], int, int]:
    nodes = read_nodes(basename.with_suffix(".n"))
    elements = read_elements(basename.with_suffix(".e"))
    edges: set[tuple[int, int]] = set()
    for a, b, c in elements:
        edges.update(
            {
                (min(a, b), max(a, b)),
                (min(b, c), max(b, c)),
                (min(c, a), max(c, a)),
            }
        )
    return [[nodes[a], nodes[b]] for a, b in edges], len(nodes), len(elements)


def format_axis(axis: plt.Axes) -> None:
    axis.set_xlim(-0.025, 2.425)
    axis.set_ylim(-0.025, 1.025)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)


def add_mesh(axis: plt.Axes, basename: Path, color: str) -> tuple[int, int]:
    segments, node_count, triangle_count = read_mesh_segments(basename)
    axis.add_collection(
        LineCollection(
            segments,
            colors=color,
            linewidths=0.38,
            rasterized=True,
        )
    )
    format_axis(axis)
    return node_count, triangle_count


def save_single_mesh(
    output_dir: Path,
    figure_dir: Path,
    basename: str,
    output_name: str,
    color: str,
) -> None:
    figure, axis = plt.subplots(figsize=(9.2, 3.85))
    add_mesh(axis, output_dir / basename, color)
    figure.subplots_adjust(left=0.005, right=0.995, bottom=0.01, top=0.99)
    figure.savefig(
        figure_dir / output_name,
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.015,
    )
    plt.close(figure)


def save_overview(output_dir: Path, figure_dir: Path) -> None:
    panels = (
        ("T_tian", "Mesh A: TIAN", COLORS["tian"]),
        ("T_yuan", "Mesh B: YUAN", COLORS["yuan"]),
        ("T_common", "Merged mesh: TIANYUAN", COLORS["common"]),
    )
    figure, axes = plt.subplots(3, 1, figsize=(10.8, 10.2))
    for axis, (basename, title, color) in zip(axes, panels):
        node_count, triangle_count = add_mesh(
            axis, output_dir / basename, color
        )
        axis.set_title(title, fontsize=19, fontweight="bold", pad=5)
        axis.text(
            0.99,
            0.02,
            f"{node_count:,} nodes; {triangle_count:,} EasyMesh triangles",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=10.5,
            color="#24313d",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
        )
    figure.subplots_adjust(
        left=0.015,
        right=0.985,
        bottom=0.02,
        top=0.97,
        hspace=0.16,
    )
    figure.savefig(
        figure_dir / "tianyuan_merge_overview.png",
        dpi=280,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(figure)


def read_distance_field(
    filename: Path,
) -> tuple[list[list[tuple[float, float]]], np.ndarray, np.ndarray]:
    polygons: list[list[tuple[float, float]]] = []
    inside_tian: list[float] = []
    inside_yuan: list[float] = []
    with filename.open(newline="") as stream:
        for row in csv.DictReader(stream):
            polygons.append(
                [
                    (float(row["x0"]), float(row["y0"])),
                    (float(row["x1"]), float(row["y1"])),
                    (float(row["x2"]), float(row["y2"])),
                ]
            )
            inside_tian.append(float(row["inside_tian"]))
            inside_yuan.append(float(row["inside_yuan"]))
    return polygons, np.asarray(inside_tian), np.asarray(inside_yuan)


def save_distance_masks(output_dir: Path, figure_dir: Path) -> None:
    polygons, inside_tian, inside_yuan = read_distance_field(
        output_dir / "glyph_distance_field.csv"
    )
    figure, axes = plt.subplots(2, 1, figsize=(10.0, 6.5))
    cmap = ListedColormap(["#f4f1e6", "#d95d4f"])
    for axis, values, title in (
        (axes[0], inside_tian, r"Tian mask: $d(x_K,S_{\rm Tian})=0$"),
        (axes[1], inside_yuan, r"Yuan mask: $d(x_K,S_{\rm Yuan})=0$"),
    ):
        collection = PolyCollection(
            polygons,
            array=values,
            cmap=cmap,
            clim=(0.0, 1.0),
            edgecolors="#7890a6",
            linewidths=0.20,
        )
        axis.add_collection(collection)
        format_axis(axis)
        axis.set_title(title, fontsize=18, fontweight="bold", pad=4)
    figure.subplots_adjust(
        left=0.015,
        right=0.985,
        bottom=0.02,
        top=0.96,
        hspace=0.18,
    )
    figure.savefig(
        figure_dir / "glyph_distance_zero_masks.png",
        dpi=280,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    arguments = parser.parse_args()
    figure_dir = arguments.output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    for basename, output_name, color in (
        ("T_root", "root_mesh.png", COLORS["root"]),
        ("T_tian", "tian_mesh.png", COLORS["tian"]),
        ("T_yuan", "yuan_mesh.png", COLORS["yuan"]),
        ("T_common", "common_mesh.png", COLORS["common"]),
    ):
        save_single_mesh(
            arguments.output_dir,
            figure_dir,
            basename,
            output_name,
            color,
        )
    save_overview(arguments.output_dir, figure_dir)
    save_distance_masks(arguments.output_dir, figure_dir)
    print(f"Figures are in {figure_dir}")


if __name__ == "__main__":
    main()
