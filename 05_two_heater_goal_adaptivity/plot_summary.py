#!/usr/bin/env python3
"""Create compact PNG summaries from the numerical comparison table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import PowerNorm


COLORS = {
    "residual": "#d97706",
    "dual": "#0f766e",
    "dwr": "#6d42c2",
}
LABELS = {
    "residual": "residual",
    "dual": "dual magnitude",
    "dwr": "DWR",
}


def read_comparison(path: Path) -> tuple[float, dict[str, np.ndarray]]:
    reference = None
    rows: dict[str, list[list[float]]] = {name: [] for name in COLORS}
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# spectral_reference "):
                reference = float(line.split()[2])
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            rows[fields[0]].append([float(value) for value in fields[1:]])
    if reference is None:
        raise ValueError(f"missing spectral reference in {path}")
    return reference, {
        name: np.asarray(values, dtype=float) for name, values in rows.items()
    }


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 12,
            "axes.labelsize": 12,
            "legend.fontsize": 11,
            "lines.linewidth": 2.4,
            "lines.markersize": 7,
        }
    )


def read_dx_mesh(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    point_count = int(re.search(r"item (\d+)", lines[0]).group(1))
    points = np.asarray(
        [[float(value) for value in line.split()] for line in lines[1 : point_count + 1]],
        dtype=float,
    )
    connection_header = point_count + 1
    while not lines[connection_header].startswith("object 2 "):
        connection_header += 1
    element_count = int(
        re.search(r"item (\d+)", lines[connection_header]).group(1)
    )
    triangles = np.asarray(
        [
            [int(value) for value in line.split()]
            for line in lines[
                connection_header + 1 : connection_header + 1 + element_count
            ]
        ],
        dtype=int,
    )
    return points, triangles


def plot_meshes(output_dir: Path, path: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for axis, name in zip(axes, ("residual", "dual", "dwr")):
        points, triangles = read_dx_mesh(output_dir / "meshes" / name / "final.dx")
        axis.triplot(
            points[:, 0],
            points[:, 1],
            triangles,
            color=COLORS[name],
            linewidth=0.38,
        )
        axis.set_title(LABELS[name])
        axis.set_aspect("equal")
        axis.set_axis_off()
    figure.savefig(path, dpi=200, facecolor="white")
    plt.close(figure)


def latest_round(directory: Path, stem: str) -> int:
    rounds = [
        int(match.group(1))
        for path in directory.glob(f"{stem}_round_*.dat")
        if (match := re.search(r"_round_(\d+)\.dat$", path.name))
    ]
    if not rounds:
        raise FileNotFoundError(f"no {stem} round data in {directory}")
    return max(rounds)


def indicator_values(name: str, data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if name == "residual":
        # New files include the mesh-size-independent strong residual before
        # the cell contribution; accept old tables as well.
        value_column = 5 if data.shape[1] >= 9 else 4
    elif name == "dual":
        value_column = 3
    else:
        value_column = 4
    return data[:, value_column], data[:, -1].astype(bool)


def plot_indicators(output_dir: Path, round_selector: str, path: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    palettes = {"residual": "Oranges", "dual": "BuGn", "dwr": "Purples"}
    gammas = {"residual": 0.55, "dual": 4.0, "dwr": 0.55}
    for axis, name in zip(axes, ("residual", "dual", "dwr")):
        directory = output_dir / "fields" / name
        round_number = 1 if round_selector == "initial" else latest_round(
            directory, "indicator" if name != "dual" else "magnitude"
        )
        stem = "magnitude" if name == "dual" else "indicator"
        data = np.loadtxt(directory / f"{stem}_round_{round_number}.dat")
        values, marked = indicator_values(name, data)
        points, triangles = read_dx_mesh(
            output_dir / "meshes" / name / f"level_{round_number - 1}.dx"
        )
        maximum = max(float(np.max(values)), np.finfo(float).tiny)
        axis.tripcolor(
            points[:, 0],
            points[:, 1],
            triangles,
            facecolors=values,
            cmap=palettes[name],
            norm=PowerNorm(gamma=gammas[name], vmin=0.0, vmax=maximum),
            edgecolors="#8090a0",
            linewidth=0.20,
        )
        if np.any(marked):
            axis.add_collection(
                PolyCollection(
                    points[triangles[marked]],
                    facecolors="none",
                    edgecolors="#111827",
                    linewidths=0.75,
                )
            )
        axis.set_title(LABELS[name])
        axis.set_aspect("equal")
        axis.set_axis_off()
    figure.savefig(path, dpi=200, facecolor="white")
    plt.close(figure)


def plot_problem_fields(output_dir: Path, path: Path) -> None:
    points, triangles = read_dx_mesh(output_dir / "meshes/residual/level_0.dx")
    fields = (
        ("source_profile.dat", "heat source", "Oranges", 0.55),
        ("sensor_weight.dat", "sensor weight", "BuGn", 4.0),
    )
    figure, axes = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)
    for axis, (filename, title, palette, gamma) in zip(axes, fields):
        values = np.loadtxt(output_dir / "fields/problem" / filename)[:, 3]
        axis.tripcolor(
            points[:, 0],
            points[:, 1],
            triangles,
            facecolors=values,
            cmap=palette,
            norm=PowerNorm(gamma=gamma, vmin=0.0, vmax=float(np.max(values))),
            edgecolors="#8090a0",
            linewidth=0.20,
        )
        axis.set_title(title)
        axis.set_aspect("equal")
        axis.set_axis_off()
    figure.savefig(path, dpi=200, facecolor="white")
    plt.close(figure)


def plot_value(reference: float, records: dict[str, np.ndarray], path: Path) -> None:
    figure, axis = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    for name, values in records.items():
        axis.plot(
            values[:, 2],
            values[:, 3],
            marker="o",
            color=COLORS[name],
            label=LABELS[name],
        )
    axis.axhline(reference, color="#27364a", linestyle="--", label="reference")
    axis.set_xlabel("degrees of freedom")
    axis.set_ylabel(r"sensor functional $J(T_h)$")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False)
    figure.savefig(path, dpi=200, facecolor="white")
    plt.close(figure)


def plot_error(records: dict[str, np.ndarray], path: Path) -> None:
    figure, axis = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    for name, values in records.items():
        axis.semilogy(
            values[:, 2],
            values[:, 4],
            marker="o",
            color=COLORS[name],
            label=LABELS[name],
        )
    axis.set_xlabel("degrees of freedom")
    axis.set_ylabel(r"target error $|J(T_h)-J_{ref}|$")
    axis.grid(alpha=0.25, which="both")
    axis.legend(frameon=False)
    figure.savefig(path, dpi=200, facecolor="white")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    comparison = args.output_dir / "summary" / "functional_comparison.dat"
    figures = args.output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    reference, records = read_comparison(comparison)
    configure_plotting()
    plot_value(reference, records, figures / "functional_value_vs_dofs.png")
    plot_error(records, figures / "functional_error_vs_dofs.png")
    plot_problem_fields(args.output_dir, figures / "problem_fields.png")
    plot_indicators(args.output_dir, "initial", figures / "indicators_initial.png")
    plot_indicators(args.output_dir, "final", figures / "indicators_final.png")
    plot_meshes(args.output_dir, figures / "meshes_final.png")


if __name__ == "__main__":
    main()
