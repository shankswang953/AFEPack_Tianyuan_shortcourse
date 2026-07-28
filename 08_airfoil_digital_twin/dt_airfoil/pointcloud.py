"""Extract airfoil boundary curves and compute mesh-independent distances.

This import-only module reads caller-supplied mesh node/element CSV paths and
returns NumPy arrays. It has no command-line interface and creates no files.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class BoundaryCurves:
    upper: np.ndarray
    lower: np.ndarray

    @property
    def all_points(self) -> np.ndarray:
        return np.vstack((self.upper, self.lower[1:-1]))


def _read_mesh(
    nodes_csv: Path,
    elements_csv: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    node_rows: list[tuple[int, float, float, int]] = []
    with nodes_csv.open(newline="") as stream:
        for row in csv.DictReader(stream):
            node_rows.append(
                (
                    int(row["index"]),
                    float(row["x"]),
                    float(row["y"]),
                    int(row["boundary_mark"]),
                )
            )
    node_rows.sort()
    points = np.asarray(
        [[x, y] for _, x, y, _ in node_rows],
        dtype=np.float64,
    )
    marks = np.asarray(
        [mark for _, _, _, mark in node_rows],
        dtype=np.int32,
    )

    elements: list[list[int]] = []
    with elements_csv.open(newline="") as stream:
        for row in csv.DictReader(stream):
            elements.append(
                [int(row["v0"]), int(row["v1"]), int(row["v2"])]
            )
    return points, marks, np.asarray(elements, dtype=np.int64)


def extract_boundary_curves(
    nodes_csv: Path,
    elements_csv: Path,
    boundary_mark: int = 3,
) -> BoundaryCurves:
    """Recover upper/lower paths using only mesh vertices and connectivity.

    No current-to-target point correspondence or vertex numbering is used.
    """

    points, marks, elements = _read_mesh(nodes_csv, elements_csv)
    edge_counts: Counter[tuple[int, int]] = Counter()
    for triangle in elements:
        for first, second in (
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ):
            edge_counts[tuple(sorted((int(first), int(second))))] += 1

    adjacency: dict[int, list[int]] = defaultdict(list)
    for (first, second), count in edge_counts.items():
        if (
            count == 1
            and marks[first] == boundary_mark
            and marks[second] == boundary_mark
        ):
            adjacency[first].append(second)
            adjacency[second].append(first)

    boundary_vertices = sorted(adjacency)
    if len(boundary_vertices) < 8:
        raise ValueError("could not recover the marked internal boundary")
    bad_degrees = {
        vertex: len(adjacency[vertex])
        for vertex in boundary_vertices
        if len(adjacency[vertex]) != 2
    }
    if bad_degrees:
        raise ValueError(
            f"internal boundary is not a simple closed loop: {bad_degrees}"
        )

    leading = min(boundary_vertices, key=lambda index: points[index, 0])
    trailing = max(boundary_vertices, key=lambda index: points[index, 0])

    def walk(first_neighbor: int) -> list[int]:
        path = [leading, first_neighbor]
        previous = leading
        current = first_neighbor
        while current != trailing:
            candidates = [
                vertex
                for vertex in adjacency[current]
                if vertex != previous
            ]
            if len(candidates) != 1:
                raise ValueError("ambiguous boundary traversal")
            previous, current = current, candidates[0]
            if current == leading or len(path) > len(boundary_vertices) + 1:
                raise ValueError("boundary path did not reach trailing edge")
            path.append(current)
        return path

    first_path = walk(adjacency[leading][0])
    second_path = walk(adjacency[leading][1])
    first_points = points[np.asarray(first_path)]
    second_points = points[np.asarray(second_path)]
    if float(np.mean(first_points[:, 1])) >= float(
        np.mean(second_points[:, 1])
    ):
        upper, lower = first_points, second_points
    else:
        upper, lower = second_points, first_points
    return BoundaryCurves(upper=upper, lower=lower)
