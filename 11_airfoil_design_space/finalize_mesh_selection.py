#!/usr/bin/env python3

"""Replace candidates whose EasyMesh output is invalid with ranked reserves."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


def read_mesh(case_directory: Path) -> tuple[np.ndarray, np.ndarray]:
    node_file = (case_directory / "airfoil.n").read_text().splitlines()
    element_file = (case_directory / "airfoil.e").read_text().splitlines()
    node_count = int(node_file[0].split()[0])
    element_count = int(element_file[0].split()[0])
    if element_count > 3 * node_count:
        raise ValueError("implausible_element_count")
    nodes = np.array(
        [
            [float(line.split()[1]), float(line.split()[2])]
            for line in node_file[1:1 + node_count]
        ],
        dtype=float,
    )
    elements = np.array(
        [
            [int(line.split()[1]), int(line.split()[2]), int(line.split()[3])]
            for line in element_file[1:1 + element_count]
        ],
        dtype=int,
    )
    return nodes, elements


def minimum_quality(nodes: np.ndarray, elements: np.ndarray) -> float:
    if len(nodes) == 0 or len(elements) == 0 or not np.all(np.isfinite(nodes)):
        raise ValueError("empty_or_nonfinite_mesh")
    if np.min(elements) < 0 or np.max(elements) >= len(nodes):
        raise ValueError("invalid_vertex_index")
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
    quality = 2.0 * np.sqrt(3.0) * doubled_area / edge_square_sum
    if not np.all(np.isfinite(quality)):
        raise ValueError("nonfinite_triangle_quality")
    return float(np.min(quality))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: finalize_mesh_selection.py OUTPUT_DIRECTORY")
    output = Path(sys.argv[1]).resolve()
    requested = int((output / "selection_request.txt").read_text().strip())
    ranked_ids = (output / "ranked_ids.txt").read_text().split()

    with (output / "ranked.csv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fieldnames = list(reader.fieldnames or [])
        ranked_rows = list(reader)
    row_by_id = {
        f"airfoil_{int(row['candidate_id']):04d}": row for row in ranked_rows
    }

    selected: list[str] = []
    rejected: list[tuple[str, str]] = []
    for sample_id in ranked_ids:
        try:
            nodes, elements = read_mesh(output / "cases" / sample_id)
            quality = minimum_quality(nodes, elements)
            if quality < 0.40:
                raise ValueError(f"minimum_quality={quality:.6f}")
        except (OSError, ValueError, IndexError) as error:
            rejected.append((sample_id, str(error)))
            continue
        if len(selected) < requested:
            selected.append(sample_id)

    if len(selected) < requested:
        raise RuntimeError(
            f"only {len(selected)} meshable candidates; requested {requested}; "
            "increase --reserve"
        )

    (output / "selected_ids.txt").write_text(
        "".join(f"{sample_id}\n" for sample_id in selected), encoding="utf-8"
    )
    with (output / "selected.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row_by_id[sample_id] for sample_id in selected)
    with (output / "mesh_rejections.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(("sample_id", "reason"))
        writer.writerows(rejected)

    print(f"Meshable selections: {len(selected)}")
    for sample_id, reason in rejected:
        print(f"  replaced {sample_id}: {reason}")


if __name__ == "__main__":
    main()
