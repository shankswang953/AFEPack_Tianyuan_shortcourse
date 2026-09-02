#!/usr/bin/env python3

"""Sample, validate, and select a dispersed low-dimensional airfoil design."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from airfoil_parameterization import (
    ENGINEERING_BOUNDS,
    PARAMETER_NAMES,
    PROPOSAL_BOUNDS,
    AirfoilParameters,
    generate_airfoil,
    validate_airfoil,
    write_uiuc,
)


def farthest_point_selection(points: np.ndarray, count: int) -> list[int]:
    if count > len(points):
        raise ValueError("not enough valid candidates for requested selection")
    center = np.full(points.shape[1], 0.5)
    selected = [int(np.argmin(np.linalg.norm(points - center, axis=1)))]
    nearest = np.linalg.norm(points - points[selected[0]], axis=1)
    while len(selected) < count:
        nearest[selected] = -1.0
        next_index = int(np.argmax(nearest))
        selected.append(next_index)
        nearest = np.minimum(
            nearest,
            np.linalg.norm(points - points[next_index], axis=1),
        )
    return selected


def write_table(filename: Path, rows: list[dict[str, object]]) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "candidate_id",
        *PARAMETER_NAMES,
        "valid",
        "reasons",
        "signed_area",
        "minimum_gap",
        "maximum_curvature",
    ]
    with filename.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("output"))
    parser.add_argument("--candidates", type=int, default=2048)
    parser.add_argument("--selected", type=int, default=6)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    sampler = qmc.LatinHypercube(d=len(PARAMETER_NAMES), seed=arguments.seed)
    unit_samples = sampler.random(arguments.candidates)
    physical = qmc.scale(
        unit_samples,
        PROPOSAL_BOUNDS[:, 0],
        PROPOSAL_BOUNDS[:, 1],
    )

    all_rows: list[dict[str, object]] = []
    valid_rows: list[dict[str, object]] = []
    valid_units: list[np.ndarray] = []
    valid_geometries = []
    reason_counts: Counter[str] = Counter()

    engineering_span = ENGINEERING_BOUNDS[:, 1] - ENGINEERING_BOUNDS[:, 0]
    for candidate_id, values in enumerate(physical):
        parameters = AirfoilParameters.from_array(values)
        geometry = generate_airfoil(parameters)
        result = validate_airfoil(geometry, parameters)
        row: dict[str, object] = {
            "candidate_id": candidate_id,
            **dict(zip(PARAMETER_NAMES, values)),
            "valid": int(result.valid),
            "reasons": ";".join(result.reasons),
            "signed_area": result.signed_area,
            "minimum_gap": result.minimum_gap,
            "maximum_curvature": result.maximum_curvature,
        }
        all_rows.append(row)
        if result.valid:
            valid_rows.append(row)
            valid_units.append((values - ENGINEERING_BOUNDS[:, 0]) / engineering_span)
            valid_geometries.append(geometry)
        else:
            reason_counts.update(result.reasons)

    if len(valid_rows) < arguments.selected:
        raise RuntimeError(
            f"only {len(valid_rows)} valid candidates; "
            f"cannot select {arguments.selected}"
        )

    selected_local = farthest_point_selection(
        np.asarray(valid_units), arguments.selected
    )
    selected_rows = [valid_rows[index] for index in selected_local]

    write_table(output / "candidates.csv", all_rows)
    write_table(output / "accepted.csv", valid_rows)
    write_table(
        output / "rejected.csv",
        [row for row in all_rows if not row["valid"]],
    )
    write_table(output / "selected.csv", selected_rows)

    data_directory = output / "data"
    selected_ids: list[str] = []
    for local_index, row in zip(selected_local, selected_rows):
        candidate_id = int(row["candidate_id"])
        sample_name = f"airfoil_{candidate_id:04d}"
        selected_ids.append(sample_name)
        write_uiuc(
            data_directory / f"{sample_name}.dat",
            valid_geometries[local_index],
            f"DESIGN SAMPLE {candidate_id:04d}",
        )
    (output / "selected_ids.txt").write_text(
        "".join(f"{name}\n" for name in selected_ids), encoding="utf-8"
    )

    with (output / "rejection_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(("reason", "count"))
        writer.writerows(reason_counts.most_common())

    print(f"Candidates proposed : {len(all_rows)}")
    print(f"Candidates accepted : {len(valid_rows)}")
    print(f"Candidates rejected : {len(all_rows) - len(valid_rows)}")
    print(f"Dispersed selections: {len(selected_rows)}")
    for name, row in zip(selected_ids, selected_rows):
        values = " ".join(f"{key}={float(row[key]):.4f}" for key in PARAMETER_NAMES)
        print(f"  {name}: {values}")


if __name__ == "__main__":
    main()
