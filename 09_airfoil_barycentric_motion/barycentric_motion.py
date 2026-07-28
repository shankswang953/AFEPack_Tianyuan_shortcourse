#!/usr/bin/env python3

"""Move one AFEPack mesh along a fixed-topology barycentric airfoil path.

Command-line usage:
    python3 barycentric_motion.py [--steps N] [--smooth-iterations N]
        [--boundary-quality-iterations N] [--relaxation FLOAT]
        [--quality-floor FLOAT] [--max-halvings N] [--reset] [--output DIR]

Inputs are `data/initial_circle.dat` and `target_naca0012.dat`. The default
output is `output/`, containing continuation CSV/JSON data, persistent and
trial meshes, snapshots, and per-step `dat` files. Configure EASYMESH_BIN and
EASYMESH2MESH_BIN when their `$HOME/bin` defaults are not valid.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Airfoil:
    title: str
    upper: list[tuple[float, float]]
    lower: list[tuple[float, float]]


def read_airfoil(filename: Path) -> Airfoil:
    with filename.open() as stream:
        title = stream.readline().strip()
        counts = stream.readline().split()
        if len(counts) < 2:
            raise ValueError(f"invalid point-count line in {filename}")
        upper_count = int(round(float(counts[0].rstrip("."))))
        lower_count = int(round(float(counts[1].rstrip("."))))
        points = [
            (float(fields[0]), float(fields[1]))
            for line in stream
            if len(fields := line.split()) >= 2
        ]
    if len(points) != upper_count + lower_count:
        raise ValueError(f"point count does not match header in {filename}")
    return Airfoil(
        title=title,
        upper=points[:upper_count],
        lower=points[upper_count:],
    )


def validate_pair(initial: Airfoil, target: Airfoil) -> None:
    for name in ("upper", "lower"):
        initial_points = getattr(initial, name)
        target_points = getattr(target, name)
        if len(initial_points) != len(target_points):
            raise ValueError(f"{name} point counts differ")
        for index, (first, second) in enumerate(
            zip(initial_points, target_points)
        ):
            if abs(first[0] - second[0]) > 1.0e-12:
                raise ValueError(
                    f"{name} x coordinates differ at index {index}"
                )


def barycenter(initial: Airfoil, target: Airfoil, theta: float) -> Airfoil:
    """Euclidean barycenter on the common upper/lower x grid."""

    validate_pair(initial, target)

    def interpolate(
        first: list[tuple[float, float]],
        second: list[tuple[float, float]],
    ) -> list[tuple[float, float]]:
        return [
            (x0, (1.0 - theta) * y0 + theta * y1)
            for (x0, y0), (x1, y1) in zip(first, second)
            if abs(x0 - x1) <= 1.0e-12
        ]

    return Airfoil(
        title=f"BARYCENTRIC AIRFOIL theta={theta:.10f}",
        upper=interpolate(initial.upper, target.upper),
        lower=interpolate(initial.lower, target.lower),
    )


def write_airfoil(filename: Path, airfoil: Airfoil) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w") as stream:
        stream.write(f" {airfoil.title}\n")
        stream.write(
            f"      {len(airfoil.upper)}.       "
            f"{len(airfoil.lower)}.\n\n"
        )
        for x, y in airfoil.upper:
            stream.write(f" {x:.10f} {y:.10f}\n")
        stream.write("\n")
        for x, y in airfoil.lower:
            stream.write(f" {x:.10f} {y:.10f}\n")


def run(
    command: list[str],
    *,
    cwd: Path | None = None,
    allow_failure: bool = False,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0 and not allow_failure:
        raise RuntimeError(
            f"command failed ({result.returncode}): "
            f"{' '.join(command)}\n{result.stdout}"
        )
    return result


def read_quality(filename: Path) -> dict[str, float | int]:
    with filename.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["stage"] == "smoothed":
                return {
                    "min_area": float(row["min_area"]),
                    "min_quality": float(row["min_shape_quality"]),
                    "inverted": int(row["inverted_elements"]),
                }
    raise RuntimeError(f"smoothed row is missing from {filename}")


def read_connectivity(filename: Path) -> list[tuple[int, int, int]]:
    with filename.open(newline="") as stream:
        return [
            (int(row["v0"]), int(row["v1"]), int(row["v2"]))
            for row in csv.DictReader(stream)
        ]


def connectivity_digest(
    connectivity: list[tuple[int, int, int]],
) -> str:
    digest = hashlib.sha256()
    for triangle in connectivity:
        digest.update(
            f"{triangle[0]},{triangle[1]},{triangle[2]}\n".encode()
        )
    return digest.hexdigest()


def csv_row_count(filename: Path) -> int:
    with filename.open() as stream:
        return max(0, sum(1 for _ in stream) - 1)


def copy_snapshot(
    source: Path,
    destination: Path,
    index: int,
    data_file: Path,
) -> tuple[int, int, str]:
    prefix = destination / f"stage_{index:03d}"
    for suffix, source_name in (
        (".mesh", "mesh_smoothed.mesh"),
        ("_nodes.csv", "mesh_smoothed_nodes.csv"),
        ("_elements.csv", "mesh_smoothed_elements.csv"),
        ("_before.mesh", "mesh_moved_unsmoothed.mesh"),
        ("_before_nodes.csv", "mesh_moved_unsmoothed_nodes.csv"),
        ("_before_elements.csv", "mesh_moved_unsmoothed_elements.csv"),
        ("_quality.csv", "quality_summary.csv"),
        (".dat", None),
    ):
        source_file = data_file if source_name is None else source / source_name
        shutil.copy2(source_file, Path(str(prefix) + suffix))
    nodes = csv_row_count(Path(str(prefix) + "_nodes.csv"))
    elements_file = Path(str(prefix) + "_elements.csv")
    connectivity = read_connectivity(elements_file)
    return nodes, len(connectivity), connectivity_digest(connectivity)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=48)
    parser.add_argument(
        "--smooth-iterations",
        type=int,
        default=0,
        help="unrestricted Laplacian sweeps before smart smoothing",
    )
    parser.add_argument(
        "--boundary-quality-iterations",
        type=int,
        default=20,
        help="equilateral-target and quality-aware smoothing sweeps",
    )
    parser.add_argument("--relaxation", type=float, default=0.45)
    parser.add_argument("--quality-floor", type=float, default=0.40)
    parser.add_argument("--max-halvings", type=int, default=8)
    parser.add_argument("--reset", action="store_true")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    if arguments.steps < 2:
        raise SystemExit("--steps must be at least two")
    if arguments.smooth_iterations < 0:
        raise SystemExit("--smooth-iterations must be nonnegative")
    if arguments.boundary_quality_iterations < 0:
        raise SystemExit("--boundary-quality-iterations must be nonnegative")
    if not 0.0 < arguments.relaxation <= 1.0:
        raise SystemExit("--relaxation must lie in (0, 1]")
    if arguments.quality_floor < 0.0:
        raise SystemExit("--quality-floor must be nonnegative")
    if arguments.max_halvings < 0:
        raise SystemExit("--max-halvings must be nonnegative")

    root = Path(__file__).resolve().parent
    backend = root / "backend"
    output = (
        arguments.output.resolve()
        if arguments.output
        else root / "output"
    )
    if arguments.reset and output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)

    initial_file = root / "data" / "initial_circle.dat"
    target_file = root / "data" / "target_naca0012.dat"
    initial = read_airfoil(initial_file)
    target = read_airfoil(target_file)
    validate_pair(initial, target)

    generate = backend / "generate_airfoil_geometry"
    move_and_smooth = backend / "move_and_smooth"
    easymesh = Path(
        os.environ.get("EASYMESH_BIN", str(Path.home() / "bin/easymesh"))
    )
    converter = Path(
        os.environ.get(
            "EASYMESH2MESH_BIN",
            str(Path.home() / "bin/easymesh2mesh"),
        )
    )
    run(["make", "-C", str(backend)])

    setup = output / "initial_setup"
    persistent = output / "persistent"
    trial = output / "trial"
    snapshots = output / "snapshots"
    data_steps = output / "data_steps"
    for directory in (setup, persistent, snapshots, data_steps):
        directory.mkdir(parents=True, exist_ok=True)

    step_zero_data = data_steps / "theta_0000000000.dat"
    write_airfoil(step_zero_data, barycenter(initial, target, 0.0))
    run(
        [
            str(generate),
            str(step_zero_data),
            str(setup),
            "96",
            "0.35",
            "0.12",
            "0.0",
            "0.0",
            str(step_zero_data),
        ]
    )
    run([str(easymesh), "airfoil.d"], cwd=setup, allow_failure=True)
    if not all((setup / f"airfoil.{suffix}").exists() for suffix in "nse"):
        raise RuntimeError("EasyMesh did not create airfoil.[nse]")
    run([str(converter), "airfoil", "airfoil.mesh"], cwd=setup)
    run(
        [
            str(move_and_smooth),
            str(setup / "airfoil.mesh"),
            str(setup / "boundary_initial.dat"),
            str(setup / "boundary_initial.dat"),
            str(setup),
            str(arguments.smooth_iterations),
            str(arguments.relaxation),
            str(arguments.boundary_quality_iterations),
        ]
    )

    current_mesh = persistent / "mesh_current.mesh"
    current_boundary = persistent / "boundary_current.dat"
    current_data = persistent / "data_current.dat"
    shutil.copy2(setup / "mesh_smoothed.mesh", current_mesh)
    shutil.copy2(setup / "boundary_moved.dat", current_boundary)
    shutil.copy2(step_zero_data, current_data)

    initial_quality = read_quality(setup / "quality_summary.csv")
    nodes, elements, topology_hash = copy_snapshot(
        setup,
        snapshots,
        0,
        step_zero_data,
    )
    history = [
        {
            "stage": 0,
            "theta": 0.0,
            "delta_theta": 0.0,
            **initial_quality,
            "nodes": nodes,
            "elements": elements,
            "connectivity_sha256": topology_hash,
        }
    ]

    nominal_delta = 1.0 / arguments.steps
    minimum_delta = nominal_delta / (2**arguments.max_halvings)
    delta = nominal_delta
    theta = 0.0
    stage = 0
    rejected_trials = 0

    while theta < 1.0 - 1.0e-13:
        next_theta = min(1.0, theta + delta)
        candidate_data = data_steps / (
            f"theta_{round(next_theta * 1.0e10):010d}.dat"
        )
        write_airfoil(
            candidate_data,
            barycenter(initial, target, next_theta),
        )
        if trial.exists():
            shutil.rmtree(trial)
        trial.mkdir(parents=True)
        run(
            [
                str(generate),
                str(current_data),
                str(trial),
                "96",
                "0.35",
                "0.12",
                "0.0",
                "0.0",
                str(candidate_data),
            ]
        )
        result = run(
            [
                str(move_and_smooth),
                str(current_mesh),
                str(current_boundary),
                str(trial / "boundary_moved.dat"),
                str(trial),
                str(arguments.smooth_iterations),
                str(arguments.relaxation),
                str(arguments.boundary_quality_iterations),
            ],
            allow_failure=True,
        )
        quality = (
            read_quality(trial / "quality_summary.csv")
            if (trial / "quality_summary.csv").exists()
            else None
        )
        valid = (
            result.returncode == 0
            and quality is not None
            and quality["inverted"] == 0
            and quality["min_quality"] >= arguments.quality_floor
        )
        if not valid:
            rejected_trials += 1
            delta *= 0.5
            if delta < minimum_delta - 1.0e-15:
                detail = (
                    result.stdout.strip()
                    if result.stdout.strip()
                    else f"quality={quality}"
                )
                raise RuntimeError(
                    "fixed-topology continuation could not satisfy the "
                    f"quality floor near theta={next_theta:.8f}\n{detail}"
                )
            print(
                f"reject theta={next_theta:.6f}; "
                f"retry with delta_theta={delta:.6f}"
            )
            continue

        connectivity = read_connectivity(
            trial / "mesh_smoothed_elements.csv"
        )
        if connectivity_digest(connectivity) != topology_hash:
            raise RuntimeError("element connectivity changed")

        accepted_delta = next_theta - theta
        theta = next_theta
        stage += 1
        shutil.copy2(trial / "mesh_smoothed.mesh", current_mesh)
        shutil.copy2(trial / "boundary_moved.dat", current_boundary)
        shutil.copy2(candidate_data, current_data)
        nodes, elements, current_hash = copy_snapshot(
            trial,
            snapshots,
            stage,
            candidate_data,
        )
        history.append(
            {
                "stage": stage,
                "theta": theta,
                "delta_theta": accepted_delta,
                **quality,
                "nodes": nodes,
                "elements": elements,
                "connectivity_sha256": current_hash,
            }
        )
        print(
            f"accept stage={stage:03d} theta={theta:.6f} "
            f"delta={accepted_delta:.6f} "
            f"quality={quality['min_quality']:.6f}"
        )
        delta = min(nominal_delta, delta * 1.5)

    history_file = output / "continuation.csv"
    with history_file.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(history[0]))
        writer.writeheader()
        writer.writerows(history)

    summary = {
        "path": "D(theta) = (1-theta) D0 + theta D1",
        "easy_mesh_calls": 1,
        "accepted_stages": stage,
        "rejected_trials": rejected_trials,
        "nodes": nodes,
        "elements": elements,
        "connectivity_sha256": topology_hash,
        "connectivity_unchanged": True,
        "minimum_quality_over_path": min(
            float(row["min_quality"]) for row in history
        ),
        "unrestricted_laplacian_iterations_per_stage":
            arguments.smooth_iterations,
        "boundary_quality_iterations_per_stage":
            arguments.boundary_quality_iterations,
        "relaxation": arguments.relaxation,
        "quality_floor": arguments.quality_floor,
        "final_theta": theta,
    }
    with (output / "summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    print("\nBarycentric fixed-topology continuation completed")
    print(f"  EasyMesh calls:       {summary['easy_mesh_calls']}")
    print(f"  accepted stages:      {summary['accepted_stages']}")
    print(f"  rejected trials:      {summary['rejected_trials']}")
    print(f"  nodes:                {summary['nodes']}")
    print(f"  elements:             {summary['elements']}")
    print(
        "  minimum path quality: "
        f"{summary['minimum_quality_over_path']:.6f}"
    )
    print("  connectivity changed: no")
    print(f"  history:              {history_file}")


if __name__ == "__main__":
    main()
