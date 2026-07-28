#!/usr/bin/env python3

"""Apply one safe airfoil action or continue smoothing a persistent mesh.

Command-line usage:
    python3 airfoil_step.py CENTER {U,L} SHIFT [OPTIONS]
    python3 airfoil_step.py --smooth-iterations N
    python3 airfoil_step.py --reset

Use `python3 airfoil_step.py --help` for `--width`, shift limits, smoothing,
reset, and display options. The script reads `data/naca0012*.dat`, delegates
mesh work to `run.sh`, and stores persistent state, CSV data, meshes, and
figures below `output/`. Configure external tools with EASYMESH_BIN,
EASYMESH2MESH_BIN, and PLOT_PYTHON when their defaults are not valid.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import platform
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Point:
    x: float
    y: float


@dataclass
class Airfoil:
    title: str
    upper: list[Point]
    lower: list[Point]


def read_airfoil(filename: Path) -> Airfoil:
    with filename.open() as stream:
        title = stream.readline().strip()
        counts = stream.readline().split()
        if len(counts) < 2:
            raise ValueError(f"invalid point-count line in {filename}")
        upper_count = int(round(float(counts[0])))
        lower_count = int(round(float(counts[1].rstrip("."))))
        values: list[Point] = []
        for line in stream:
            fields = line.split()
            if len(fields) >= 2:
                values.append(Point(float(fields[0]), float(fields[1])))

    expected = upper_count + lower_count
    if len(values) != expected:
        raise ValueError(
            f"{filename} declares {expected} points but contains {len(values)}"
        )
    return Airfoil(
        title=title,
        upper=values[:upper_count],
        lower=values[upper_count:],
    )


def write_airfoil(filename: Path, airfoil: Airfoil) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w") as stream:
        stream.write(f" {airfoil.title}\n")
        stream.write(f"      {len(airfoil.upper)}.       {len(airfoil.lower)}.\n\n")
        for point in airfoil.upper:
            stream.write(f" {point.x:.7f} {point.y:.7f}\n")
        stream.write("\n")
        for point in airfoil.lower:
            stream.write(f" {point.x:.7f} {point.y:.7f}\n")


def copy_airfoil(airfoil: Airfoil) -> Airfoil:
    return Airfoil(
        title=airfoil.title,
        upper=[Point(point.x, point.y) for point in airfoil.upper],
        lower=[Point(point.x, point.y) for point in airfoil.lower],
    )


def clip(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def ensure_matching_x(airfoil: Airfoil) -> None:
    if len(airfoil.upper) != len(airfoil.lower):
        raise ValueError("upper and lower surfaces must have the same point count")
    for index, (upper, lower) in enumerate(zip(airfoil.upper, airfoil.lower)):
        if abs(upper.x - lower.x) > 1.0e-12:
            raise ValueError(
                f"upper/lower x coordinates differ at data index {index}"
            )


def safe_shift(
    airfoil: Airfoil,
    original: Airfoil,
    surface: str,
    requested_shift: float,
    weights: list[float],
    hard_limit: float,
) -> tuple[float, str | None]:
    effective = clip(requested_shift, -hard_limit, hard_limit)
    reason = None
    if effective != requested_shift:
        reason = f"hard limit ±{hard_limit:g}"

    # Moving U downward or L upward reduces local thickness. Retain at least
    # 30% of the thickness that existed before this action.
    reduces_thickness = (
        surface == "U" and effective < 0.0
    ) or (
        surface == "L" and effective > 0.0
    )
    if reduces_thickness:
        admissible = hard_limit
        for upper, lower, original_upper, original_lower, weight in zip(
            airfoil.upper,
            airfoil.lower,
            original.upper,
            original.lower,
            weights,
        ):
            if weight <= 1.0e-10:
                continue
            thickness = upper.y - lower.y
            if thickness <= 0.0:
                raise ValueError("working airfoil already has crossed surfaces")
            original_thickness = original_upper.y - original_lower.y
            minimum_thickness = 0.30 * original_thickness
            available_thickness = thickness - minimum_thickness
            admissible = min(
                admissible,
                max(0.0, available_thickness) / weight,
            )
        clipped_magnitude = min(abs(effective), admissible)
        thickness_limited = math.copysign(clipped_magnitude, effective)
        if thickness_limited != effective:
            effective = thickness_limited
            reason = "local-thickness safety limit"
    return effective, reason


def sanitize_working_geometry(
    working: Airfoil,
    original: Airfoil,
    total_limit: float,
) -> tuple[Airfoil, int]:
    ensure_matching_x(working)
    ensure_matching_x(original)
    if len(working.upper) != len(original.upper):
        raise ValueError("working and original airfoil sizes differ")

    sanitized = copy_airfoil(working)
    clipped_points = 0
    for current_surface, original_surface in (
        (sanitized.upper, original.upper),
        (sanitized.lower, original.lower),
    ):
        for current, reference in zip(current_surface, original_surface):
            if abs(current.x - reference.x) > 1.0e-12:
                raise ValueError("working data changed an original x coordinate")
            if not math.isfinite(current.y):
                raise ValueError("working data contains a non-finite y coordinate")
            bounded = clip(
                current.y,
                reference.y - total_limit,
                reference.y + total_limit,
            )
            if bounded != current.y:
                current.y = bounded
                clipped_points += 1

    # Endpoints stay exactly at the original data values.
    for current_surface, original_surface in (
        (sanitized.upper, original.upper),
        (sanitized.lower, original.lower),
    ):
        current_surface[0].y = original_surface[0].y
        current_surface[-1].y = original_surface[-1].y

    for index, (upper, lower, original_upper, original_lower) in enumerate(
        zip(
            sanitized.upper,
            sanitized.lower,
            original.upper,
            original.lower,
        )
    ):
        original_thickness = original_upper.y - original_lower.y
        minimum_thickness = 0.30 * original_thickness
        if upper.y - lower.y + 1.0e-14 < minimum_thickness:
            raise ValueError(
                "persistent geometry is too thin at data index "
                f"{index}; run with --reset"
            )
    return sanitized, clipped_points


def endpoint_taper(x: float, taper_length: float = 0.05) -> float:
    """C1 taper that is zero at both endpoints and one in the interior."""
    if x <= 0.0 or x >= 1.0:
        return 0.0
    if x < taper_length:
        return 0.5 * (1.0 - math.cos(math.pi * x / taper_length))
    if x > 1.0 - taper_length:
        return 0.5 * (
            1.0
            - math.cos(math.pi * (1.0 - x) / taper_length)
        )
    return 1.0


def write_update_csv(
    filename: Path,
    before: Airfoil,
    after: Airfoil,
    selected_surface: str,
    center: float,
    width: float,
    requested_shift: float,
    effective_shift: float,
) -> None:
    with filename.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "surface",
                "index",
                "x",
                "y_before",
                "y_after",
                "delta_y",
                "selected_surface",
                "center",
                "width",
                "requested_shift",
                "effective_shift",
            ]
        )
        for surface, old_points, new_points in (
            ("U", before.upper, after.upper),
            ("L", before.lower, after.lower),
        ):
            for index, (old, new) in enumerate(zip(old_points, new_points)):
                writer.writerow(
                    [
                        surface,
                        index,
                        f"{old.x:.16g}",
                        f"{old.y:.16g}",
                        f"{new.y:.16g}",
                        f"{new.y - old.y:.16g}",
                        selected_surface,
                        f"{center:.16g}",
                        f"{width:.16g}",
                        f"{requested_shift:.16g}",
                        f"{effective_shift:.16g}",
                    ]
                )


def show_figure(filename: Path) -> None:
    system = platform.system()
    if system == "Darwin":
        subprocess.run(["open", str(filename)], check=False)
    elif system == "Linux" and shutil.which("xdg-open"):
        subprocess.run(["xdg-open", str(filename)], check=False)
    elif system == "Windows":
        os.startfile(filename)  # type: ignore[attr-defined]
    else:
        print(f"Open this figure manually: {filename}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Update one airfoil surface with a Gaussian y-displacement, "
            "then refit, move, smooth, and plot the AFEPack mesh."
        )
    )
    parser.add_argument(
        "center",
        type=float,
        nargs="?",
        help="Gaussian center x/c; clipped to [0.05, 0.95]",
    )
    parser.add_argument(
        "surface",
        nargs="?",
        type=str.upper,
        choices=("U", "L"),
        help="surface to update: U (upper) or L (lower)",
    )
    parser.add_argument(
        "shift",
        type=float,
        nargs="?",
        help="requested peak y/c shift; clipped to a safe range",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=0.12,
        help="Gaussian standard deviation in x/c; clipped to [0.04, 0.30]",
    )
    parser.add_argument(
        "--max-shift",
        type=float,
        default=0.10,
        help="hard absolute shift limit; itself clipped to [0.002, 0.10]",
    )
    parser.add_argument(
        "--max-total-shift",
        type=float,
        default=0.10,
        help=(
            "maximum accumulated y displacement from the original data; "
            "clipped to [0.01, 0.10]"
        ),
    )
    parser.add_argument(
        "--smooth-iterations",
        type=int,
        default=None,
        help=(
            "number of Laplacian mesh-smoothing sweeps; clipped to [1, 1000]. "
            "When supplied without center/surface/shift, smooth the current "
            "mesh without changing the airfoil."
        ),
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help="restore the persistent working data from the original backup",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="open the generated overview figure after the run",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / "data"
    output_dir = script_dir / "output"
    original_file = data_dir / "naca0012_original.dat"
    compatibility_original = data_dir / "naca0012.dat"
    working_file = data_dir / "naca0012_working.dat"

    if not original_file.exists():
        shutil.copy2(compatibility_original, original_file)

    if arguments.reset:
        shutil.copy2(original_file, working_file)
        for persistent_name in (
            "mesh_current.mesh",
            "boundary_current.dat",
            "mesh_update_mode.txt",
        ):
            persistent_file = output_dir / persistent_name
            if persistent_file.exists():
                persistent_file.unlink()
        print(f"Working data restored from {original_file}")
        print("Persistent mesh cleared; the next action creates it once.")
        print(f"Persistent working file: {working_file}")
        return

    geometry_values = (
        arguments.center,
        arguments.surface,
        arguments.shift,
    )
    has_any_geometry_argument = any(
        value is not None for value in geometry_values
    )
    has_all_geometry_arguments = all(
        value is not None for value in geometry_values
    )
    if has_any_geometry_argument and not has_all_geometry_arguments:
        raise SystemExit(
            "center, surface, and shift must be supplied together"
        )

    if not has_any_geometry_argument and arguments.smooth_iterations is None:
        raise SystemExit(
            "supply center surface shift, --smooth-iterations, or --reset"
        )

    requested_smooth_iterations = (
        50
        if arguments.smooth_iterations is None
        else arguments.smooth_iterations
    )
    smooth_iterations = max(
        1,
        min(1000, requested_smooth_iterations),
    )

    if not has_any_geometry_argument:
        current_mesh = output_dir / "mesh_current.mesh"
        current_boundary = output_dir / "boundary_current.dat"
        if not current_mesh.exists() or not current_boundary.exists():
            raise SystemExit(
                "no persistent mesh is available; run one geometry action "
                "before using smoothing-only mode"
            )

        subprocess.run(
            [
                str(script_dir / "run.sh"),
                "",
                str(output_dir),
                "",
                str(smooth_iterations),
                "smooth-only",
            ],
            cwd=script_dir,
            env=os.environ.copy(),
            check=True,
        )

        print()
        print("Smoothing-only action completed")
        print("  airfoil geometry: unchanged")
        print("  EasyMesh:         not called")
        print(f"  smoothing sweeps: {smooth_iterations}")
        if smooth_iterations != requested_smooth_iterations:
            print("  smoothing iterations were clipped to [1, 1000]")
        print(f"  current mesh:     {current_mesh}")
        print(
            "  overview figure:  "
            f"{output_dir / 'figures/mesh_motion_overview.png'}"
        )
        if arguments.show:
            show_figure(
                output_dir / "figures" / "mesh_motion_overview.png"
            )
        return

    if not working_file.exists():
        shutil.copy2(original_file, working_file)

    original = read_airfoil(original_file)
    working = read_airfoil(working_file)
    total_limit = clip(abs(arguments.max_total_shift), 0.01, 0.10)
    before, repaired_points = sanitize_working_geometry(
        working,
        original,
        total_limit,
    )

    center = clip(arguments.center, 0.05, 0.95)
    width = clip(arguments.width, 0.04, 0.30)
    hard_limit = clip(abs(arguments.max_shift), 0.002, 0.10)
    weights = [
        math.exp(-0.5 * ((point.x - center) / width) ** 2)
        * endpoint_taper(point.x)
        for point in before.upper
    ]

    effective_shift, clip_reason = safe_shift(
        before,
        original,
        arguments.surface,
        arguments.shift,
        weights,
        hard_limit,
    )
    after = copy_airfoil(before)
    selected = after.upper if arguments.surface == "U" else after.lower
    original_selected = (
        original.upper if arguments.surface == "U" else original.lower
    )
    total_clipped_points = 0
    for point, reference, weight in zip(
        selected,
        original_selected,
        weights,
    ):
        requested_y = point.y + effective_shift * weight
        bounded_y = clip(
            requested_y,
            reference.y - total_limit,
            reference.y + total_limit,
        )
        if bounded_y != requested_y:
            total_clipped_points += 1
        point.y = bounded_y
    selected[0].y = original_selected[0].y
    selected[-1].y = original_selected[-1].y

    output_dir.mkdir(parents=True, exist_ok=True)
    previous_file = output_dir / "airfoil_previous.dat"
    candidate_file = output_dir / "airfoil_candidate.dat"
    update_file = output_dir / "last_update.csv"
    write_airfoil(previous_file, before)
    write_airfoil(candidate_file, after)
    write_update_csv(
        update_file,
        before,
        after,
        arguments.surface,
        center,
        width,
        arguments.shift,
        effective_shift,
    )

    environment = os.environ.copy()
    subprocess.run(
        [
            str(script_dir / "run.sh"),
            str(previous_file),
            str(output_dir),
            str(candidate_file),
            str(smooth_iterations),
        ],
        cwd=script_dir,
        env=environment,
        check=True,
    )

    # Commit the new geometry only after the full mesh workflow succeeds.
    shutil.copy2(candidate_file, working_file)

    closest_index = min(
        range(len(before.upper)),
        key=lambda index: abs(before.upper[index].x - center),
    )
    old_selected = before.upper if arguments.surface == "U" else before.lower
    new_selected = after.upper if arguments.surface == "U" else after.lower
    actual_peak = (
        new_selected[closest_index].y - old_selected[closest_index].y
    )
    print()
    print("Persistent Gaussian airfoil action completed")
    print(f"  surface:          {arguments.surface}")
    print(f"  requested center: {arguments.center:.6g}")
    print(f"  effective center: {center:.6g}")
    print(f"  requested shift:  {arguments.shift:.6g}")
    print(f"  effective shift:  {effective_shift:.6g}")
    print(
        "  closest data x:   "
        f"{before.upper[closest_index].x:.7f}, "
        f"actual delta y = {actual_peak:.7g}"
    )
    if center != arguments.center:
        print("  center was clipped to [0.05, 0.95]")
    if width != arguments.width:
        print("  width was clipped to [0.04, 0.30]")
    if smooth_iterations != requested_smooth_iterations:
        print("  smoothing iterations were clipped to [1, 1000]")
    if clip_reason:
        print(f"  shift clipped by: {clip_reason}")
    print(f"  smoothing sweeps: {smooth_iterations}")
    if repaired_points:
        print(
            "  repaired legacy points outside accumulated limit: "
            f"{repaired_points}"
        )
    if total_clipped_points:
        print(
            "  points clipped by accumulated deformation limit: "
            f"{total_clipped_points}"
        )
    print(f"  working data:     {working_file}")
    print(f"  overview figure:  {output_dir / 'figures/mesh_motion_overview.png'}")

    if arguments.show:
        show_figure(output_dir / "figures" / "mesh_motion_overview.png")


if __name__ == "__main__":
    main()
