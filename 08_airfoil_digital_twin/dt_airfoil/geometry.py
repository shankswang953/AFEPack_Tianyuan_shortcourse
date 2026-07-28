"""Provide airfoil `dat` I/O, loss/state calculations, and safe actions.

This import-only library has no command-line interface and no fixed output
path. Functions read or write only the `Path` objects supplied by callers.
Use the parent-directory scripts to run an experiment.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Point:
    x: float
    y: float


@dataclass
class Airfoil:
    title: str
    upper: list[Point]
    lower: list[Point]


@dataclass(frozen=True)
class Action:
    """One boundary update.

    ``surface`` is ``"U"`` or ``"L"``.  The Gaussian width is intentionally
    fixed by the experiment rather than learned.
    """

    surface: str
    center: float
    shift: float

    def normalized(
        self,
        center_min: float = 0.0,
        center_max: float = 1.0,
        shift_max: float = 0.08,
    ) -> "Action":
        surface = self.surface.upper()
        if surface not in {"U", "L"}:
            raise ValueError("surface must be U or L")
        return Action(
            surface=surface,
            center=float(np.clip(self.center, center_min, center_max)),
            shift=float(np.clip(self.shift, -shift_max, shift_max)),
        )


def read_airfoil(filename: Path) -> Airfoil:
    with filename.open() as stream:
        title = stream.readline().strip()
        counts = stream.readline().split()
        if len(counts) < 2:
            raise ValueError(f"invalid point-count line in {filename}")
        upper_count = int(round(float(counts[0].rstrip("."))))
        lower_count = int(round(float(counts[1].rstrip("."))))
        values: list[Point] = []
        for line in stream:
            fields = line.split()
            if len(fields) >= 2:
                values.append(Point(float(fields[0]), float(fields[1])))

    if len(values) != upper_count + lower_count:
        raise ValueError(
            f"{filename} declares {upper_count + lower_count} points "
            f"but contains {len(values)}"
        )
    airfoil = Airfoil(
        title=title,
        upper=values[:upper_count],
        lower=values[upper_count:],
    )
    ensure_matching_x(airfoil)
    return airfoil


def write_airfoil(filename: Path, airfoil: Airfoil) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w") as stream:
        stream.write(f" {airfoil.title}\n")
        stream.write(
            f"      {len(airfoil.upper)}.       "
            f"{len(airfoil.lower)}.\n\n"
        )
        for point in airfoil.upper:
            stream.write(f" {point.x:.10f} {point.y:.10f}\n")
        stream.write("\n")
        for point in airfoil.lower:
            stream.write(f" {point.x:.10f} {point.y:.10f}\n")


def copy_airfoil(airfoil: Airfoil) -> Airfoil:
    return Airfoil(
        title=airfoil.title,
        upper=[Point(point.x, point.y) for point in airfoil.upper],
        lower=[Point(point.x, point.y) for point in airfoil.lower],
    )


def ensure_matching_x(airfoil: Airfoil) -> None:
    if len(airfoil.upper) != len(airfoil.lower):
        raise ValueError("upper and lower point counts differ")
    for index, (upper, lower) in enumerate(
        zip(airfoil.upper, airfoil.lower)
    ):
        if abs(upper.x - lower.x) > 1.0e-10:
            raise ValueError(
                f"upper/lower x coordinates differ at index {index}"
            )


def make_circle_from_x_grid(target: Airfoil) -> Airfoil:
    """Use the target x-grid to make a radius-0.5 circle.

    The shared endpoints remain (0, 0) and (1, 0), so every later action can
    keep the leading and trailing edges fixed.
    """

    ensure_matching_x(target)
    upper: list[Point] = []
    lower: list[Point] = []
    for target_upper, target_lower in zip(target.upper, target.lower):
        x = 0.5 * (target_upper.x + target_lower.x)
        radius_squared = max(0.0, 0.25 - (x - 0.5) ** 2)
        y = math.sqrt(radius_squared)
        upper.append(Point(x, y))
        lower.append(Point(x, -y))
    return Airfoil("CIRCULAR INITIAL SHAPE", upper, lower)


def _action_weights(
    airfoil: Airfoil,
    center: float,
    width: float,
) -> list[float]:
    """Gaussian weights on movable points; only the two endpoints are fixed."""

    weights = []
    for index, point in enumerate(airfoil.upper):
        if index == 0 or index == len(airfoil.upper) - 1:
            weights.append(0.0)
            continue
        weight = math.exp(
            -0.5 * ((point.x - center) / width) ** 2
        )
        weights.append(0.0 if weight < 1.0e-2 else weight)
    return weights


def minimum_thickness(x: float) -> float:
    """Small positivity guard; it does not prescribe the target thickness."""

    chord_factor = math.sqrt(max(0.0, 4.0 * x * (1.0 - x)))
    return 2.0e-3 * chord_factor


def surface_ordinate(
    airfoil: Airfoil,
    surface: str,
    x: float,
) -> float:
    """Interpolate one data-curve ordinate at a chordwise position."""

    normalized_surface = surface.upper()
    if normalized_surface not in {"U", "L"}:
        raise ValueError("surface must be U or L")
    points = (
        airfoil.upper
        if normalized_surface == "U"
        else airfoil.lower
    )
    coordinates = sorted((point.x, point.y) for point in points)
    return float(
        np.interp(
            float(np.clip(x, coordinates[0][0], coordinates[-1][0])),
            [coordinate[0] for coordinate in coordinates],
            [coordinate[1] for coordinate in coordinates],
        )
    )


def action_points_toward_target(
    current: Airfoil,
    target: Airfoil,
    action: Action,
    *,
    width: float = 0.12,
    tolerance: float = 1.0e-10,
) -> bool:
    """Return whether the action is a descent direction for data-file MSE.

    This is only a direction mask.  The full data-file MSE still decides
    whether the proposed action is accepted.
    """

    if not math.isclose(width, 0.12):
        raise ValueError("this example fixes the Gaussian width at 0.12")
    ensure_matching_x(current)
    ensure_matching_x(target)
    normalized = action.normalized()
    current_points = (
        current.upper
        if normalized.surface == "U"
        else current.lower
    )
    target_points = (
        target.upper
        if normalized.surface == "U"
        else target.lower
    )
    weights = _action_weights(current, normalized.center, width)
    weighted_error = sum(
        weight * (target_point.y - current_point.y)
        for current_point, target_point, weight in zip(
            current_points,
            target_points,
            weights,
        )
    )
    if abs(weighted_error) <= tolerance:
        return False
    return normalized.shift * weighted_error > 0.0


def data_shape_mse(current: Airfoil, target: Airfoil) -> float:
    """Mean squared ordinate difference on the shared interior x-grid."""

    ensure_matching_x(current)
    ensure_matching_x(target)
    if len(current.upper) != len(target.upper):
        raise ValueError("current and target point counts differ")
    for index, (current_point, target_point) in enumerate(
        zip(current.upper, target.upper)
    ):
        if abs(current_point.x - target_point.x) > 1.0e-10:
            raise ValueError(
                f"current and target x coordinates differ at index {index}"
            )

    errors = np.asarray(
        [
            point.y - target_point.y
            for point, target_point in zip(
                current.upper[1:-1],
                target.upper[1:-1],
            )
        ]
        + [
            point.y - target_point.y
            for point, target_point in zip(
                current.lower[1:-1],
                target.lower[1:-1],
            )
        ],
        dtype=np.float64,
    )
    if errors.size == 0:
        raise ValueError("airfoil has no movable interior points")
    return float(np.mean(errors**2))


def apply_action(
    airfoil: Airfoil,
    action: Action,
    *,
    width: float = 0.12,
    shift_max: float = 0.08,
    thickness_reference: Airfoil | None = None,
) -> tuple[Airfoil, Action]:
    """Apply one safe Gaussian update and return its effective action."""

    if not math.isclose(width, 0.12):
        raise ValueError("this example fixes the Gaussian width at 0.12")
    ensure_matching_x(airfoil)
    if thickness_reference is not None:
        ensure_matching_x(thickness_reference)
        if len(thickness_reference.upper) != len(airfoil.upper):
            raise ValueError("airfoil and thickness reference sizes differ")
        for index, (point, reference_point) in enumerate(
            zip(airfoil.upper, thickness_reference.upper)
        ):
            if abs(point.x - reference_point.x) > 1.0e-10:
                raise ValueError(
                    "airfoil and thickness reference x coordinates "
                    f"differ at index {index}"
                )
    requested = action.normalized(shift_max=shift_max)
    weights = _action_weights(airfoil, requested.center, width)

    effective_shift = requested.shift
    reduces_thickness = (
        requested.surface == "U" and effective_shift < 0.0
    ) or (
        requested.surface == "L" and effective_shift > 0.0
    )
    if reduces_thickness:
        admissible = shift_max
        for upper, lower, weight in zip(
            airfoil.upper,
            airfoil.lower,
            weights,
        ):
            if weight <= 1.0e-12:
                continue
            thickness_floor = minimum_thickness(upper.x)
            remaining = (
                upper.y - lower.y - thickness_floor
            )
            admissible = min(admissible, max(0.0, remaining) / weight)
        effective_shift = math.copysign(
            min(abs(effective_shift), admissible),
            effective_shift,
        )

    moved = copy_airfoil(airfoil)
    selected = moved.upper if requested.surface == "U" else moved.lower
    original_selected = (
        airfoil.upper if requested.surface == "U" else airfoil.lower
    )
    reference_selected = (
        None
        if thickness_reference is None
        else (
            thickness_reference.upper
            if requested.surface == "U"
            else thickness_reference.lower
        )
    )
    for index, (point, original_point, weight) in enumerate(
        zip(selected, original_selected, weights)
    ):
        proposed_y = original_point.y + effective_shift * weight
        if reference_selected is not None:
            target_y = reference_selected[index].y
            error_before = original_point.y - target_y
            error_after = proposed_y - target_y
            if (
                abs(error_before) <= 1.0e-14
                or error_before * error_after < 0.0
            ):
                proposed_y = target_y
        point.y = proposed_y

    realized_shifts = [
        (point.y - original_point.y) / weight
        for point, original_point, weight in zip(
            selected,
            original_selected,
            weights,
        )
        if weight > 1.0e-12
    ]
    effective_shift = (
        max(realized_shifts, key=abs)
        if realized_shifts
        else 0.0
    )

    # The action window is already zero at both ends.  Assign the values
    # explicitly so accumulated floating-point noise can never move them.
    selected[0].y = (
        airfoil.upper[0].y
        if requested.surface == "U"
        else airfoil.lower[0].y
    )
    selected[-1].y = (
        airfoil.upper[-1].y
        if requested.surface == "U"
        else airfoil.lower[-1].y
    )

    for index, (upper, lower) in enumerate(
        zip(moved.upper[1:-1], moved.lower[1:-1]),
        start=1,
    ):
        if upper.y <= lower.y:
            raise ValueError(f"action crosses the surfaces at index {index}")

    effective = Action(
        surface=requested.surface,
        center=requested.center,
        shift=effective_shift,
    )
    return moved, effective


def state_vector(airfoil: Airfoil) -> np.ndarray:
    """The 128 trainable y values; endpoint coordinates are omitted."""

    return np.asarray(
        [point.y for point in airfoil.upper[1:-1]]
        + [point.y for point in airfoil.lower[1:-1]],
        dtype=np.float32,
    )


def action_vector(action: Action, shift_max: float = 0.08) -> np.ndarray:
    """Continuous encoding used by the digital-twin network."""

    surface = 1.0 if action.surface == "U" else -1.0
    center = 2.0 * action.center - 1.0
    shift = action.shift / shift_max
    return np.asarray([surface, center, shift], dtype=np.float32)
