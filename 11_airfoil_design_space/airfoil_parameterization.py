#!/usr/bin/env python3

"""Five-parameter, NACA-inspired airfoil geometry and validity checks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


PARAMETER_NAMES = ("m", "x_c", "t", "x_t", "t_te")

# The proposal box is deliberately wider than the engineering envelope.  This
# makes rejection part of the example rather than assuming every proposal is
# already suitable for CFD.
PROPOSAL_BOUNDS = np.array(
    [
        [-0.060, 0.140],  # signed maximum camber / chord
        [0.050, 0.900],   # location of maximum camber / chord
        [0.030, 0.320],   # maximum thickness / chord
        [0.120, 0.580],   # location of maximum thickness / chord
        [-0.020, 0.050],  # trailing-edge thickness / chord
    ],
    dtype=float,
)

ENGINEERING_BOUNDS = np.array(
    [
        [-0.030, 0.100],
        [0.120, 0.750],
        [0.060, 0.240],
        [0.200, 0.480],
        [0.000, 0.025],
    ],
    dtype=float,
)


@dataclass(frozen=True)
class AirfoilParameters:
    m: float
    x_c: float
    t: float
    x_t: float
    t_te: float

    @classmethod
    def from_array(cls, values: np.ndarray) -> "AirfoilParameters":
        return cls(*(float(value) for value in values))

    def as_array(self) -> np.ndarray:
        return np.array([self.m, self.x_c, self.t, self.x_t, self.t_te])


@dataclass
class AirfoilGeometry:
    upper: np.ndarray
    lower: np.ndarray

    def contour(self) -> np.ndarray:
        # Conventional orientation: upper LE -> TE, then lower TE -> LE.
        return np.vstack((self.upper, self.lower[-2:0:-1]))


@dataclass(frozen=True)
class ValidationResult:
    valid: bool
    reasons: tuple[str, ...]
    signed_area: float
    minimum_gap: float
    maximum_curvature: float


def cosine_grid(points_per_surface: int = 66) -> np.ndarray:
    theta = np.linspace(0.0, np.pi, points_per_surface)
    return 0.5 * (1.0 - np.cos(theta))


def _camber_line(x: np.ndarray, maximum: float, location: float) -> np.ndarray:
    before = maximum / location**2 * (2.0 * location * x - x**2)
    after = maximum / (1.0 - location) ** 2 * (
        (1.0 - 2.0 * location) + 2.0 * location * x - x**2
    )
    return np.where(x <= location, before, after)


def _reference_thickness(x: np.ndarray, thickness: float) -> np.ndarray:
    # Closed-trailing-edge NACA four-digit thickness polynomial.  It returns
    # half thickness; the coefficient -0.1036 closes the trailing edge.
    return 5.0 * thickness * (
        0.2969 * np.sqrt(np.maximum(x, 0.0))
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1036 * x**4
    )


def _warp_for_thickness_location(x: np.ndarray, target: float) -> np.ndarray:
    # The reference polynomial reaches maximum thickness near 0.30c.  A
    # monotone piecewise-linear coordinate warp moves that station to x_t.
    reference = 0.30
    left = reference * x / target
    right = reference + (1.0 - reference) * (x - target) / (1.0 - target)
    return np.where(x <= target, left, right)


def generate_airfoil(
    parameters: AirfoilParameters,
    points_per_surface: int = 66,
) -> AirfoilGeometry:
    x = cosine_grid(points_per_surface)
    camber = _camber_line(x, parameters.m, parameters.x_c)
    warped_x = _warp_for_thickness_location(x, parameters.x_t)
    half_thickness = _reference_thickness(warped_x, parameters.t)
    # Add a localized trailing-edge gap without perturbing the forward part.
    half_thickness += 0.5 * parameters.t_te * x**8
    upper = np.column_stack((x, camber + half_thickness))
    lower = np.column_stack((x, camber - half_thickness))
    return AirfoilGeometry(upper=upper, lower=lower)


def signed_area(polygon: np.ndarray) -> float:
    following = np.roll(polygon, -1, axis=0)
    return 0.5 * float(
        np.sum(polygon[:, 0] * following[:, 1]
               - following[:, 0] * polygon[:, 1])
    )


def _orientation(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return float(
        (b[0] - a[0]) * (c[1] - a[1])
        - (b[1] - a[1]) * (c[0] - a[0])
    )


def _segments_intersect(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    d: np.ndarray,
    tolerance: float = 1.0e-12,
) -> bool:
    o1 = _orientation(a, b, c)
    o2 = _orientation(a, b, d)
    o3 = _orientation(c, d, a)
    o4 = _orientation(c, d, b)
    return o1 * o2 < -tolerance and o3 * o4 < -tolerance


def has_self_intersection(polygon: np.ndarray) -> bool:
    count = len(polygon)
    for first in range(count):
        a = polygon[first]
        b = polygon[(first + 1) % count]
        for second in range(first + 1, count):
            if second in (first, first + 1):
                continue
            if first == 0 and second == count - 1:
                continue
            c = polygon[second]
            d = polygon[(second + 1) % count]
            if _segments_intersect(a, b, c, d):
                return True
    return False


def maximum_curvature(geometry: AirfoilGeometry) -> float:
    values: list[float] = []
    for surface in (geometry.upper, geometry.lower):
        x = surface[:, 0]
        y = surface[:, 1]
        first = np.gradient(y, x, edge_order=2)
        second = np.gradient(first, x, edge_order=2)
        curvature = np.abs(second) / np.power(1.0 + first**2, 1.5)
        interior = (x >= 0.02) & (x <= 0.98)
        values.append(float(np.max(curvature[interior])))
    return max(values)


def validate_airfoil(
    geometry: AirfoilGeometry,
    parameters: AirfoilParameters | None = None,
) -> ValidationResult:
    reasons: list[str] = []

    if parameters is not None:
        values = parameters.as_array()
        for index, name in enumerate(PARAMETER_NAMES):
            low, high = ENGINEERING_BOUNDS[index]
            if not low <= values[index] <= high:
                reasons.append(f"{name}_outside_engineering_range")

    for name, surface in (("upper", geometry.upper), ("lower", geometry.lower)):
        if not np.all(np.isfinite(surface)):
            reasons.append(f"{name}_nonfinite")
        if np.any(np.diff(surface[:, 0]) <= 0.0):
            reasons.append(f"{name}_x_not_increasing")

    if len(geometry.upper) != len(geometry.lower):
        reasons.append("surface_sizes_differ")
        gap = np.array([-np.inf])
    elif not np.allclose(geometry.upper[:, 0], geometry.lower[:, 0]):
        reasons.append("surface_x_grids_differ")
        gap = np.array([-np.inf])
    else:
        gap = geometry.upper[:, 1] - geometry.lower[:, 1]
        if np.any(gap[1:-1] <= 1.0e-5):
            reasons.append("upper_lower_crossing")

    contour = geometry.contour()
    area = signed_area(contour)
    if area >= 0.0:
        reasons.append("reversed_boundary_orientation")
    if not 0.035 <= abs(area) <= 0.260:
        reasons.append("area_outside_range")
    # With two strictly x-monotone surfaces on one x grid, a positive gap is
    # already a complete no-crossing certificate.  Fall back to the general
    # segment test for externally supplied or malformed data.
    monotone_common_grid = (
        len(geometry.upper) == len(geometry.lower)
        and np.all(np.diff(geometry.upper[:, 0]) > 0.0)
        and np.all(np.diff(geometry.lower[:, 0]) > 0.0)
        and np.allclose(geometry.upper[:, 0], geometry.lower[:, 0])
        and len(gap) > 2
        and np.all(gap[1:-1] > 1.0e-5)
    )
    if not monotone_common_grid and has_self_intersection(contour):
        reasons.append("self_intersection")

    if (
        np.all(np.diff(geometry.upper[:, 0]) != 0.0)
        and np.all(np.diff(geometry.lower[:, 0]) != 0.0)
    ):
        curvature = maximum_curvature(geometry)
    else:
        curvature = float("inf")
    if curvature > 85.0:
        reasons.append("excessive_curvature")
    if np.max(np.abs(contour[:, 1])) > 0.30:
        reasons.append("vertical_extent_outside_range")

    return ValidationResult(
        valid=not reasons,
        reasons=tuple(reasons),
        signed_area=area,
        minimum_gap=(
            float(np.min(gap[1:-1])) if len(gap) > 2 else float("-inf")
        ),
        maximum_curvature=curvature,
    )


def write_uiuc(filename: Path, geometry: AirfoilGeometry, title: str) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", encoding="utf-8") as stream:
        stream.write(f"{title}\n")
        stream.write(f"{len(geometry.upper)}. {len(geometry.lower)}.\n\n")
        for x, y in geometry.upper:
            stream.write(f" {x:.10f} {y:.10f}\n")
        stream.write("\n")
        for x, y in geometry.lower:
            stream.write(f" {x:.10f} {y:.10f}\n")


def read_uiuc(filename: Path) -> AirfoilGeometry:
    lines = [line.strip() for line in filename.read_text().splitlines()]
    nonempty = [line for line in lines if line]
    if len(nonempty) < 3:
        raise ValueError(f"{filename}: incomplete UIUC data")
    counts = nonempty[1].replace(".", "").split()
    upper_count, lower_count = int(counts[0]), int(counts[1])
    values = np.array(
        [[float(value) for value in line.split()[:2]] for line in nonempty[2:]],
        dtype=float,
    )
    if len(values) != upper_count + lower_count:
        raise ValueError(f"{filename}: point count does not match header")
    return AirfoilGeometry(
        upper=values[:upper_count],
        lower=values[upper_count:],
    )
