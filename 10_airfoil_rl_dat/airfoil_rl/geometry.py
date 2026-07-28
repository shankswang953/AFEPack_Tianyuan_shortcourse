"""Provide airfoil `dat` I/O and local Gaussian boundary actions.

This import-only module has no command-line interface or fixed output path.
Its read/write functions use only caller-supplied paths.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Action:
    surface: str
    center: float
    shift: float


@dataclass
class Airfoil:
    title: str
    x: np.ndarray
    upper_y: np.ndarray
    lower_y: np.ndarray

    def copy(self, *, title: str | None = None) -> "Airfoil":
        return Airfoil(
            title=self.title if title is None else title,
            x=self.x.copy(),
            upper_y=self.upper_y.copy(),
            lower_y=self.lower_y.copy(),
        )


def read_airfoil(filename: Path) -> Airfoil:
    with filename.open() as stream:
        title = stream.readline().strip()
        counts = stream.readline().split()
        if len(counts) < 2:
            raise ValueError(f"invalid point-count line in {filename}")
        upper_count = int(round(float(counts[0].rstrip("."))))
        lower_count = int(round(float(counts[1].rstrip("."))))
        values = [
            (float(fields[0]), float(fields[1]))
            for line in stream
            if len(fields := line.split()) >= 2
        ]
    if len(values) != upper_count + lower_count:
        raise ValueError(
            f"{filename} declares {upper_count + lower_count} points "
            f"but contains {len(values)}"
        )
    upper = np.asarray(values[:upper_count], dtype=np.float64)
    lower = np.asarray(values[upper_count:], dtype=np.float64)
    if upper_count != lower_count:
        raise ValueError("upper and lower point counts differ")
    if not np.allclose(upper[:, 0], lower[:, 0], atol=1.0e-10):
        raise ValueError("upper and lower x grids differ")
    return Airfoil(
        title=title,
        x=upper[:, 0],
        upper_y=upper[:, 1],
        lower_y=lower[:, 1],
    )


def write_airfoil(filename: Path, airfoil: Airfoil) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w") as stream:
        stream.write(f" {airfoil.title}\n")
        stream.write(
            f"      {airfoil.x.size}.       {airfoil.x.size}.\n\n"
        )
        for x_value, y_value in zip(airfoil.x, airfoil.upper_y):
            stream.write(f" {x_value:.10f} {y_value:.10f}\n")
        stream.write("\n")
        for x_value, y_value in zip(airfoil.x, airfoil.lower_y):
            stream.write(f" {x_value:.10f} {y_value:.10f}\n")


def make_circle_from_x_grid(target: Airfoil) -> Airfoil:
    radius_squared = np.maximum(
        0.0,
        0.25 - (target.x - 0.5) ** 2,
    )
    radius = np.sqrt(radius_squared)
    radius[0] = 0.0
    radius[-1] = 0.0
    return Airfoil(
        title="CIRCULAR INITIAL SHAPE",
        x=target.x.copy(),
        upper_y=radius.copy(),
        lower_y=-radius,
    )


def data_mse(current: Airfoil, target: Airfoil) -> float:
    errors = np.concatenate(
        (
            current.upper_y[1:-1] - target.upper_y[1:-1],
            current.lower_y[1:-1] - target.lower_y[1:-1],
        )
    )
    return float(np.mean(errors**2))


def error_state(
    current: Airfoil,
    target: Airfoil,
    *,
    scale: float = 0.5,
) -> np.ndarray:
    return (
        np.concatenate(
            (
                current.upper_y[1:-1] - target.upper_y[1:-1],
                current.lower_y[1:-1] - target.lower_y[1:-1],
            )
        )
        / scale
    ).astype(np.float32)


def gaussian_weights(
    x_grid: np.ndarray,
    center: float,
    *,
    width: float = 0.12,
) -> np.ndarray:
    weights = np.exp(-0.5 * ((x_grid - center) / width) ** 2)
    weights[weights < 1.0e-2] = 0.0
    weights[0] = 0.0
    weights[-1] = 0.0
    return weights


def minimum_thickness(x_grid: np.ndarray) -> np.ndarray:
    chord_factor = np.sqrt(
        np.maximum(0.0, 4.0 * x_grid * (1.0 - x_grid))
    )
    return 2.0e-3 * chord_factor


def apply_action(
    airfoil: Airfoil,
    action: Action,
    *,
    width: float = 0.12,
    ordinate_limit: float = 0.55,
    target_reference: Airfoil | None = None,
) -> tuple[Airfoil, float]:
    """Apply a clipped Gaussian update and return its realized maximum shift."""

    surface = action.surface.upper()
    if surface not in {"U", "L"}:
        raise ValueError("surface must be U or L")
    center = float(np.clip(action.center, 0.0, 1.0))
    shift = float(np.clip(action.shift, -0.08, 0.08))
    weights = gaussian_weights(airfoil.x, center, width=width)
    moved = airfoil.copy()
    selected = moved.upper_y if surface == "U" else moved.lower_y
    selected += shift * weights
    np.clip(selected, -ordinate_limit, ordinate_limit, out=selected)
    if target_reference is not None:
        target_values = (
            target_reference.upper_y
            if surface == "U"
            else target_reference.lower_y
        )
        original_values = (
            airfoil.upper_y
            if surface == "U"
            else airfoil.lower_y
        )
        error_before = original_values - target_values
        error_after = selected - target_values
        crossing = (np.abs(error_before) <= 1.0e-14) | (
            error_before * error_after < 0.0
        )
        selected[crossing] = target_values[crossing]

    floor = minimum_thickness(moved.x)
    if surface == "U":
        selected[:] = np.maximum(selected, moved.lower_y + floor)
    else:
        selected[:] = np.minimum(selected, moved.upper_y - floor)
    selected[0] = airfoil.upper_y[0] if surface == "U" else airfoil.lower_y[0]
    selected[-1] = (
        airfoil.upper_y[-1] if surface == "U" else airfoil.lower_y[-1]
    )

    original = airfoil.upper_y if surface == "U" else airfoil.lower_y
    realized = float(np.max(np.abs(selected - original)))
    return moved, realized


def interpolate_shapes(
    target: Airfoil,
    initial: Airfoil,
    alpha: float,
) -> Airfoil:
    clipped = float(np.clip(alpha, 0.0, 1.0))
    return Airfoil(
        title=f"TRAINING RESET alpha={clipped:.3f}",
        x=target.x.copy(),
        upper_y=target.upper_y
        + clipped * (initial.upper_y - target.upper_y),
        lower_y=target.lower_y
        + clipped * (initial.lower_y - target.lower_y),
    )


def smooth_random_perturbation(
    airfoil: Airfoil,
    rng: np.random.Generator,
    *,
    amplitude: float,
) -> Airfoil:
    """Add a few smooth Gaussian bumps while preserving valid surfaces."""

    if amplitude <= 0.0:
        return airfoil.copy()
    perturbed = airfoil.copy()
    for surface in ("U", "L"):
        values = (
            perturbed.upper_y
            if surface == "U"
            else perturbed.lower_y
        )
        for _ in range(2):
            center = float(rng.uniform(0.05, 0.95))
            width = float(rng.uniform(0.08, 0.20))
            outward_sign = 1.0 if surface == "U" else -1.0
            shift = outward_sign * abs(float(rng.normal(0.0, amplitude)))
            values += shift * gaussian_weights(
                perturbed.x,
                center,
                width=width,
            )
    floor = minimum_thickness(perturbed.x)
    midpoint = 0.5 * (perturbed.upper_y + perturbed.lower_y)
    half_thickness = np.maximum(
        0.5 * (perturbed.upper_y - perturbed.lower_y),
        0.5 * floor,
    )
    perturbed.upper_y = np.clip(
        midpoint + half_thickness,
        -0.55,
        0.55,
    )
    perturbed.lower_y = np.clip(
        midpoint - half_thickness,
        -0.55,
        0.55,
    )
    perturbed.upper_y[[0, -1]] = airfoil.upper_y[[0, -1]]
    perturbed.lower_y[[0, -1]] = airfoil.lower_y[[0, -1]]
    return perturbed


def action_grid() -> list[Action]:
    centers = np.linspace(0.0, 1.0, 21)
    return [
        Action(surface, float(center), float(sign * magnitude))
        for surface, sign in (("U", -1.0), ("L", 1.0))
        for center in centers
        for magnitude in (0.01, 0.02, 0.04)
    ]
