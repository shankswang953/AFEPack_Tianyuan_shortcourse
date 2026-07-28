"""Test digital-twin airfoil geometry, actions, safety, and state construction.

Run from `08_airfoil_digital_twin/` with:
    python3 -m unittest tests.test_geometry -v

The tests read checked-in files below `data/`, print results to the terminal,
and create no persistent output files.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from dt_airfoil.geometry import (
    Action,
    action_points_toward_target,
    apply_action,
    copy_airfoil,
    data_shape_mse,
    make_circle_from_x_grid,
    read_airfoil,
    state_vector,
)


ROOT = Path(__file__).resolve().parents[1]


class GeometryTest(unittest.TestCase):
    def test_data_shape_mse_uses_matching_interior_points(self) -> None:
        target = read_airfoil(ROOT / "data" / "target_naca0012.dat")
        current = copy_airfoil(target)
        self.assertEqual(data_shape_mse(current, target), 0.0)
        current.upper[10].y += 0.02
        movable_values = 2 * (len(target.upper) - 2)
        self.assertAlmostEqual(
            data_shape_mse(current, target),
            0.02**2 / movable_values,
            places=16,
        )

    def test_circle_and_action_keep_endpoints(self) -> None:
        target = read_airfoil(ROOT / "data" / "target_naca0012.dat")
        circle = make_circle_from_x_grid(target)
        moved, effective = apply_action(
            circle,
            Action("U", 0.5, -0.03),
        )
        self.assertEqual(state_vector(circle).size, 128)
        self.assertAlmostEqual(moved.upper[0].y, circle.upper[0].y)
        self.assertAlmostEqual(moved.upper[-1].y, circle.upper[-1].y)
        self.assertLess(effective.shift, 0.0)
        self.assertLess(moved.upper[32].y, circle.upper[32].y)
        far_index = min(
            range(len(circle.upper)),
            key=lambda index: abs(circle.upper[index].x - 0.95),
        )
        local_moved, _ = apply_action(
            circle,
            Action("U", 0.0, -0.03),
        )
        self.assertLess(
            local_moved.upper[1].y,
            circle.upper[1].y,
        )
        self.assertEqual(
            local_moved.upper[far_index].y,
            circle.upper[far_index].y,
        )

        for center in (0.0, 1.0):
            edge_moved, edge_action = apply_action(
                circle,
                Action("L", center, 0.03),
            )
            self.assertEqual(edge_action.center, center)
            self.assertAlmostEqual(
                edge_moved.lower[0].y,
                circle.lower[0].y,
            )
            self.assertAlmostEqual(
                edge_moved.lower[-1].y,
                circle.lower[-1].y,
            )

    def test_direction_mask_reverses_after_local_overshoot(self) -> None:
        target = read_airfoil(ROOT / "data" / "target_naca0012.dat")
        circle = make_circle_from_x_grid(target)
        self.assertTrue(
            action_points_toward_target(
                circle,
                target,
                Action("U", 0.5, -0.01),
            )
        )
        self.assertFalse(
            action_points_toward_target(
                circle,
                target,
                Action("U", 0.5, 0.01),
            )
        )
        self.assertTrue(
            action_points_toward_target(
                circle,
                target,
                Action("U", 0.0, -0.01),
            )
        )
        self.assertFalse(
            action_points_toward_target(
                circle,
                target,
                Action("U", 1.0, 0.01),
            )
        )
        overshot, _ = apply_action(
            target,
            Action("U", 0.5, -0.01),
        )
        self.assertTrue(
            action_points_toward_target(
                overshot,
                target,
                Action("U", 0.5, 0.01),
            )
        )
        self.assertFalse(
            action_points_toward_target(
                overshot,
                target,
                Action("U", 0.5, -0.01),
            )
        )

    def test_target_envelope_prevents_surface_overshoot(self) -> None:
        target = read_airfoil(ROOT / "data" / "target_naca0012.dat")
        shape = make_circle_from_x_grid(target)
        for _ in range(30):
            shape, _ = apply_action(
                shape,
                Action("U", 0.05, -0.04),
                thickness_reference=target,
            )
        for upper, lower, target_upper, target_lower in zip(
            shape.upper[1:-1],
            shape.lower[1:-1],
            target.upper[1:-1],
            target.lower[1:-1],
        ):
            self.assertGreaterEqual(
                upper.y + 1.0e-12,
                target_upper.y,
            )
            self.assertLessEqual(
                lower.y - 1.0e-12,
                target_lower.y,
            )


if __name__ == "__main__":
    unittest.main()
