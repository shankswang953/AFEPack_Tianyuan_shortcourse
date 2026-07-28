"""Test the pure-`dat` airfoil environment, actions, and safety constraints.

Run from `10_airfoil_rl_dat/` with:
    python3 -m unittest tests.test_environment -v

The tests read `data/target_naca0012.dat`, print results to the terminal, and
create no persistent output files.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from airfoil_rl.environment import DatAirfoilEnvironment
from airfoil_rl.geometry import (
    Action,
    apply_action,
    data_mse,
    make_circle_from_x_grid,
    read_airfoil,
)


ROOT = Path(__file__).resolve().parents[1]


class EnvironmentTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.target = read_airfoil(
            ROOT / "data" / "target_naca0012.dat"
        )
        cls.circle = make_circle_from_x_grid(cls.target)

    def test_action_grid_and_state_dimension(self) -> None:
        environment = DatAirfoilEnvironment(
            self.target,
            self.circle,
        )
        self.assertEqual(environment.action_count, 126)
        self.assertEqual(environment.state_dimension, 128)
        self.assertTrue(environment.valid_action_mask().all())

    def test_target_clipping_prevents_overshoot(self) -> None:
        moved = self.circle
        for _ in range(30):
            moved, _ = apply_action(
                moved,
                Action("U", 0.5, -0.04),
                target_reference=self.target,
            )
        self.assertTrue(
            all(
                upper >= target
                for upper, target in zip(
                    moved.upper_y[1:-1],
                    self.target.upper_y[1:-1],
                )
            )
        )

    def test_correct_direction_reduces_data_mse(self) -> None:
        before = data_mse(self.circle, self.target)
        moved, realized = apply_action(
            self.circle,
            Action("U", 0.5, -0.04),
        )
        self.assertGreater(realized, 0.0)
        self.assertLess(data_mse(moved, self.target), before)

    def test_wrong_direction_increases_data_mse(self) -> None:
        before = data_mse(self.circle, self.target)
        moved, _ = apply_action(
            self.circle,
            Action("U", 0.5, 0.04),
        )
        self.assertGreater(data_mse(moved, self.target), before)

    def test_surface_order_is_preserved(self) -> None:
        moved = self.circle
        for _ in range(30):
            moved, _ = apply_action(
                moved,
                Action("U", 0.5, -0.04),
            )
        self.assertTrue(
            all(
                upper > lower
                for upper, lower in zip(
                    moved.upper_y[1:-1],
                    moved.lower_y[1:-1],
                )
            )
        )


if __name__ == "__main__":
    unittest.main()
