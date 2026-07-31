"""Test ML-iteration action banks and policy checkpoint round trips.

Run from `08_airfoil_ml_iteration/` with:
    python3 -m unittest tests.test_models -v

Temporary checkpoints are created in an automatically removed temporary
directory; no project output is retained.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from ml_airfoil.geometry import Action
from ml_airfoil.learning import (
    PolicyPredictor,
    discrete_actions,
    discrete_thickness_reducing_actions,
    train_policy,
)
from ml_airfoil.replay import PolicyExample


class ModelTest(unittest.TestCase):
    def test_discrete_action_bank_covers_full_chord(self) -> None:
        actions = discrete_thickness_reducing_actions(center_count=3)
        self.assertEqual(len(actions), 30)
        self.assertEqual(
            {action.center for action in actions},
            {0.0, 0.5, 1.0},
        )
        self.assertTrue(
            all(
                action.shift < 0.0
                if action.surface == "U"
                else action.shift > 0.0
                for action in actions
            )
        )
        bidirectional = discrete_actions(center_count=3)
        self.assertEqual(len(bidirectional), 60)
        self.assertEqual(
            {abs(action.shift) for action in bidirectional},
            {0.005, 0.01, 0.02, 0.03, 0.04},
        )
        for surface in ("U", "L"):
            self.assertEqual(
                {
                    int(action.shift > 0.0)
                    for action in bidirectional
                    if action.surface == surface
                },
                {0, 1},
            )

    def test_policy_checkpoint_round_trip(self) -> None:
        examples = [
            PolicyExample.create(
                np.full(128, value, dtype=np.float32),
                Action(
                    "U" if index % 2 == 0 else "L",
                    0.2 + 0.15 * index,
                    -0.02 if index % 2 == 0 else 0.02,
                ),
            )
            for index, value in enumerate((0.0, 0.1, 0.2, 0.3))
        ]
        with tempfile.TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "policy.pt"
            train_policy(examples, checkpoint, epochs=10)
            action = PolicyPredictor(checkpoint).predict(
                np.zeros(128, dtype=np.float32)
            )
            self.assertIn(action.surface, {"U", "L"})
            self.assertGreaterEqual(action.center, 0.0)
            self.assertLessEqual(action.center, 1.0)
            self.assertGreaterEqual(action.shift, -0.04)
            self.assertLessEqual(action.shift, 0.04)


if __name__ == "__main__":
    unittest.main()
