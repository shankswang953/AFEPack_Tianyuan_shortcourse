"""Test ML-iteration boundary-curve point-cloud behavior.

Run from `08_airfoil_ml_iteration/` with:
    python3 -m unittest tests.test_pointcloud -v

The test uses in-memory NumPy arrays and creates no output files.
"""

from __future__ import annotations

import unittest

import numpy as np

from ml_airfoil.pointcloud import BoundaryCurves


class PointCloudTest(unittest.TestCase):
    def test_closed_boundary_omits_duplicate_lower_endpoints(self) -> None:
        upper = np.asarray([[0.0, 0.0], [0.5, 0.1], [1.0, 0.0]])
        lower = np.asarray([[0.0, 0.0], [0.5, -0.1], [1.0, 0.0]])
        curves = BoundaryCurves(upper=upper, lower=lower)
        self.assertEqual(curves.all_points.shape, (4, 2))


if __name__ == "__main__":
    unittest.main()
