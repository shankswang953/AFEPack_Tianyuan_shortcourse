"""Public imports for the AFEPack airfoil ML-iteration package.

This is a library module and has no command-line interface or direct outputs.
Import shared geometry/action types from `ml_airfoil`; run the scripts in the
parent directory for complete experiments.
"""

from .geometry import (
    Action,
    Airfoil,
    Point,
    data_shape_mse,
    read_airfoil,
    write_airfoil,
)
from .pointcloud import BoundaryCurves

__all__ = [
    "Action",
    "Airfoil",
    "BoundaryCurves",
    "Point",
    "data_shape_mse",
    "read_airfoil",
    "write_airfoil",
]
