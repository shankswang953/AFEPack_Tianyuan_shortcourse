#!/usr/bin/env python3

"""Build the backend and prepare independent initial/target reference meshes.

Command-line usage:
    python3 prepare_case.py [--force]

`--force` rebuilds both EasyMesh reference cases. Generated reference meshes
and the persistent working case are written below `output/reference/` and
`output/current/`; `data/working.dat` is the mutable shape. Configure
EASYMESH_BIN and EASYMESH2MESH_BIN on systems with non-default tool paths.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ml_airfoil.environment import AirfoilMeshEnvironment
from ml_airfoil.pointcloud import extract_boundary_curves


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force",
        action="store_true",
        help="rebuild both EasyMesh reference meshes",
    )
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parent
    environment = AirfoilMeshEnvironment(root)
    environment.prepare(force=arguments.force)
    state = environment.reset()
    initial = extract_boundary_curves(
        environment.initial_reference / "mesh_smoothed_nodes.csv",
        environment.initial_reference / "mesh_smoothed_elements.csv",
    )
    target = extract_boundary_curves(
        environment.target_reference / "mesh_smoothed_nodes.csv",
        environment.target_reference / "mesh_smoothed_elements.csv",
    )
    print("AFEPack ML-iteration case is ready")
    print(f"  state dimension:       {state.size}")
    print(f"  initial boundary:      {initial.all_points.shape[0]} mesh points")
    print(f"  target boundary:       {target.all_points.shape[0]} mesh points")
    print(f"  initial shape loss:    {environment.current_loss:.8e}")
    print(f"  smoothing iterations: {environment.smoothing_iterations}")
    print(f"  Gaussian width:       {environment.gaussian_width}")
    print(f"  working data:         {environment.working_dat}")


if __name__ == "__main__":
    main()

