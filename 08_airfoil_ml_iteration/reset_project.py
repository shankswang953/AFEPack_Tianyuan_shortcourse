#!/usr/bin/env python3

"""Remove generated ML-iteration experiment state while preserving source data.

Command-line usage:
    python3 reset_project.py [--clean-build]

The default reset deletes `output/` and `data/working.dat`.
`--clean-build` also deletes the compiled backend objects and executables.
The initial/target `dat` files and all source files are preserved.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def reset_generated_state(root: Path, *, clean_build: bool = False) -> None:
    generated_paths = [
        root / "output",
        root / "data" / "working.dat",
    ]
    if clean_build:
        generated_paths.extend(
            [
                root / "backend" / "generate_airfoil_geometry",
                root / "backend" / "generate_airfoil_geometry.o",
                root / "backend" / "move_and_smooth",
                root / "backend" / "move_and_smooth.o",
            ]
        )

    for path in generated_paths:
        if path.is_dir():
            shutil.rmtree(path)
            print(f"removed directory: {path}")
        elif path.exists():
            path.unlink()
            print(f"removed file:      {path}")

    print("reset complete")
    print("preserved:")
    print(f"  {root / 'data' / 'initial_circle.dat'}")
    print(f"  {root / 'data' / 'target_naca0012.dat'}")
    print("  all source files")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--clean-build",
        action="store_true",
        help="also remove compiled C++ objects and executables",
    )
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parent
    reset_generated_state(root, clean_build=arguments.clean_build)


if __name__ == "__main__":
    main()
