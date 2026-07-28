# Airfoil input and working data

- `naca0012_original.dat`: immutable backup used by `airfoil_step.py --reset`.
- `naca0012.dat`: default input used by the lower-level launcher.
- `naca0012_working.dat`: mutable persistent shape used across Python actions.

The files use UIUC-style upper/lower airfoil coordinates. Generated meshes,
diagnostics, and figures are written to `../output/`, not this directory.
