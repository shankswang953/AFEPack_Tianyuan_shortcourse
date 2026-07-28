# Digital-twin airfoil data

- `initial_circle.dat`: checked-in initial geometry;
- `target_naca0012.dat`: checked-in target geometry;
- `working.dat`: generated mutable state for the current run.

`reset_project.py` deletes only `working.dat` and preserves both source
shapes. Meshes, replay data, models, and figures are written to `../output/`.
