# Barycentric path endpoints

`initial_circle.dat` and `target_naca0012.dat` are the checked-in endpoints
of the pointwise barycentric path. They must have identical upper/lower
chordwise x-grids.

The scripts treat these files as immutable inputs. Interpolated `dat` files,
meshes, histories, snapshots, figures, and animations are written to
`../output/` (or the directory selected with `--output`).
