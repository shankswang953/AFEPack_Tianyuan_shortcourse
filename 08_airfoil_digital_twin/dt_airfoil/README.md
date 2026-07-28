# `dt_airfoil` Python package

These modules implement geometry/actions, the real mesh environment, replay
serialization, PyTorch models/training, point-cloud diagnostics, and trajectory
rendering. They are imported by the entry scripts in the parent directory and
are not intended to be run directly.

The environment reads `../data/` and maintains generated state in
`../output/`. Define `EASYMESH_BIN` and `EASYMESH2MESH_BIN` when the mesh
tools are not under `$HOME/bin`. Model code requires NumPy and PyTorch;
rendering additionally requires Matplotlib and ImageIO/Pillow support.
