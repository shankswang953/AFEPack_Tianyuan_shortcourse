# `dt_airfoil` Python package

These modules implement geometry/actions, the real mesh environment, replay
serialization, PyTorch models/training, point-cloud diagnostics, and trajectory
rendering. They are imported by the entry scripts in the parent directory and
are not intended to be run directly.

The environment reads `../data/` and maintains generated state in
`../output/`. Mesh tools are selected through the top-level optional
`course_config.local`. Model code requires NumPy and PyTorch; rendering
additionally requires Matplotlib and ImageIO/Pillow support.
