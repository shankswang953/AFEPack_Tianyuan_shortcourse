# README animations

This directory contains selected high-resolution result animations embedded in
the repository README files. Unlike reproducible files below each example's
ignored `output/` directory, these presentation copies are tracked by Git.

| animation | source example | resolution | regeneration command |
|---|---|---:|---|
| `digital_twin_shape_evolution.gif` | `08_airfoil_digital_twin` | 1600 x 1200 | `python3 render_trajectory.py --fps 2 --gif-dpi 200` |
| `barycentric_fixed_topology.gif` | `09_airfoil_barycentric_motion` | 1380 x 750 | `python3 plot_results.py output --gif-dpi 150` |
| `barycentric_smoothing_mechanism.gif` | `09_airfoil_barycentric_motion` | 2325 x 660 | `python3 plot_results.py output --gif-dpi 150` |
| `rl_shape_evolution.gif` | `10_airfoil_rl_dat` | 1440 x 972 | `python3 train_dqn.py --evaluate-only --max-steps 140 --seed 2026 --gif-dpi 180` |

After regeneration, copy the four files from their documented `output/`
locations into this directory. The two barycentric GIFs checked in here use
GIF frame-delta optimization to reduce repository size without changing their
pixel dimensions.
