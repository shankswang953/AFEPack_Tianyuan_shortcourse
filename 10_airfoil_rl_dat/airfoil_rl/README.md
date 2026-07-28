# `airfoil_rl` Python package

This import-only package contains:

- `geometry.py`: `dat` I/O and Gaussian actions;
- `environment.py`: deterministic shape-matching environment;
- `dqn.py`: Double-DQN agent, replay buffer, and checkpoints;
- `plotting.py`: CSV, `dat`, PNG, and GIF output helpers.

Run `../train_dqn.py` or `../run.sh` for the complete workflow. The package
has no machine-specific AFEPack or EasyMesh paths.
