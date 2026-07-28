#!/usr/bin/env bash

# Purpose: train or evaluate the pure-airfoil-data Double-DQN example.
# Usage: ./run.sh [TRAIN_DQN_OPTIONS...]
# All arguments are forwarded to train_dqn.py; use --help for the full list.
# Output: checkpoints, configuration, training history, rollout data, PNGs,
# and GIFs below output/.
# Configuration: set PYTHON to an interpreter with NumPy, PyTorch, Matplotlib,
# and ImageIO/Pillow when python3 is not the correct environment.

set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
python_executable="${PYTHON:-python3}"

cd "$script_dir"
"$python_executable" train_dqn.py "$@"
