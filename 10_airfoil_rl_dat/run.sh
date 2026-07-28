#!/usr/bin/env bash

# Purpose: train or evaluate the pure-airfoil-data Double-DQN example.
# Usage: ./run.sh [TRAIN_DQN_OPTIONS...]
# All arguments are forwarded to train_dqn.py; use --help for the full list.
# Output: checkpoints, configuration, training history, rollout data, PNGs,
# and GIFs below output/.
# Configuration: edit ../course_config.local once for non-default paths.

set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
AFEPACK_EXAMPLES_ROOT="$(cd "$script_dir/.." && pwd)"
source "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
python_executable="${PYTHON:-$PYTHON_BIN}"
course_require_executable Python "$python_executable" PYTHON_BIN

cd "$script_dir"
"$python_executable" train_dqn.py "$@"
