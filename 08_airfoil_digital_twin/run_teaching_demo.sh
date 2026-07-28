#!/usr/bin/env bash

# Purpose: launch the deterministic, clean, full digital-twin teaching run.
# Usage: ./run_teaching_demo.sh
# The fixed recipe resets state and passes seed/training/controller parameters
# to run_experiment.py; use that Python entry point to change the parameters.
# Output: configuration, replay, models, meshes, histories, and animation under
# output/.
# Configuration: set PYTHON_BIN, EASYMESH_BIN, EASYMESH2MESH_BIN, and the
# documented backend build variables when their defaults are not valid.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${PYTHON_BIN:-}" ]]; then
  python_bin="$PYTHON_BIN"
elif [[ -x "$HOME/anaconda3/bin/python" ]]; then
  python_bin="$HOME/anaconda3/bin/python"
else
  python_bin="python3"
fi

cd "$script_dir"
exec "$python_bin" run_experiment.py \
  --reset \
  --seed 2026 \
  --warmup-episodes 10 \
  --warmup-steps 10 \
  --twin-epochs 1000 \
  --controller-twin-epochs 30 \
  --optimization-episodes 6 \
  --optimization-steps 12 \
  --controller-steps 180 \
  --candidates 512 \
  --max-real-trials 12 \
  --min-improvement 1e-8 \
  --minimum-mesh-quality 0.40 \
  --exploration-bonus 1e-3 \
  --twin-safeguard-ratio 0.25 \
  --fine-action-threshold 0.01
