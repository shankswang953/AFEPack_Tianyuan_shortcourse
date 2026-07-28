#!/usr/bin/env bash

# Purpose: launch the deterministic, clean, full digital-twin teaching run.
# Usage: ./run_teaching_demo.sh
# The fixed recipe resets state and passes seed/training/controller parameters
# to run_experiment.py; use that Python entry point to change the parameters.
# Output: configuration, replay, models, meshes, histories, and animation under
# output/.
# Configuration: edit ../course_config.local once for non-default paths.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AFEPACK_EXAMPLES_ROOT="$(cd "$script_dir/.." && pwd)"
source "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
course_require_executable Python "$PYTHON_BIN" PYTHON_BIN
course_require_executable EasyMesh "$EASYMESH_BIN" EASYMESH_BIN
course_require_executable easymesh2mesh "$EASYMESH2MESH_BIN" EASYMESH2MESH_BIN

cd "$script_dir"
exec "$PYTHON_BIN" run_experiment.py \
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
