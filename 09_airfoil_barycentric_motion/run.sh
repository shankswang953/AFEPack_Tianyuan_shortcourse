#!/usr/bin/env bash

# Purpose: reset and run barycentric mesh continuation, then render figures.
# Usage: ./run.sh [BARYCENTRIC_MOTION_OPTIONS...]
# All arguments are forwarded to barycentric_motion.py; run it with --help for
# steps, smoothing, quality, relaxation, retry, and output options.
# Output: continuation data, meshes, snapshots, PNGs, and GIFs under output/.
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
"$python_bin" -B barycentric_motion.py --reset "$@"
MPLCONFIGDIR="$script_dir/output/.matplotlib" \
  XDG_CACHE_HOME="$script_dir/output/.cache" \
  "$python_bin" -B plot_results.py "$script_dir/output"

echo "Barycentric continuation and figures are in $script_dir/output"
