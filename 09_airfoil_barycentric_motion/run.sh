#!/usr/bin/env bash

# Purpose: reset and run barycentric mesh continuation, then render figures.
# Usage: ./run.sh [BARYCENTRIC_MOTION_OPTIONS...]
# All arguments are forwarded to barycentric_motion.py; run it with --help for
# steps, smoothing, quality, relaxation, retry, and output options.
# Output: continuation data, meshes, snapshots, PNGs, and GIFs under output/.
# Configuration: edit ../course_config.local once for non-default paths.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AFEPACK_EXAMPLES_ROOT="$(cd "$script_dir/.." && pwd)"
source "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
course_require_executable Python "$PYTHON_BIN" PYTHON_BIN
course_require_executable EasyMesh "$EASYMESH_BIN" EASYMESH_BIN
course_require_executable easymesh2mesh "$EASYMESH2MESH_BIN" EASYMESH2MESH_BIN

cd "$script_dir"
"$PYTHON_BIN" -B barycentric_motion.py --reset "$@"
if course_python_can_plot; then
  MPLCONFIGDIR="$script_dir/output/.matplotlib" \
    XDG_CACHE_HOME="$script_dir/output/.cache" \
    "$PLOT_PYTHON" -B plot_results.py "$script_dir/output"
else
  echo "Python plotting is unavailable or disabled; continuation data were generated." >&2
fi

echo "Barycentric continuation output is in $script_dir/output"
