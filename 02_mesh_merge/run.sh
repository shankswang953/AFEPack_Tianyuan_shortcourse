#!/bin/sh

# Purpose: refine Tian/Yuan peer meshes, merge them, and render the result.
# Usage: ./run.sh [ROUNDS] [STROKE_HALF_WIDTH]
# Defaults: 3 rounds and a 0.035 stroke half-width.
# Output: mesh/CSV data in output/ and PNG files in output/figures/.
# Configuration: set EASYMESH_BIN when EasyMesh is not at $HOME/bin/easymesh;
# set PLOT_PYTHON and the documented build variables when needed.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
rounds=${1:-3}
stroke_half_width=${2:-0.035}
easymesh_bin=${EASYMESH_BIN:-"$HOME/bin/easymesh"}

if [ ! -x "$easymesh_bin" ]; then
  echo "EasyMesh executable not found: $easymesh_bin" >&2
  echo "Set EASYMESH_BIN to the executable location." >&2
  exit 1
fi

make -C "$script_dir"
mkdir -p "$script_dir/output"
cp "$script_dir/tianyuan_rectangle.d" "$script_dir/output/T_root.d"

(
  cd "$script_dir/output"
  "$easymesh_bin" T_root.d || {
    # EasyMesh commonly returns status 1 even after writing a valid grid.
    test -s T_root.n && test -s T_root.e && test -s T_root.s
  }
  "$script_dir/tianyuan_mesh_merge" \
    T_root "$rounds" "$stroke_half_width"
)

plot_python=${PLOT_PYTHON:-}
if [ -z "$plot_python" ]; then
  for candidate in "$HOME/anaconda3/bin/python" python3 python; do
    if command -v "$candidate" >/dev/null 2>&1; then
      plot_python=$candidate
      break
    fi
  done
fi

if [ -n "$plot_python" ]; then
  "$plot_python" "$script_dir/plot_tianyuan_merge.py" \
    "$script_dir/output"
else
  echo "Python was not found; mesh data were generated without PNG figures." >&2
fi

echo "Mesh data: $script_dir/output"
echo "Figures:   $script_dir/output/figures"
