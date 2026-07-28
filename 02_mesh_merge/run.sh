#!/bin/sh

# Purpose: refine Tian/Yuan peer meshes, merge them, and render the result.
# Usage: ./run.sh [ROUNDS] [STROKE_HALF_WIDTH]
# Defaults: 3 rounds and a 0.035 stroke half-width.
# Output: mesh/CSV data in output/ and PNG files in output/figures/.
# Configuration: edit ../course_config.local once for non-default paths.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
AFEPACK_EXAMPLES_ROOT=$(CDPATH= cd -- "$script_dir/.." && pwd)
. "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
rounds=${1:-3}
stroke_half_width=${2:-0.035}

course_require_executable EasyMesh "$EASYMESH_BIN" EASYMESH_BIN

make -C "$script_dir"
mkdir -p "$script_dir/output"
cp "$script_dir/tianyuan_rectangle.d" "$script_dir/output/T_root.d"

(
  cd "$script_dir/output"
  "$EASYMESH_BIN" T_root.d || {
    # EasyMesh commonly returns status 1 even after writing a valid grid.
    test -s T_root.n && test -s T_root.e && test -s T_root.s
  }
  "$script_dir/tianyuan_mesh_merge" \
    T_root "$rounds" "$stroke_half_width"
)

if course_python_can_plot; then
  "$PLOT_PYTHON" "$script_dir/plot_tianyuan_merge.py" \
    "$script_dir/output"
else
  echo "Python plotting is unavailable or disabled; mesh data were generated." >&2
fi

echo "Mesh data: $script_dir/output"
echo "Figures:   $script_dir/output/figures"
