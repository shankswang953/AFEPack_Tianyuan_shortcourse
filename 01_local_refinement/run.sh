#!/bin/sh

# Purpose: run the independent/combined local-indicator refinement study.
# Usage: ./run.sh [ROOT_MESH] [ROUNDS] [LOCAL_FRACTION] [SCALE_B]
# Defaults: ../common/mesh/unit_square/D, 3, 0.05, and 50.
# Output: meshes and CSV files in output/, plus PNG files in output/figures/.
# Configuration: edit ../course_config.local once for non-default paths.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
AFEPACK_EXAMPLES_ROOT=$(CDPATH= cd -- "$script_dir/.." && pwd)
. "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
root_mesh=${1:-"$script_dir/../common/mesh/unit_square/D"}
rounds=${2:-3}
local_fraction=${3:-0.05}
scale_b=${4:-50}

case "$root_mesh" in
  /*) ;;
  *) root_mesh="$(pwd)/$root_mesh" ;;
esac

make -C "$script_dir"
mkdir -p "$script_dir/output"

"$script_dir/indicator_study" \
  "$root_mesh" \
  "$script_dir/output" \
  "$rounds" \
  "$local_fraction" \
  "$scale_b"

if course_python_can_plot; then
  MPLCONFIGDIR="$script_dir/output/.matplotlib" \
    "$PLOT_PYTHON" "$script_dir/plot_indicator_study.py" \
    "$script_dir/output" --scale-b "$scale_b"
else
  echo "Python plotting is unavailable or disabled; CSV/mesh data were generated." >&2
fi

echo "Outputs are in $script_dir/output"
