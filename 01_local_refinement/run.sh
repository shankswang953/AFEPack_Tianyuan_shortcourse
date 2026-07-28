#!/bin/sh

# Purpose: run the independent/combined local-indicator refinement study.
# Usage: ./run.sh [ROOT_MESH] [ROUNDS] [LOCAL_FRACTION] [SCALE_B]
# Defaults: ../common/mesh/unit_square/D, 3, 0.05, and 50.
# Output: meshes and CSV files in output/, plus PNG files in output/figures/.
# Configuration: set AFEPACK_PREFIX/OPENBLAS_PREFIX/BOOST_INCLUDE for the
# build and PLOT_PYTHON for a Python interpreter with Matplotlib and NumPy.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
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

if [ -n "${PLOT_PYTHON:-}" ]; then
  plot_python=$PLOT_PYTHON
elif [ -x "$HOME/anaconda3/bin/python" ]; then
  plot_python=$HOME/anaconda3/bin/python
else
  plot_python=python3
fi

MPLCONFIGDIR="$script_dir/output/.matplotlib" \
  "$plot_python" "$script_dir/plot_indicator_study.py" \
  "$script_dir/output" --scale-b "$scale_b"

echo "Outputs are in $script_dir/output"
