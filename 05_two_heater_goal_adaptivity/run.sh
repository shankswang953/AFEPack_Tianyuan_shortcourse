#!/bin/sh

# Purpose: compare residual, dual-magnitude, and DWR adaptive refinements.
# Usage: ./run.sh [ROOT_MESH] [--rounds N | --comparison-rounds R [U] D |
#                                   --uniform-reference LEVEL]
# ROOT_MESH defaults to ../common/mesh/unit_square/D.
# Output: output/runs/<configuration> with summary, fields, meshes, and figures.
# Configuration: edit ../course_config.local once for non-default paths.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
AFEPACK_EXAMPLES_ROOT=$(CDPATH= cd -- "$script_dir/.." && pwd)
. "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
root_mesh="$script_dir/../common/mesh/unit_square/D"
if [ "$#" -gt 0 ]; then
  case "$1" in
    --rounds|--comparison-rounds|--uniform-reference)
      ;;
    *)
      root_mesh=$1
      shift
      ;;
  esac
fi

case "$root_mesh" in
  /*) ;;
  *) root_mesh="$(pwd)/$root_mesh" ;;
esac

make -C "$script_dir"

run_name=default
if [ "$#" -eq 2 ] && [ "$1" = "--rounds" ]; then
  run_name="rounds_$2"
fi
if [ "$#" -eq 3 ] && [ "$1" = "--comparison-rounds" ]; then
  run_name="comparison_residual_$2_dual_$3_dwr_$3"
fi
if [ "$#" -eq 4 ] && [ "$1" = "--comparison-rounds" ]; then
  run_name="comparison_residual_$2_dual_$3_dwr_$4"
fi
if [ "$#" -eq 2 ] && [ "$1" = "--uniform-reference" ]; then
  run_name="uniform_reference_level_$2"
fi

output_dir="$script_dir/output/runs/$run_name"
mkdir -p "$output_dir"
cd "$output_dir"

"$script_dir/main" "$root_mesh" "$@"

if [ -f "$output_dir/summary/functional_comparison.dat" ] &&
   course_python_can_plot; then
  plot_cache=${MPLCONFIGDIR:-"${TMPDIR:-/tmp}/afepack-matplotlib"}
  mkdir -p "$plot_cache"
  if ! MPLCONFIGDIR="$plot_cache" \
       "$PLOT_PYTHON" "$script_dir/plot_summary.py" "$output_dir"; then
    echo "Warning: numerical outputs succeeded, but PNG plotting failed." >&2
  fi
fi

echo "Output directory: $output_dir"
echo "  summary/  convergence histories and reference values"
echo "  fields/   solution, residual, dual, and DWR data"
echo "  meshes/   OpenDX meshes for each strategy"
if [ -d "$output_dir/figures" ]; then
  echo "  figures/  optional PNG summaries"
fi
