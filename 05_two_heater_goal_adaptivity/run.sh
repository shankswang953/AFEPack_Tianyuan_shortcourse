#!/bin/sh

# Purpose: compare residual, dual-magnitude, and DWR adaptive refinements.
# Usage: ./run.sh [ROOT_MESH] [--rounds N | --comparison-rounds R [U] D]
# ROOT_MESH defaults to ../common/mesh/unit_square/D.
# Output: output/ by default; round/comparison modes use named subdirectories.
# Configuration: edit ../course_config.local once for non-default paths.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
AFEPACK_EXAMPLES_ROOT=$(CDPATH= cd -- "$script_dir/.." && pwd)
. "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
root_mesh=${1:-"$script_dir/../common/mesh/unit_square/D"}
if [ "$#" -gt 0 ]; then
  shift
fi

case "$root_mesh" in
  /*) ;;
  *) root_mesh="$(pwd)/$root_mesh" ;;
esac

make -C "$script_dir"
output_dir="$script_dir/output"
if [ "$#" -eq 2 ] && [ "$1" = "--rounds" ] && [ "$2" != "6" ]; then
  output_dir="$script_dir/output/rounds_$2"
fi
if [ "$#" -eq 3 ] && [ "$1" = "--comparison-rounds" ]; then
  output_dir="$script_dir/output/comparison_residual_$2_dwr_$3"
fi
if [ "$#" -eq 4 ] && [ "$1" = "--comparison-rounds" ]; then
  output_dir="$script_dir/output/comparison_residual_$2_dual_$3_dwr_$4"
fi
mkdir -p "$output_dir"
cd "$output_dir"

"$script_dir/main" "$root_mesh" "$@"

echo "Outputs are in $output_dir"
