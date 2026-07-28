#!/bin/sh

# Purpose: run the original two-circle peer-mesh merge and Poisson demo.
# Usage: ./run_circles.sh [ROOT_MESH]
# ROOT_MESH defaults to ../common/mesh/unit_square/D.
# Output: D_left, D_right, D_common mesh files and u_common.dx in output/.
# Configuration: edit ../course_config.local once for non-default paths.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
AFEPACK_EXAMPLES_ROOT=$(CDPATH= cd -- "$script_dir/.." && pwd)
. "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
root_mesh=${1:-"$script_dir/../common/mesh/unit_square/D"}

case "$root_mesh" in
  /*) ;;
  *) root_mesh="$(pwd)/$root_mesh" ;;
esac

make -C "$script_dir"
mkdir -p "$script_dir/output"
cd "$script_dir/output"

"$script_dir/mesh_merge_demo" "$root_mesh" 4

"$script_dir/poisson_solver" D_common u_common.dx

echo "Outputs are in $script_dir/output"
