#!/bin/sh

# Purpose: solve the manufactured boundary-flux problem and its discrete dual.
# Usage: ./run.sh [ROOT_MESH]
# ROOT_MESH defaults to ../common/mesh/unit_square/D.
# Output: primal/dual fields, residual tables, and SVG figures in output/.
# Configuration: set the documented AFEPack build/runtime variables when their
# default $HOME and /opt/local locations do not match the current system.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root_mesh=${1:-"$script_dir/../common/mesh/unit_square/D"}

case "$root_mesh" in
  /*) ;;
  *) root_mesh="$(pwd)/$root_mesh" ;;
esac

make -C "$script_dir"
mkdir -p "$script_dir/output"
cd "$script_dir/output"

AFEPACK_PATH="${AFEPACK_PATH:-$HOME/include/AFEPack}" \
AFEPACK_TEMPLATE_PATH="${AFEPACK_TEMPLATE_PATH:-$HOME/include/AFEPack/template/triangle:$HOME/include/AFEPack/template/twin_triangle}" \
  "$script_dir/main" "$root_mesh"

echo "Outputs are in $script_dir/output"
