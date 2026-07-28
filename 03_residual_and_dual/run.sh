#!/bin/sh

# Purpose: run the residual, discrete-dual, and DWR refinement example.
# Usage: ./run.sh [ROOT_MESH] [PROGRAM_ARGUMENTS...]
# ROOT_MESH defaults to ../common/mesh/unit_square/D; remaining arguments are
# passed to main (for example --uniform-reference LEVEL or source/sensor data).
# Output: all solver tables, fields, meshes, and figures are written to output/.
# Configuration: set the documented AFEPack build/runtime variables as needed.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root_mesh=${1:-"$script_dir/../common/mesh/unit_square/D"}
if [ "$#" -gt 0 ]; then
  shift
fi

case "$root_mesh" in
  /*) ;;
  *) root_mesh="$(pwd)/$root_mesh" ;;
esac

make -C "$script_dir"
mkdir -p "$script_dir/output"
cd "$script_dir/output"

AFEPACK_PATH="${AFEPACK_PATH:-$HOME/include/AFEPack}" \
AFEPACK_TEMPLATE_PATH="${AFEPACK_TEMPLATE_PATH:-$HOME/include/AFEPack/template/triangle:$HOME/include/AFEPack/template/twin_triangle}" \
  "$script_dir/main" "$root_mesh" "$@"

echo "Outputs are in $script_dir/output"
