#!/bin/sh

# Purpose: build and run the basic conforming P1 Poisson example.
# Usage: ./run.sh [ROOT_MESH]
# ROOT_MESH defaults to ../common/mesh/unit_square/D.
# Output: output/u.dx; the solver also prints the L2 error.
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

"$script_dir/main" "$root_mesh" u.dx

echo "Wrote $script_dir/output/u.dx"
