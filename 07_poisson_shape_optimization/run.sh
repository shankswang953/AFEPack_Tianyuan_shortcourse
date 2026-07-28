#!/bin/sh

# Purpose: build and run the fixed-area Poisson shape optimization.
# Usage: ./run.sh [OUTPUT_DIR] [SWEEPS]
# Defaults: output/ and 4 coordinate-search sweeps.
# Output: optimization history, meshes, fields, and per-iteration data in
# OUTPUT_DIR.
# Configuration: edit ../course_config.local once for non-default paths.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
AFEPACK_EXAMPLES_ROOT=$(CDPATH= cd -- "$script_dir/.." && pwd)
. "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
output_dir=${1:-"$script_dir/output"}
sweeps=${2:-4}

course_require_executable EasyMesh "$EASYMESH_BIN" EASYMESH_BIN
make -C "$script_dir"

"$script_dir/shape_optimization" "$output_dir" "$sweeps"
