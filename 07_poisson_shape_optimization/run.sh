#!/bin/sh

# Purpose: build and run the fixed-area Poisson shape optimization.
# Usage: ./run.sh [OUTPUT_DIR] [SWEEPS]
# Defaults: output/ and 4 coordinate-search sweeps.
# Output: optimization history, meshes, fields, and per-iteration data in
# OUTPUT_DIR.
# Configuration: define EASYMESH_BIN and the documented AFEPack build/runtime
# variables when their default $HOME and /opt/local locations are not valid.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
output_dir=${1:-"$script_dir/output"}
sweeps=${2:-4}

make -C "$script_dir"

EASYMESH_BIN="${EASYMESH_BIN:-$HOME/bin/easymesh}" \
AFEPACK_TEMPLATE_PATH="${AFEPACK_TEMPLATE_PATH:-$HOME/include/AFEPack/template/triangle}" \
  "$script_dir/shape_optimization" "$output_dir" "$sweeps"
