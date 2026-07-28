#!/bin/sh

# Purpose: generate meshes for all checked-in EasyMesh teaching geometries.
# Usage: ./run_examples.sh [OUTPUT_DIR]
# OUTPUT_DIR defaults to this directory's output/ folder.
# Input: every *.d file beside this script.
# Output: copied .d files and generated .n, .e, and .s files in OUTPUT_DIR.
# Configuration: set EASYMESH_BIN if EasyMesh is not at $HOME/bin/easymesh.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
output_dir=${1:-"$script_dir/output"}
easymesh_bin=${EASYMESH_BIN:-"$HOME/bin/easymesh"}

case "$output_dir" in
  /*) ;;
  *) output_dir="$(pwd)/$output_dir" ;;
esac

if [ ! -x "$easymesh_bin" ]; then
  echo "EasyMesh executable not found: $easymesh_bin" >&2
  echo "Set EASYMESH_BIN to the executable location." >&2
  exit 1
fi

mkdir -p "$output_dir"

for source_file in "$script_dir"/*.d; do
  example_name=$(basename "$source_file" .d)
  cp "$source_file" "$output_dir/$example_name.d"
  (
    cd "$output_dir"
    "$easymesh_bin" "$example_name.d" -g -m || {
      # EasyMesh 1.4 may return 1 after producing a valid grid.
      test -s "$example_name.n"
      test -s "$example_name.e"
      test -s "$example_name.s"
    }
  )
  echo "Generated $example_name.[nes] in $output_dir"
done
