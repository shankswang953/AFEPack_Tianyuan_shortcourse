#!/bin/sh

# Purpose: generate meshes for all checked-in EasyMesh teaching geometries.
# Usage: ./run_examples.sh [OUTPUT_DIR]
# OUTPUT_DIR defaults to this directory's output/ folder.
# Input: every *.d file beside this script.
# Output: copied .d files and generated .n, .e, and .s files in OUTPUT_DIR.
# Configuration: edit ../../../course_config.local once for non-default paths.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
AFEPACK_EXAMPLES_ROOT=$(CDPATH= cd -- "$script_dir/../../.." && pwd)
. "$AFEPACK_EXAMPLES_ROOT/course_config.sh"
output_dir=${1:-"$script_dir/output"}

case "$output_dir" in
  /*) ;;
  *) output_dir="$(pwd)/$output_dir" ;;
esac

course_require_executable EasyMesh "$EASYMESH_BIN" EASYMESH_BIN

mkdir -p "$output_dir"

for source_file in "$script_dir"/*.d; do
  example_name=$(basename "$source_file" .d)
  cp "$source_file" "$output_dir/$example_name.d"
  (
    cd "$output_dir"
    "$EASYMESH_BIN" "$example_name.d" -g -m || {
      # EasyMesh 1.4 may return 1 after producing a valid grid.
      test -s "$example_name.n"
      test -s "$example_name.e"
      test -s "$example_name.s"
    }
  )
  echo "Generated $example_name.[nes] in $output_dir"
done
