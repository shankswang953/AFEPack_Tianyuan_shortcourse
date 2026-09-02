#!/bin/sh

# Purpose: sample a five-parameter airfoil family, reject invalid geometry,
# select a dispersed subset, and generate one EasyMesh grid per selection.
# Usage: ./run.sh [--candidates N] [--selected N] [--seed N]

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
AFEPACK_EXAMPLES_ROOT=$(CDPATH= cd -- "$script_dir/.." && pwd)
. "$AFEPACK_EXAMPLES_ROOT/course_config.sh"

course_require_executable Python "$PYTHON_BIN" PYTHON_BIN
course_require_executable EasyMesh "$EASYMESH_BIN" EASYMESH_BIN
course_require_executable easymesh2mesh "$EASYMESH2MESH_BIN" EASYMESH2MESH_BIN

output_dir="$script_dir/output"
generator_dir="$script_dir/../06_airfoil_mesh_motion"
generator="$generator_dir/generate_airfoil_geometry"

make -C "$generator_dir" generate_airfoil_geometry

"$PYTHON_BIN" "$script_dir/sample_design_space.py" \
  --output "$output_dir" "$@"

while IFS= read -r sample_id; do
  test -n "$sample_id" || continue
  data_file="$output_dir/data/$sample_id.dat"
  case_dir="$output_dir/cases/$sample_id"
  mkdir -p "$case_dir"

  "$generator" \
    "$data_file" "$case_dir" 96 0.35 0.14 0.0 0.0 "$data_file"

  (
    cd "$case_dir"
    "$EASYMESH_BIN" airfoil.d || {
      test -s airfoil.n && test -s airfoil.e && test -s airfoil.s
    }
    "$EASYMESH2MESH_BIN" airfoil airfoil.mesh
  )
done < "$output_dir/selected_ids.txt"

if course_python_can_plot; then
  "$PLOT_PYTHON" "$script_dir/plot_results.py" "$output_dir"
else
  echo "Python plotting is disabled; parameter tables, data, and meshes were generated."
fi

echo "Output directory: $output_dir"
echo "  candidates.csv / accepted.csv / rejected.csv"
echo "  selected.csv and data/*.dat"
echo "  cases/*/airfoil.mesh"
echo "  figures/*.png"
