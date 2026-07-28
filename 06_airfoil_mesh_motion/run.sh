#!/bin/sh

# Purpose: generate, move, validate, persist, and plot an AFEPack airfoil mesh.
# Usage: ./run.sh [DATA_FILE] [OUTPUT_DIR] [MOVED_DATA_FILE]
#                 [SMOOTH_ITERATIONS] [move|smooth-only]
# Defaults: data/naca0012.dat, output/, no moved data, 50 iterations, and move.
# Output: persistent/current meshes, CSV/DX diagnostics, and figures in OUTPUT_DIR.
# Configuration: define EASYMESH_BIN, EASYMESH2MESH_BIN, PLOT_PYTHON, and the
# documented build variables when their default locations are not valid.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
data_file=${1:-"$script_dir/data/naca0012.dat"}
output_dir=${2:-"$script_dir/output"}
moved_data_file=${3:-}
smoothing_iterations=${4:-50}
mode=${5:-move}
current_mesh="$output_dir/mesh_current.mesh"
current_boundary="$output_dir/boundary_current.dat"

easymesh_bin=${EASYMESH_BIN:-"$HOME/bin/easymesh"}
easymesh2mesh_bin=${EASYMESH2MESH_BIN:-"$HOME/bin/easymesh2mesh"}

make -C "$script_dir"
mkdir -p "$output_dir"

initial_mesh_created=0
case "$mode" in
  smooth-only)
    if ! test -s "$current_mesh" || ! test -s "$current_boundary"; then
      echo "No persistent mesh is available for smoothing." >&2
      echo "Run one geometry action first." >&2
      exit 2
    fi
    # The current boundary is both the source and target: only interior
    # boundary-mark-0 vertices will move during Laplacian smoothing.
    cp "$current_boundary" "$output_dir/boundary_initial.dat"
    cp "$current_boundary" "$output_dir/boundary_moved.dat"
    # Keep the established figure layout while showing unchanged geometry.
    if test -s "$output_dir/fit_moved.csv"; then
      cp "$output_dir/fit_moved.csv" "$output_dir/fit_initial.csv"
    fi
    update_kind=fixed_topology_smoothing_only
    ;;
  move)
    if test -n "$moved_data_file"; then
      "$script_dir/generate_airfoil_geometry" \
        "$data_file" "$output_dir" 96 0.35 0.14 0.060 0.015 \
        "$moved_data_file"
    else
      "$script_dir/generate_airfoil_geometry" \
        "$data_file" "$output_dir" 96 0.35 0.14 0.060 0.015
    fi

    # The fixed-topology workflow needs EasyMesh only once. A standalone
    # built-in test restarts from its prescribed initial geometry.
    reset_mesh=${RESET_MESH:-0}
    if test -z "$moved_data_file"; then
      reset_mesh=1
    fi
    if test "$reset_mesh" = "1"; then
      rm -f "$current_mesh" "$current_boundary"
    fi

    if ! test -s "$current_mesh" || ! test -s "$current_boundary"; then
      (
        cd "$output_dir"
        "$easymesh_bin" airfoil.d || {
          test -s airfoil.n
          test -s airfoil.e
          test -s airfoil.s
        }
        "$easymesh2mesh_bin" airfoil airfoil.mesh
      )
      cp "$output_dir/airfoil.mesh" "$current_mesh"
      cp "$output_dir/boundary_initial.dat" "$current_boundary"
      initial_mesh_created=1
    elif ! cmp -s \
        "$current_boundary" "$output_dir/boundary_initial.dat"; then
      echo "The persistent mesh boundary does not match the current airfoil." >&2
      echo "Run 'python3 airfoil_step.py --reset' before continuing." >&2
      exit 2
    fi
    update_kind=fixed_topology_boundary_motion
    ;;
  *)
    echo "Unknown run mode: $mode" >&2
    exit 2
    ;;
esac

"$script_dir/move_and_smooth" \
  "$current_mesh" \
  "$current_boundary" \
  "$output_dir/boundary_moved.dat" \
  "$output_dir" \
  "$smoothing_iterations" 0.45

# Commit the mesh only after movement, smoothing, and validity checks succeed.
cp "$output_dir/mesh_smoothed.mesh" "$current_mesh.tmp"
mv "$current_mesh.tmp" "$current_mesh"
if test "$mode" = "move"; then
  cp "$output_dir/boundary_moved.dat" "$current_boundary.tmp"
  mv "$current_boundary.tmp" "$current_boundary"
fi

{
  echo "update=$update_kind"
  echo "initial_mesh_created=$initial_mesh_created"
  echo "smoothing_iterations=$smoothing_iterations"
  echo "persistent_mesh=$current_mesh"
} > "$output_dir/mesh_update_mode.txt"

plot_python=${PLOT_PYTHON:-}
if test -z "$plot_python"; then
  for candidate in python3 "$HOME/anaconda3/bin/python"; do
    if command -v "$candidate" >/dev/null 2>&1 \
        && "$candidate" -c "import matplotlib" >/dev/null 2>&1; then
      plot_python=$candidate
      break
    fi
  done
fi

if test -n "$plot_python"; then
  mkdir -p "$output_dir/.matplotlib"
  MPLCONFIGDIR="$output_dir/.matplotlib" \
    XDG_CACHE_HOME="$output_dir/.cache" \
    "$plot_python" "$script_dir/visualize_results.py" "$output_dir"
else
  echo "Matplotlib was not found; CSV and OpenDX output are still available."
  echo "Set PLOT_PYTHON to a Python executable with matplotlib to create PNG files."
fi

if test "$mode" = "smooth-only"; then
  echo "Smoothed the persistent AFEPack mesh; geometry and topology are unchanged."
elif test "$initial_mesh_created" = "1"; then
  echo "EasyMesh created the initial mesh."
else
  echo "Reused the persistent AFEPack mesh; EasyMesh was not called."
fi
echo "Outputs are in $output_dir"
