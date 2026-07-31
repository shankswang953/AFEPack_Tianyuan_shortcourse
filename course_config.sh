#!/bin/sh

# Shared runtime configuration for all teaching examples.
#
# The calling run.sh sets AFEPACK_EXAMPLES_ROOT before sourcing this file.
# A machine-specific course_config.local is optional.  Copy
# course_config.local.example to course_config.local and edit only the paths
# that differ on the local machine.

if [ -z "${AFEPACK_EXAMPLES_ROOT:-}" ]; then
  echo "AFEPACK_EXAMPLES_ROOT must be set before sourcing course_config.sh." >&2
  return 2 2>/dev/null || exit 2
fi

course_local_config="$AFEPACK_EXAMPLES_ROOT/course_config.local"
if [ -f "$course_local_config" ]; then
  # shellcheck disable=SC1090
  . "$course_local_config"
fi

# AFEPack: finite-element assembly, solves, refinement, and mesh operations.
: "${AFEPACK_PREFIX:=$HOME}"
: "${OPENBLAS_PREFIX:=/opt/local}"
: "${BOOST_INCLUDE:=/opt/local/libexec/boost/1.81/include}"
: "${AFEPACK_PATH:=$AFEPACK_PREFIX/include/AFEPack}"
: "${AFEPACK_TEMPLATE_PATH:=$AFEPACK_PATH/template/triangle:$AFEPACK_PATH/template/twin_triangle}"

course_find_executable()
{
  preferred=$1
  fallback_name=$2

  if [ -n "$preferred" ] && command -v "$preferred" >/dev/null 2>&1; then
    command -v "$preferred"
    return 0
  fi
  if command -v "$fallback_name" >/dev/null 2>&1; then
    command -v "$fallback_name"
    return 0
  fi
  return 1
}

# EasyMesh: initial triangulation and topology reconstruction after remeshing.
if [ -z "${EASYMESH_BIN:-}" ]; then
  EASYMESH_BIN=$(course_find_executable "$HOME/bin/easymesh" easymesh || true)
fi
if [ -z "${EASYMESH2MESH_BIN:-}" ]; then
  EASYMESH2MESH_BIN=$(course_find_executable "$HOME/bin/easymesh2mesh" easymesh2mesh || true)
fi

# Python: optional plotting/helper layer.  C++ AFEPack/EasyMesh examples still
# run when it is unavailable.  Examples 08--10 use Python for their algorithms
# and therefore check it explicitly in their own launchers.
if [ -z "${PYTHON_BIN:-}" ]; then
  PYTHON_BIN=$(course_find_executable "$HOME/anaconda3/bin/python" python3 || true)
  if [ -z "$PYTHON_BIN" ]; then
    PYTHON_BIN=$(course_find_executable "" python || true)
  fi
fi
: "${ENABLE_PYTHON_PLOTS:=auto}"
case "$ENABLE_PYTHON_PLOTS" in
  0|no|false|off)
    PLOT_PYTHON=
    ;;
  *)
    : "${PLOT_PYTHON:=$PYTHON_BIN}"
    ;;
esac
: "${PYTHON:=$PYTHON_BIN}"

# ImageMagick: rasterize the C++-generated teaching figures when an example
# publishes PNG output without requiring Python.
if [ -z "${IMAGE_CONVERT_BIN:-}" ]; then
  IMAGE_CONVERT_BIN=$(course_find_executable \
    "/opt/local/bin/convert" convert || true)
fi

course_require_executable()
{
  tool_name=$1
  executable=$2
  variable_name=$3

  if [ -z "$executable" ] || ! command -v "$executable" >/dev/null 2>&1; then
    echo "$tool_name executable was not found." >&2
    echo "Set $variable_name in $AFEPACK_EXAMPLES_ROOT/course_config.local." >&2
    return 1
  fi
}

course_python_can_plot()
{
  [ -n "${PLOT_PYTHON:-}" ] &&
    command -v "$PLOT_PYTHON" >/dev/null 2>&1 &&
    "$PLOT_PYTHON" -c "import matplotlib, numpy" >/dev/null 2>&1
}

export AFEPACK_PREFIX OPENBLAS_PREFIX BOOST_INCLUDE
export AFEPACK_PATH AFEPACK_TEMPLATE_PATH
export EASYMESH_BIN EASYMESH2MESH_BIN
export PYTHON_BIN PLOT_PYTHON PYTHON ENABLE_PYTHON_PLOTS
export IMAGE_CONVERT_BIN
