# AFEPack course examples

**High-resolution course slides:** [View or download from Google Drive](https://drive.google.com/file/d/1hIARgAdCFNDT0fIGOSmqJkdJqEACcKSq/view?usp=sharing)

The slides introduce the concepts behind the examples below and include the
commands recommended for live classroom demonstrations.

<!-- script-interface -->
## First-time setup: configure paths once

If the built-in macOS/MacPorts defaults already match your machine, no setup
is needed. Otherwise, before running any example, create one machine-local
configuration file at the repository root:

```bash
cp course_config.local.example course_config.local
```

Open `course_config.local` and uncomment only the paths that differ on your
machine. For example:

```text
AFEPACK_PREFIX=/opt/afepack
OPENBLAS_PREFIX=/opt/openblas
BOOST_INCLUDE=/opt/boost/include
EASYMESH_BIN=/opt/easymesh/bin/easymesh
EASYMESH2MESH_BIN=/opt/easymesh/bin/easymesh2mesh
PYTHON_BIN=/path/to/python
```

Keep the `NAME=value` form with no spaces around `=`. Do not edit
`course_config.sh`, `course_config.mk`, or `course_config.py`: they are the
shared loaders. `course_config.local` is ignored by Git, while
`course_config.local.example` remains the portable template distributed with
the course.

Every official `run.sh` and the Python mesh entry points automatically read
`course_config.local`; every C++ Makefile reads the same file through
`course_config.mk`.

| layer | role | main settings |
|---|---|---|
| AFEPack computation | finite-element assembly, solves, refinement, and mesh operations | `AFEPACK_PREFIX`, `AFEPACK_PATH`, `AFEPACK_TEMPLATE_PATH`, `OPENBLAS_PREFIX`, `BOOST_INCLUDE` |
| EasyMesh generation | initial triangulation and topology reconstruction after remeshing | `EASYMESH_BIN`, `EASYMESH2MESH_BIN` |
| Python helper/visualization | figures, animations, and the Python-based examples 08--10 | `PYTHON_BIN`, `PLOT_PYTHON`, `ENABLE_PYTHON_PLOTS` |

### Install the Python packages

Python 3.10 or newer is recommended. After setting `PYTHON_BIN` in
`course_config.local`, install every package used in the course with:

```bash
AFEPACK_EXAMPLES_ROOT="$PWD"
. ./course_config.sh
"$PYTHON_BIN" -m pip install -r requirements.txt
```

The packages are grouped by purpose:

| use | required packages |
|---|---|
| numerical arrays and indicator plots | `numpy` |
| PNG figures and GIF animations | `matplotlib`, `pillow` |
| digital-twin neural networks (`08`) | `torch` |
| reinforcement learning and its animation (`10`) | `torch`, `imageio` |

For a smaller installation, example `08` and example `10` also provide their
own requirement files:

```bash
AFEPACK_EXAMPLES_ROOT="$PWD"
. ./course_config.sh
"$PYTHON_BIN" -m pip install -r 08_airfoil_digital_twin/requirements.txt
"$PYTHON_BIN" -m pip install -r 10_airfoil_rl_dat/requirements.txt
```

Python plotting is optional for the C++ examples. Set
`ENABLE_PYTHON_PLOTS=0` to keep numerical and mesh output while skipping
PNG/GIF generation. Python remains required for the algorithms in examples
08--10.

| variable | default | define when |
|---|---|---|
| `CXX` | Make's C++ compiler | a different C++20 compiler is required |
| `AFEPACK_PREFIX` | `$HOME` | AFEPack headers/libraries use another prefix |
| `OPENBLAS_PREFIX` | `/opt/local` | OpenBLAS is not installed by MacPorts |
| `BOOST_INCLUDE` | `/opt/local/libexec/boost/1.81/include` | Boost headers are elsewhere |
| `AFEPACK_PATH` | `$HOME/include/AFEPack` | AFEPack runtime data are elsewhere |
| `AFEPACK_TEMPLATE_PATH` | under `$HOME/include/AFEPack/template` | AFEPack templates are elsewhere |
| `EASYMESH_BIN` | `$HOME/bin/easymesh`, then `PATH` | EasyMesh is elsewhere |
| `EASYMESH2MESH_BIN` | `$HOME/bin/easymesh2mesh`, then `PATH` | the converter is elsewhere |
| `PLOT_PYTHON` | derived from `PYTHON_BIN` | plotting packages are in another interpreter |
| `PYTHON_BIN` | `$HOME/anaconda3/bin/python`, then `PATH` | another Python environment is needed |
| `ENABLE_PYTHON_PLOTS` | `auto` | use `0` to skip Python figures |

Root-mesh paths are command-line parameters, not machine constants. Keep the
provided `common/mesh/unit_square/D` basename or pass another EasyMesh basename
without the `.n`, `.e`, or `.s` suffix.

The examples are ordered by the concepts introduced in the short course:

1. `00_poisson_basic`: assemble and solve a conforming P1 Poisson problem.
2. `01_local_refinement`: compare independent and combined local indicators,
   and show how a scale coefficient can dominate a shared marking budget.
3. `02_mesh_merge`: independently refine the thick-stroke distance-zero sets
   for Tian and Yuan on a rectangle, then merge the two refinement histories
   into one Tianyuan mesh.
4. `03_residual_and_dual`: visualize the strong residual, perform
   residual-based primal refinement, solve a fully discrete dual on an
   independent mesh, and refine by dual magnitude.
5. `04_boundary_flux_dual`: place the heat source above a local cooled-boundary
   flux target and visualize the corresponding fully discrete dual.
6. `05_two_heater_goal_adaptivity`: compare residual and adjoint-weighted
   refinement for a smooth local-temperature functional with a spectral
   reference value.
7. `06_airfoil_mesh_motion`: fit airfoil data with regularized Bezier curves,
   generate an EasyMesh triangulation, move the airfoil boundary, and smooth
   the interior AFEPack mesh without changing connectivity.
8. `07_poisson_shape_optimization`: optimize a fixed-area Fourier boundary for
   a steady Poisson heat problem; a C++ coordinate search transforms an
   irregular multi-lobed shape into the optimal circle.
9. `08_airfoil_digital_twin`: learn and validate local upper/lower airfoil
   updates, move a persistent mesh, and remesh only when its quality becomes
   unacceptable.
10. `09_airfoil_barycentric_motion`: construct a constant-speed pointwise
    barycentric path from a circular profile to NACA0012, then move and smooth
    one fixed-topology mesh along that path without remeshing.
11. `10_airfoil_rl_dat`: train and evaluate a deterministic Double-DQN
    controller entirely in airfoil `dat` space, without calling AFEPack,
    EasyMesh, smoothing, or a CFD solver.

Examples `00`--`05` use the EasyMesh input
`common/mesh/unit_square/D.[nse]`. Examples `06`--`10` generate their own
geometry and mesh files. Generated outputs stay inside each example's own
`output/` directory.

## Selected result animations

The animations below are high-resolution presentation copies. Click an
animation to open the corresponding example and its reproducible workflow.

| digital-twin boundary control (`08`) | digital-twin mesh evolution (`08`) |
|---|---|
| [![A circular obstacle is progressively controlled toward a NACA0012 target by the online digital twin.](assets/animations/digital_twin_shape_evolution.gif)](08_airfoil_digital_twin/README.md) | [![The real AFEPack triangular mesh follows accepted digital-twin boundary updates and remeshing events.](assets/animations/digital_twin_mesh_evolution.gif)](08_airfoil_digital_twin/README.md) |

| fixed-topology barycentric motion (`09`) | quality-aware boundary smoothing (`09`) |
|---|---|
| [![One fixed-connectivity triangular mesh moves from a circle to a NACA0012 airfoil.](assets/animations/barycentric_fixed_topology.gif)](09_airfoil_barycentric_motion/README.md) | [![Boundary motion, equilateral targets, and quality-aware smoothing are compared during continuation.](assets/animations/barycentric_smoothing_mechanism.gif)](09_airfoil_barycentric_motion/README.md) |

| pure-dat reinforcement learning (`10`) |
|---|
| [![A greedy Double-DQN rollout changes a circular profile into a NACA0012 airfoil.](assets/animations/rl_shape_evolution.gif)](10_airfoil_rl_dat/README.md) |

## Verified lecture commands

Run each command from the named example directory. The commands below were
retested on 2026-07-27.

| example | recommended command | purpose |
|---|---|---|
| `00_poisson_basic` | `./run.sh` | short AFEPack installation check |
| `01_local_refinement` | `./run.sh` | compare independent, equal-sum, and scale-dominated refinement |
| `02_mesh_merge` | `./run.sh 3 0.035` | build the Tian and Yuan peer meshes and merge them |
| `03_residual_and_dual` | `./run.sh` | residual, discrete dual, and one DWR step |
| `04_boundary_flux_dual` | `./run.sh` | manufactured boundary-flux target |
| `05_two_heater_goal_adaptivity` | `./run.sh ../common/mesh/unit_square/D --comparison-rounds 1 1 1` | fast three-strategy classroom comparison |
| `06_airfoil_mesh_motion` | `python3 airfoil_step.py --reset`, then `python3 airfoil_step.py 0.50 U -0.02 --width 0.12 --smooth-iterations 200` | one safe fixed-topology boundary action |
| `07_poisson_shape_optimization` | `./run.sh output 4` | four coordinate-search sweeps with remeshing |
| `08_airfoil_digital_twin` | `./run_teaching_demo.sh` | reset and run the complete seed-2026 digital-twin experiment |
| `09_airfoil_barycentric_motion` | `./run.sh --steps 48 --smooth-iterations 0 --boundary-quality-iterations 20 --quality-floor 0.40` | fixed-topology disk-to-NACA0012 path |
| `10_airfoil_rl_dat` | `./run.sh --evaluate-only --max-steps 140 --seed 2026` | reload and evaluate the saved RL checkpoint |

Full RL retraining uses
`./run.sh --episodes 400 --max-steps 140 --seed 2026`. Select the Python
environment once in the first-time setup above.
