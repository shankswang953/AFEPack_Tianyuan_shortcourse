# Residual, discrete dual, and DWR refinement

This reproducible AFEPack example solves the steady heat-conduction problem

```text
-div(grad T) = Q(x)  in the unit square,
T = 0               on the bottom boundary,
grad T . n = 0      on the other three boundaries.
```

`Q` is a Gaussian heat source centered at `(0.40, 0.80)` with width `0.075`.
Boundary mark `1` is the cooled bottom edge. The other boundaries are
insulated. The source and sensor are deliberately separated so that residual,
dual, and DWR refinement have visibly different meanings.

## Residual refinement

The residual-only sequence marks with

```text
eta_cell,K^2 = h_K^2 ||Q + div(grad T_h)||_K^2.
```

For P1 elements with constant conductivity, `div(grad T_h) = 0` inside each
triangle. The strong cell residual therefore follows the localized heat
source. Four rounds mark `20`, `40`, `80`, and `160` elements.

The code also exports the full P1 residual estimator, including interior
normal-flux jumps and natural-boundary fluxes. `strong_residual.svg` omits the
mesh-size factor and is the cleanest figure for explaining the physical
residual distribution.

## Quantity of interest

The functional is the temperature reported by a finite-size sensor centered at
`x_s = (0.78, 0.22)`:

```text
J(T) = integral_Omega w_s(x) T(x) dx,

w_s(x) = exp(-|x-x_s|^2 / (2 sigma_s^2)) / Z,
sigma_s = 0.12,   integral_Omega w_s(x) dx = 1.
```

This is a local average rather than a singular point value. The larger sensor
footprint is intentional: an extremely narrow sensor makes the enriched dual
correction almost singular and causes a nominal DWR indicator to be visually
indistinguishable from dual-magnitude refinement.

The numerical reference is evaluated in the exact mixed-boundary Laplacian
eigenbasis using 80 modes per direction and composite Simpson quadrature:

```text
J_ref = 0.914110005950467.
```

Repeating with 60 modes changes the result by about `2.12e-12`.

An independent finite-element reference can be computed by constructing the
fifth globally refined mesh directly and solving only on that mesh:

```text
./run.sh ../common/mesh/unit_square/D --uniform-reference 5
```

This mode does not solve levels 0--4 and does not run the adaptive, dual, or
DWR sequences. It also avoids writing the nearly one-million-element
temperature field. The resulting mesh has 952,320 triangles and 477,441 P1
DOFs. On the current machine the complete run takes about 40 seconds and gives

```text
J_5 = 0.9141100027274,
|J_5 - J_spectral| = 3.22e-9.
```

The summary is written to `output/uniform_reference_level_5.dat`.

## Fully discrete dual

With the assembled residual `R_h(U) = A_h U - b_h` and
`g_h = grad_U J_h(U)`, the code solves

```text
(dR_h/dU)^T psi_h = -g_h,
or, here, A_h^T psi_h = -g_h.
```

No continuous dual PDE is needed. The independent dual-magnitude sequence is
retained only to visualize sensitivity near the sensor; dual magnitude alone
is not a DWR indicator.

## DWR: residual weighted by dual error

The DWR calculation uses a P1 primal solution on `T_h` and a fully discrete P1
dual on one uniformly h-refined copy `T_H`. The fine dual is projected to the
coarse nodes and prolonged back:

```text
delta Z = Z_H - I_h^H Pi_h Z_H.
```

The signed functional-error estimate is

```text
eta = R_H(I_h^H U_h)(delta Z) = sum_K eta_K.
```

For marking, the code exposes the two factors rather than balancing positive
and negative local signs:

```text
eta_mark,K = eta_residual,K * omega_dual,K,

eta_residual,K = full primal residual norm,
omega_dual,K   = RMS_K(delta Z).
```

This is the local Cauchy--Schwarz form of the DWR action. It is a product at the
same element, not a union of a residual mesh and a dual mesh. An element is
important only when its primal defect and its influence on the sensor are both
significant.

The example performs one DWR marking step with a budget of 20 elements. Of
these, 12 lie near the source and 8 near the sensor footprint, so the result is
not merely dual-magnitude refinement. Small adaptive steps are used because
signed goal errors need not decrease monotonically under repeated fixed-budget
refinement.

A separate forced six-step run is stored in `continued_experiment/`. It shows
that the signed DWR estimate changes sign after the first successful step; if
the same 20-element budget is nevertheless repeated, the raw goal error grows.
This is why the default teaching run stops after one step.

Measured values are:

| mesh/strategy | DOFs | absolute goal error |
|---|---:|---:|
| root mesh | 506 | `3.221e-6` |
| residual, one step | 547 | `1.231e-6` |
| dual magnitude, one step | 542 | `5.356e-6` |
| DWR product, one step | 551 | `5.879e-7` |

At essentially the same DOF count, the DWR mesh has about half the goal error
of residual-only refinement and about one fifth of the root-mesh error. On the
root mesh, the signed DWR correction also reduces the error from `3.221e-6` to
`8.309e-7`.

## Build and run

```bash
make
chmod +x run.sh
./run.sh
```

The default root mesh is `../common/mesh/unit_square/D`. Alternative source and
sensor parameters can be supplied as

```text
./run.sh ROOT_MESH SOURCE_X SOURCE_Y SOURCE_SIGMA \
         SENSOR_X SENSOR_Y SENSOR_SIGMA
```

Important outputs:

- `strong_residual.svg`: physical strong residual around the heat source;
- `dual_magnitude.svg`: independent discrete dual magnitude;
- `dwr_primal_residual_factor.svg`: first DWR marking factor;
- `dwr_dual_correction_factor.svg`: second DWR marking factor;
- `dwr_indicator.svg`: their elementwise product and the selected elements;
- `dwr_indicator.dat`: signed DWR action, both marking factors, their product,
  and the mark flag;
- `dwr_mesh_after.svg`: mesh after the DWR step;
- `functional_comparison.dat`: residual, dual-magnitude, and DWR functional
  errors versus DOFs;
- `dwr_history.dat`: raw and DWR-corrected functional values;
- `functional_error_vs_dofs.svg`: comparison plot.

<!-- script-interface -->
## Launcher interface and portability

Use `./run.sh [ROOT_MESH] [PROGRAM_ARGUMENTS...]`. The root mesh defaults to
`../common/mesh/unit_square/D`; all remaining arguments are passed to the C++
program, including `--uniform-reference LEVEL` and the documented
source/sensor positional values. All generated tables, fields, meshes, and
figures go to `output/`.

Other users should set the AFEPack/OpenBLAS/Boost build variables and the
`AFEPACK_PATH`/`AFEPACK_TEMPLATE_PATH` runtime variables listed in the
parent README when the defaults do not match their installation.
