# Two-heater goal-oriented adaptivity

This example introduces residual refinement, a fully discrete dual, and DWR
refinement through a steady heat-conduction problem.  The construction
deliberately separates the largest equation defect from the region that
controls the requested measurement.

## Physical problem

Let `Omega = (0,1)^2`.  The temperature satisfies

```text
-Delta T = f_A + f_B  in Omega,
       T = 0          on the whole boundary.
```

The boundary is held at zero temperature.  The two volumetric heaters are

```text
f_A(x) = 240 exp(-|x-(0.15,0.85)|^2 / (2 * 0.050^2)),
f_B(x) =  25 exp(-|x-(0.55,0.40)|^2 / (2 * 0.080^2)).
```

Heater A is the strong, localized distractor.  It dominates the ordinary
cell residual but is far from the measurement region and close to a cooled
boundary.  Heater B is weaker.  Heat from B must propagate toward a sensor in
the opposite, lower-right part of the plate.

## Quantity of interest and reference value

The output is the normalized average temperature reported by a finite-size
sensor centered at `x_s = (0.82,0.18)`:

```text
J(T) = integral_Omega w(x) T(x) dx,

w(x) = exp(-|x-x_s|^2 / (2 * 0.120^2)) / Z,
Z    = integral_Omega exp(-|x-x_s|^2 / (2 * 0.120^2)) dx.
```

This is a local average rather than a singular point evaluation.  It also has
a reproducible reference value.  With homogeneous Dirichlet conditions, the
Laplacian eigenfunctions are

```text
phi_mn(x,y) = 2 sin(m pi x) sin(n pi y),
lambda_mn   = pi^2 (m^2+n^2),
```

so the code evaluates

```text
J_ref = sum_(m,n>=1) <f,phi_mn> <w,phi_mn> / lambda_mn.
```

Using 160 modes per direction gives

```text
J_ref = 0.0742880436830910.
```

Changing the truncation from 120 to 160 modes changes the result by only
`2.64e-14`.

## Residual refinement

The residual-only sequence marks with the P1 cell term

```text
eta_R,K^2 = h_K^2 ||f + Delta T_h||^2_(L2(K)).
```

For affine P1 elements and constant conductivity, `Delta T_h = 0` inside each
triangle.  The indicator therefore follows the two physical heaters, with the
strong heater A taking most of the budget.  Interior flux jumps are computed
and exported as diagnostics but are not used for marking in this introductory
comparison.

Each round marks the largest 5% of the active-element indicators.  The first
round has 47 marks: 38 lie near heater A and 9 near heater B.  By round ten,
159 of 186 marks are near A and 27 near B.  Residual refinement improves the
global finite-element solution, but it does not know which part of that
solution controls `J`.

## Fully discrete dual and localized DWR indicator

Write the assembled primal residual and the discrete functional gradient as

```text
R_h(U) = A_h U - b_h,
g_h    = grad_U J_h(U).
```

The algebraic dual is defined by

```text
(dR_h/dU)^T psi_h = -g_h.
```

The sign follows from `L(U,psi) = J(U) + psi^T R(U)`.  No continuous dual PDE
is needed by this implementation.

For DWR marking, the primal remains P1 on the current adaptive mesh.  The
dual is solved with P1 elements on one uniformly h-refined copy, and its
coarse projection is removed.  The code localizes the residual action into
volume and face terms and marks with

```text
eta_DWR,K = | rho_K(T_h)(psi_h^+ - I_h psi_h^+) |.
```

This is different from multiplying a cell-residual norm by `|psi_h|`.  The
projection difference supplies the unresolved dual scale, while the localized
residual action retains both magnitude and spatial coupling.

In the first DWR round, only 8 of 47 marks lie near distractor A; 9 lie near
heater B, 18 near the sensor, and the remaining marks begin to resolve the
region coupling B to the sensor.  The final DWR mesh continues to refine that
influence path rather than only the strongest source.

## Fair numerical comparison

Both strategies start from the same 506-DOF mesh.  Because mesh closure makes
one DWR round add more cells than one residual round, the recommended run uses
different round counts to compare approximately the same final primal budget:

```bash
./run.sh ../common/mesh/unit_square/D --comparison-rounds 10 10 7
```

| strategy | rounds | final DOFs | `J(T_h)` | `|J(T_h)-J_ref|` |
|---|---:|---:|---:|---:|
| residual refinement | 10 | 2322 | 0.0742546148 | `3.3429e-5` |
| dual-magnitude refinement | 10 | 2155 | 0.0742254651 | `6.2579e-5` |
| localized DWR refinement | 7 | 2107 | 0.0742900990 | `2.0553e-6` |

The DWR error is about 16.3 times smaller even though its final primal mesh
has fewer degrees of freedom.  At the nearest lower residual budget, 1987
DOFs, the residual error is `2.9957e-5`, still about 14.6 times larger.

The raw target error is not guaranteed to decrease after every refinement.
Galerkin energy error has a monotonicity structure for this symmetric Poisson
problem, but a general output `J(T_h)` does not.  The meaningful claim is that
DWR spends a comparable budget more effectively for the requested output,
not that every plotted point must be monotone.

The independent dual-magnitude sequence is retained as a diagnostic control.
It resolves the sensor sensitivity but then plateaus because the dual alone
contains no primal-defect information.

After ten dual-magnitude rounds it has 2155 DOFs, comparable with the other
two final meshes.  Its target error remains `6.2579e-5`.  The indicator
distribution explains the plateau: all 172 elements marked in round ten have
centroids within distance `0.15` of the sensor, and their mean distance from
the sensor is only `0.058`.  The sequence keeps resolving the already smooth
sensor neighborhood while ignoring whether the primal equation has a defect
there.

## Build and run

```bash
make
./run.sh ../common/mesh/unit_square/D --comparison-rounds 10 10 7
```

The default root mesh is `../common/mesh/unit_square/D`.  Equal round counts
remain available with `--rounds N`.  To solve only on a uniformly refined mesh,
use, for example,

```bash
./run.sh ../common/mesh/unit_square/D --uniform-reference 5
```

The recommended comparison is written to
`output/comparison_residual_10_dual_10_dwr_7/`.  Important files are:

- `TEACHING.md`: a short four-part classroom narrative and a map of the
  text-free slide assets;
- `teaching/`: clean SVG and PNG figures with no embedded labels.  It includes
  the problem fields, first/final indicator triptychs, and final-mesh
  triptych;

- `source_profile.svg` and `sensor_weight.svg`: the PDE source and functional;
- `residual_indicator_round_1.svg` and `dwr_indicator_round_1.svg`: the first
  marking decisions;
- `mesh_after.svg` and `dwr_mesh_after.svg`: the two final meshes;
- `dual_magnitude.svg`: the discrete dual magnitude;
- `dual_magnitude_round_10.svg`: the final dual-magnitude indicator and marks;
- `functional_error_vs_dofs.svg`: the target-error comparison;
- `functional_comparison.dat`: all values used in the comparison plot;
- `dwr_indicator.dat`: signed and absolute localized DWR contributions plus
  diagnostic residual and dual factors;
- `dwr_history.dat`: raw and DWR-corrected functional values.

<!-- script-interface -->
## Launcher interface and output selection

```text
./run.sh [ROOT_MESH]
./run.sh [ROOT_MESH] --rounds N
./run.sh [ROOT_MESH] --comparison-rounds RESIDUAL_ROUNDS DWR_ROUNDS
./run.sh [ROOT_MESH] --comparison-rounds RESIDUAL_ROUNDS DUAL_ROUNDS DWR_ROUNDS
```

The default output is `output/`. Non-default round counts use
`output/rounds_N/` or a `output/comparison_.../` directory. The root mesh
defaults to `../common/mesh/unit_square/D`. Define the shared build and
AFEPack runtime variables from the parent README when needed.
