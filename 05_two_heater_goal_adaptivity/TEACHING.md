# Teaching sequence

This version keeps the 506-DOF initial mesh and the verified
`10 / 10 / 7` comparison.  A compact explanation works better than starting
with the DWR formula.

## 1. Begin with the question

The plate is held at zero temperature on its boundary and satisfies

```text
-Delta T = f  in Omega,       T = 0 on boundary(Omega).
```

There are two heaters.  The strong heater in the upper-left corner dominates
the equation defect.  The sensor reports a weighted average temperature near
the lower-right corner:

```text
J(T) = integral_Omega w T dx.
```

Use `problem_fields_pair.png`.  The left panel is `f`; the right panel is
`w`.  This makes the key question concrete: should the mesh resolve the
largest source, the sensor, or the way an unresolved source affects the
sensor?

## 2. Introduce the three indicators one at a time

Use `indicator_initial_triptych.png`, or place its three component images on
successive overlays.  From left to right they are:

1. cell residual: where the discrete primal equation is locally inaccurate;
2. dual magnitude: where the requested measurement is sensitive;
3. localized DWR: where unresolved primal error matters to that measurement.

The fill shows relative indicator magnitude.  A dark outline marks an element
in the largest 5% selected for refinement.  The colors have a fixed role:
orange is residual, teal is dual magnitude, and violet is DWR.

The main teaching point is visible before showing any algebra: residual and
dual answer different questions, while DWR uses both pieces of information.

## 3. Then show the discrete formulas

For affine P1 elements, the introductory residual marking term is

```text
eta_R,K^2 = h_K^2 ||f + Delta T_h||^2_L2(K).
```

The discrete dual is the transpose-Jacobian solve

```text
(dR_h/dU)^T psi_h = -grad J_h(U).
```

The localized DWR indicator is

```text
eta_DWR,K = |rho_K(T_h)(psi_h^+ - I_h psi_h^+)|.
```

Emphasize that the DWR indicator is not simply residual magnitude times dual
magnitude.  The enriched dual correction supplies the sensitivity that the
current primal space cannot represent.

## 4. Compare where the DOFs went

Use `final_meshes_triptych.png`.  From left to right: residual refinement,
dual-magnitude refinement, and localized DWR refinement.

- Residual refinement spends most of its budget around the two heaters,
  especially the strong distractor.
- Dual-magnitude refinement concentrates almost entirely at the sensor and
  eventually plateaus.
- DWR refinement resolves the source-to-sensor influence pattern.

The final primal budgets are comparable:

| strategy | rounds | final DOFs | target error |
|---|---:|---:|---:|
| residual | 10 | 2322 | `3.3429e-5` |
| dual magnitude | 10 | 2155 | `6.2579e-5` |
| localized DWR | 7 | 2107 | `2.0553e-6` |

Finish with `figures/functional_error_vs_dofs.png`.  The safe conclusion is not that
every target-error value must decrease monotonically.  It is that, at a
comparable final budget, DWR allocates the mesh much more effectively for the
chosen measurement.

## Figure assets

The two final teaching figures distributed with the examples are PNG files:

```text
../assets/results/05_final_meshes_triptych.png
../assets/results/05_functional_error_vs_dofs.png
```

Fresh runs create problem-field, indicator, final-mesh, and convergence PNGs
in `figures/`.  All detailed field and indicator values remain available
under `fields/` for another plotting tool if desired.
