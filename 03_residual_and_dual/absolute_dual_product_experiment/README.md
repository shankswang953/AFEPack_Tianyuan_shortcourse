# Residual times absolute-dual experiment

This diagnostic run marks with the practical adjoint-weighted residual

```text
eta_mark,K = eta_residual,K * RMS_K(abs(z_H)).
```

It deliberately uses the enriched dual magnitude itself rather than the dual
correction `z_H - I_h^H Pi_h z_H`. The signed DWR correction is still computed
separately for diagnostics. Six rounds mark 20 elements each.

The heat source is centered at `(0.40, 0.80)` and the sensor at `(0.78, 0.22)`.
In round 1 all 20 marked elements lie near the source and none lie near the
sensor. In round 6, 9 of the 20 marked elements remain near the source and none
lie near the sensor. The product therefore does not create two refinement
centers: the very large primal residual at the source dominates the nonzero,
spatially broad elliptic dual.

Errors below use the level-5 uniform finite-element reference
`J_5 = 0.914110002727527`.

| level | primal DOFs | raw goal error | DWR-corrected error |
|---:|---:|---:|---:|
| 0 | 506 | `3.218e-6` | `8.277e-7` |
| 1 | 545 | `1.699e-6` | `7.526e-7` |
| 2 | 586 | `1.194e-6` | `6.624e-7` |
| 3 | 629 | `8.653e-6` | `1.923e-6` |
| 4 | 667 | `6.946e-6` | `2.110e-6` |
| 5 | 708 | `1.114e-5` | `2.194e-6` |
| 6 | 761 | `2.803e-5` | -- |

The first two refinements improve the raw quantity of interest, but continued
fixed-budget refinement overshoots. The mesh remains source-dominated and does
not visibly resolve the sensor footprint.

The factor used for marking can be reconstructed from the data columns as
`primal_residual_norm * dual_solution_rms`.
