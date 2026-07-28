# Forced continuation of the DWR refinement

This auxiliary run repeats the residual-times-dual-correction marking step six
times, always marking 20 elements. It is intentionally separate from the
one-step default example.

| DWR level | DOFs | raw absolute goal error | signed DWR estimate |
|---:|---:|---:|---:|
| 0 | 506 | `3.221e-6` | `+2.390e-6` |
| 1 | 551 | `5.879e-7` | `-3.434e-7` |
| 2 | 601 | `4.290e-6` | `-4.148e-6` |
| 3 | 653 | `2.133e-5` | `-1.766e-5` |
| 4 | 715 | `4.232e-5` | `-3.255e-5` |
| 5 | 784 | `7.032e-5` | `-5.606e-5` |
| 6 | 845 | `1.168e-4` | not recomputed |

The first DWR step improves the raw goal error by a factor of about 5.5. The
signed estimate then changes from positive to negative, indicating that the
adaptive step has crossed the zero of the signed goal error. Repeating the same
fixed-budget step nevertheless refines the same error channel and increasingly
overshoots the goal value.

This is not evidence that DWR is invalid. It shows that a signed goal error is
not a norm and that fixed-count marking needs a stopping or step-control rule.
For this example, the estimator sign change is a useful warning to stop or to
reduce the marking budget substantially.

Files:

- `dwr_history.dat`: raw and corrected functional histories;
- `functional_comparison.dat`: residual, dual-magnitude, and forced DWR runs;
- `functional_error_vs_dofs.svg`: log-scale error plot;
- `dwr_indicator_round_2.svg`: indicator immediately after the sign change;
- `dwr_mesh_after.svg`: mesh after six forced steps.
