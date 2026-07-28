# Sign-balanced DWR marking experiment

This diagnostic run retains the sign of every localized DWR contribution and
splits the elements into

```text
P = {K : eta_K > 0},   N = {K : eta_K < 0}.
```

Minimal Dörfler sets with `theta = 0.5` are selected independently in `P` and
`N`, and their union is refined. This prevents the severe one-sided loss of
cancellation observed when all elements are ranked by `abs(eta_K)` together.

The experiment uses the level-5 uniform finite-element value
`J_5 = 0.914110002727527` for the errors below.

| level | primal DOFs | marked in the preceding step | raw goal error | DWR-corrected error |
|---:|---:|---:|---:|---:|
| 0 | 506 | -- | `3.218e-6` | `8.277e-7` |
| 1 | 989 | 235 | `5.000e-6` | `1.258e-6` |
| 2 | 1,537 | 252 | `7.820e-6` | `1.751e-7` |
| 3 | 2,816 | 709 | `6.213e-6` | `1.085e-6` |
| 4 | 4,667 | 998 | `4.399e-6` | `9.307e-7` |
| 5 | 7,153 | 1,432 | `2.535e-6` | `5.638e-7` |
| 6 | 13,539 | 3,786 | `1.762e-6` | -- |

The sign-balanced rule removes the catastrophic growth to `1e-4`, and the
absolute DWR sum decreases from `6.52e-4` to `4.06e-5` by the last evaluated
indicator. It does not, however, give an efficient monotone raw-functional
sequence. For comparison, two global refinements use 7,601 DOFs and give a
level-5-referenced goal error of about `2.03e-7`.

The diagnostic therefore supports two conclusions: signed balance matters in
this cancellation-dominated example, while the raw goal value remains a poor
teaching convergence curve. The DWR-corrected functional is substantially more
informative than the raw value.
