# Discrete dual for a local cooling-boundary flux

This example uses the same mixed-boundary Poisson operator as the preceding
steady heat-conduction example. The bottom boundary is cooled (`T=0`) and the
other three sides are insulated. It uses the manufactured exact temperature

```text
T_ex(x,y) = y(2-y) [1 - 0.04 cos(2 pi x)]
```

and defines the heat source from the equation:

```text
Q(x,y) = -Delta T_ex
       = 2 - 0.04 [2 + 4 pi^2 y(2-y)] cos(2 pi x).
```

Thus `Q` is positive distributed heating that is strongest in the central
vertical column above the target patch. The numerical code is still a forward
solve: it receives `Q` and computes `T_h`; `T_ex` is retained only for checking
the result.

The target functional is the outward heat flux through that patch,

```text
J_h(U) = integral_patch -k grad(T_h).n ds.
```

For the bottom boundary `n=(0,-1)`, so the integrand is `k d_y T_h`. The code
uses the target patch `[0.40,0.60] x {0}`. Its exact value is

```text
J(T_ex) = 0.4 + 0.08/pi sin(pi/5).
```

The code computes the P1 primal solution and both residual views:

```text
strong RMS residual: sqrt(|K|^-1 integral_K |Q + Delta T_h|^2),
full indicator      : h_K^2 ||Q + Delta T_h||_K^2 plus flux jumps.
```

It then
assembles `g_h = grad_U J_h` from the P1 basis gradients and solves the fully
discrete dual

```text
A_h^T psi_h = -g_h.
```

The primal residual and dual are computed on the same uniformly level-2 refined
mesh so that their spatial distributions can be compared without first making
any adaptive refinement.

Important: for this linear Poisson problem, the dual depends on the Jacobian
`A_h` and the target functional `J_h`, not on the heat source. The residual
shows where the primal equation is poorly resolved; the dual shows where those
errors can influence the target flux.

Build and run:

```bash
make
chmod +x run.sh
./run.sh
```

Main outputs:

- `output/manufactured_strong_residual.svg`: strong residual RMS;
- `output/manufactured_full_residual_indicator.svg`: cell residual plus flux
  jumps;
- `output/manufactured_residual.dat`: all residual components;
- `output/manufactured_temperature.dx`: numerical primal temperature;
- `output/boundary_flux_dual_magnitude.svg`: elementwise RMS dual magnitude
  with a linear color scale, including its global elliptic tail;
- `output/boundary_flux_dual_magnitude_focused.svg`: the same values with a
  squared color transform that makes the peak near the target patch clearer;
- `output/boundary_flux_dual_magnitude.dat`: elementwise values;
- `output/boundary_flux_dual_signed.dx`: signed dual for OpenDX;
- `output/boundary_flux_dual_signed.dat`: nodal signed values.

<!-- script-interface -->
## Launcher interface and portability

`./run.sh [ROOT_MESH]` builds the solver and writes all primal/dual fields,
residual tables, and figures to `output/`. The mesh defaults to
`../common/mesh/unit_square/D`. Set the shared AFEPack/OpenBLAS/Boost build
variables and AFEPack runtime/template paths from the parent README on systems
that do not use the supplied `$HOME` and `/opt/local` defaults.
