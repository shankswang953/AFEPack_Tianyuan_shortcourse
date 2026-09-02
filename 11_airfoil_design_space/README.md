# Low-dimensional airfoil design-space sampling

This example separates **design freedom** from **boundary resolution**. Each
airfoil is controlled by only five engineering descriptors, although its UIUC
`data` file contains 66 points on each surface and the existing Bezier pipeline
generates 96 EasyMesh boundary points.

## Five design parameters

| parameter | meaning |
|---|---|
| `m` | signed maximum camber divided by chord |
| `x_c` | chordwise location of maximum camber |
| `t` | maximum thickness divided by chord |
| `x_t` | chordwise location of maximum thickness |
| `t_te` | trailing-edge thickness divided by chord |

The camber line follows the familiar NACA four-digit construction. A monotone
coordinate warp moves the thickness maximum, and a localized term controls the
trailing-edge gap. This is therefore a **NACA-inspired teaching family**, not
an exact named NACA series.

## Sampling and rejection

The launcher first generates a deterministic Latin-hypercube design in a box
wider than the accepted engineering envelope. Every candidate is then checked
before any mesh or CFD work:

- both surfaces must use increasing `x` coordinates;
- the stored upper/lower order must give the expected boundary orientation;
- the upper surface must stay above the lower surface;
- the closed contour must not intersect itself;
- area, curvature, thickness, camber, and trailing-edge gap must stay in the
  prescribed engineering envelope.

Negative camber is not automatically invalid: it is a legitimate geometry.
"Upside down" here means reversed surface storage or the wrong closed-contour
orientation, not merely a negative value of `m`.

From the valid candidates, greedy maximin selection chooses the next point
farthest from the selected set in normalized five-dimensional parameter space.
Five additional ranked reserves are kept by default. Only these dispersed
designs are passed through the existing
`data -> Bezier fit -> boundary points -> EasyMesh` pipeline.
An EasyMesh result with invalid vertex indices, degenerate triangles, or
minimum quality below 0.40 is rejected and replaced by the next reserve.

## Run

```bash
./run.sh
```

The deterministic defaults are:

```text
2048 wide Latin-hypercube proposals
6 dispersed, meshable airfoils (+5 ranked reserves)
seed 2026
```

They can be changed explicitly:

```bash
./run.sh --candidates 4096 --selected 9 --reserve 5 --seed 2026
```

The main outputs are:

```text
output/candidates.csv          every proposal and rejection reason
output/accepted.csv            geometrically admissible proposals
output/rejected.csv            rejected proposals
output/selected.csv            dispersed designs intended for initial CFD
output/ranked.csv              selections followed by replacement candidates
output/data/*.dat              high-resolution boundary representations
output/cases/*/airfoil.mesh    corresponding EasyMesh/AFEPack meshes
output/mesh_quality.csv        node/cell counts and minimum triangle quality
output/mesh_rejections.csv     post-EasyMesh failures replaced by reserves
output/figures/design_space_selection.png
output/figures/validator_examples.png
output/figures/selected_airfoils.png
output/figures/selected_meshes.png
```

No CFD equation is solved. The selected set is the initial design of
experiments; CFD outputs can later train a surrogate and drive active learning.
