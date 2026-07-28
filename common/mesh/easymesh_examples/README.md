# Runnable EasyMesh geometries

These examples isolate mesh generation from the PDE solvers. Each `.d` file is
small enough to edit by hand and demonstrates one geometry concept.

## Run

Generate every example:

```bash
EASYMESH_BIN=/path/to/easymesh ./run_examples.sh
```

Use the default `$HOME/bin/easymesh` by omitting `EASYMESH_BIN`. Select another
output directory with:

```bash
./run_examples.sh /tmp/my-easymesh-output
```

EasyMesh outputs are never written beside the source files. The default output
directory is `output/`, which is ignored by Git.

The examples were verified with EasyMesh 1.4:

| example | nodes | triangles | sides |
|---|---:|---:|---:|
| `unit_square` | 139 | 236 | 374 |
| `l_shape` | 347 | 620 | 966 |
| `rectangle_with_hole` | 269 | 459 | 728 |

Counts can differ slightly with another EasyMesh build or changed spacing,
relaxation, or smoothing settings.

## Example 1: unit square

```text
3 (marker 4) 2
  +---------+
  |         |
  |         |
  +---------+
0 (marker 1) 1
```

`unit_square.d` uses uniform spacing `0.10` and assigns a distinct segment
marker to each side: bottom `1`, right `2`, top `3`, and left `4`.

## Example 2: nonconvex L-shape

```text
5 +---+ 4
  |   |
  |   +---+ 3
  |       |
0 +-------+ 1
          2
```

`l_shape.d` demonstrates a reentrant corner at point `3`. It is useful for
testing adaptive refinement because elliptic PDE solutions commonly have
reduced regularity at this corner. All outer sides use marker `1`.

## Example 3: rectangle with a square hole

```text
3 +-------------------+ 2
  |                   |
  |    5 +-----+ 6    |
  |      |     |      |
  |    4 +-----+ 7    |
  |                   |
0 +-------------------+ 1
```

`rectangle_with_hole.d` demonstrates multiple boundary loops. The outer loop
is counterclockwise with markers `1` and `2`; the hole is clockwise with
marker `3`. Its point spacing varies from `0.08` to `0.20`, illustrating local
mesh-size control.

## Use one result with AFEPack

After running the launcher, pass the generated basename:

```bash
cd ../../../00_poisson_basic
./run.sh ../common/mesh/easymesh_examples/output/unit_square
```

The Poisson example applies its own boundary-condition interpretation to the
markers. If a new PDE expects different marker numbers, edit the segment
markers in the `.d` file or adjust the PDE boundary-condition code.
