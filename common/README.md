# Shared resources and EasyMesh tutorial

This directory contains checked-in meshes and geometry examples shared by the
numbered AFEPack exercises. Generated solver results should stay in each
exercise's own `output/` directory.

## What EasyMesh does

EasyMesh turns a two-dimensional boundary description into a triangular mesh:

```text
geometry.d  --EasyMesh-->  geometry.n + geometry.e + geometry.s
```

The three output files contain:

| file | contents |
|---|---|
| `NAME.n` | node index, coordinates, and boundary marker |
| `NAME.e` | triangle vertices, neighbors, sides, and circumcenter data |
| `NAME.s` | edge endpoints, adjacent triangles, and boundary marker |

AFEPack normally receives the common basename without a suffix. For example,
the files `mesh/unit_square/D.n`, `D.e`, and `D.s` are passed to an AFEPack
program as `mesh/unit_square/D`.

## 1. Select the EasyMesh executable

The launchers in this repository look for EasyMesh at
`$HOME/bin/easymesh`. If it is installed elsewhere, define:

```bash
export EASYMESH_BIN=/path/to/easymesh
```

Check the installation and display the EasyMesh 1.4 options with:

```bash
"$EASYMESH_BIN"
```

The most useful noninteractive options are:

| option | meaning |
|---|---|
| `-g` | disable the graphics window |
| `-m` | suppress progress messages |
| `-r` | disable relaxation |
| `-s` | disable Laplacian smoothing |
| `-d` | read the domain but do not triangulate it |
| `+dxf` | additionally write a DXF drawing |
| `+fig` | additionally write an Xfig drawing |
| `+example` | create EasyMesh's built-in `example.d` |

## 2. Understand a `.d` file

Ignoring comments and blank lines, an EasyMesh input file has two sections:

```text
NUMBER_OF_BOUNDARY_POINTS
point_id: x y target_spacing point_marker
...
NUMBER_OF_BOUNDARY_SEGMENTS
segment_id: start_point end_point segment_marker
...
```

For example, a unit square with a different marker on each side is:

```text
4
0: 0.0 0.0 0.10 1
1: 1.0 0.0 0.10 2
2: 1.0 1.0 0.10 3
3: 0.0 1.0 0.10 4
4
0: 0 1 1
1: 1 2 2
2: 2 3 3
3: 3 0 4
```

The columns mean:

- `x y`: boundary-point coordinates;
- `target_spacing`: the requested local mesh spacing near that point; smaller
  values create a finer mesh;
- `point_marker`: an integer tag attached to that boundary point;
- `start_point end_point`: indices of the endpoints of a boundary segment;
- `segment_marker`: the boundary-condition tag carried by the generated
  boundary sides.

Markers have no universal physical meaning. The PDE code decides whether a
marker means Dirichlet, Neumann, wall, inlet, and so on. In the square example
above, markers `1`, `2`, `3`, and `4` mean bottom, right, top, and left.

List an outer boundary counterclockwise. List every hole as a separate closed
loop in clockwise order. Segments must not cross, point and segment indices
must be unique, and every loop must close.

## 3. Generate a mesh manually

Work in a scratch directory so generated files do not overwrite checked-in
data:

```bash
mkdir -p /tmp/easymesh-unit-square
cp mesh/easymesh_examples/unit_square.d /tmp/easymesh-unit-square/
cd /tmp/easymesh-unit-square
"$EASYMESH_BIN" unit_square.d -g -m
```

A successful run creates:

```text
unit_square.n
unit_square.e
unit_square.s
```

EasyMesh 1.4 can return exit status `1` even after successfully writing a
valid mesh. Therefore the repository launchers verify that all three output
files are nonempty instead of relying only on the process status.

Inspect the output counts with:

```bash
head -n 1 unit_square.n
head -n 1 unit_square.e
head -n 1 unit_square.s
```

## 4. Run all teaching geometries

The `mesh/easymesh_examples/` directory contains:

1. `unit_square.d`: four boundary markers and uniform spacing;
2. `l_shape.d`: a nonconvex domain with a reentrant corner;
3. `rectangle_with_hole.d`: an outer loop plus a clockwise inner loop, with
   locally varying spacing.

Generate all three without modifying the source files:

```bash
cd mesh/easymesh_examples
EASYMESH_BIN=/path/to/easymesh ./run_examples.sh
```

The generated `.d`, `.n`, `.e`, and `.s` files are placed in
`mesh/easymesh_examples/output/`.

## 5. Use a generated mesh with AFEPack

Pass the basename, not an individual suffix:

```bash
cd ../../../00_poisson_basic
./run.sh ../common/mesh/easymesh_examples/output/unit_square
```

For code that expects AFEPack's `.mesh` format, configure and run the optional
converter:

```bash
export EASYMESH2MESH_BIN=/path/to/easymesh2mesh
"$EASYMESH2MESH_BIN" unit_square unit_square.mesh
```

## Common problems

- **`EasyMesh executable not found`**: set `EASYMESH_BIN` to an executable
  absolute path.
- **No mesh files appear**: run without `-m`, then check point indices,
  segment indices, loop closure, and intersecting segments.
- **The hole is filled**: reverse the hole loop so that it is clockwise.
- **The mesh is too coarse or too large**: decrease the point spacing; larger
  spacing produces fewer triangles.
- **Boundary conditions are applied to the wrong side**: inspect the marker
  column in `NAME.s` and make the segment markers match the PDE code.

See [`mesh/README.md`](mesh/README.md) for the mesh catalog and
[`mesh/easymesh_examples/README.md`](mesh/easymesh_examples/README.md) for
the geometry of each runnable example.

The mesh catalog also documents
[`mesh/visualize_dx.py`](mesh/visualize_dx.py), a Matplotlib-based viewer for
AFEPack/OpenDX mesh and scalar-field files.
