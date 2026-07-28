# Shared meshes

This directory groups checked-in root meshes used by the AFEPack examples.
`unit_square/` provides the default mesh for examples 00--05.
`easymesh_examples/` contains three small `.d` geometries and a launcher that
generates their meshes in an ignored output directory.

Pass a mesh basename (for example `unit_square/D`) to a launcher. AFEPack and
EasyMesh add the `.n`, `.e`, and `.s` suffixes themselves.

| directory | purpose |
|---|---|
| `unit_square/` | checked-in `D.[nse]` default mesh plus its `D.d` source |
| `easymesh_examples/` | unit square, nonconvex L-shape, and rectangle-with-hole tutorials |

Start with the complete EasyMesh walkthrough in
[`../README.md`](../README.md).

## Visualize an OpenDX `.dx` file

`visualize_dx.py` creates an ordinary PNG, SVG, or PDF without requiring the
legacy OpenDX graphical application. It supports both AFEPack mesh-only files
and files containing a nodal or element scalar field.

From the repository root, visualize a mesh-only file:

```bash
python3 common/mesh/visualize_dx.py \
  01_local_refinement/output/D_root.dx \
  -o /tmp/D_root.png
```

Visualize the finite-element solution and its color scale:

```bash
python3 common/mesh/visualize_dx.py \
  00_poisson_basic/output/u.dx \
  -o /tmp/u.png
```

To draw only the triangulation from a solution file and show all nodes:

```bash
python3 common/mesh/visualize_dx.py \
  00_poisson_basic/output/u.dx \
  --mesh-only --show-nodes \
  -o /tmp/u_mesh.svg
```

If `-o/--output` is omitted, the figure is written next to the input file with
the same basename and a `.png` suffix. The output directory is created
automatically. Use `--no-edges` for a smoother scalar-field image, `--cmap`
to select a Matplotlib color map, `--title` to replace the title, and `--dpi`
to set raster resolution:

```bash
python3 common/mesh/visualize_dx.py --help
```

The only plotting dependency is Matplotlib. If it is installed in a
non-default Python environment, invoke that interpreter explicitly:

```bash
/path/to/python common/mesh/visualize_dx.py INPUT.dx -o OUTPUT.png
```
