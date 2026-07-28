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
