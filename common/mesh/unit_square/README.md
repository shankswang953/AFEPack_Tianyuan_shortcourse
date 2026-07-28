# Unit-square root mesh

`D.d` is the EasyMesh boundary description. `D.n`, `D.e`, and `D.s` are the
checked-in EasyMesh node, element, and side files used by examples 00--05.
The four side markers are:

| marker | side |
|---:|---|
| `1` | bottom |
| `2` | right |
| `3` | top |
| `4` | left |

Use the basename without a suffix:

```bash
./run.sh ../common/mesh/unit_square/D
```

To regenerate the mesh safely, copy `D.d` to a scratch directory:

```bash
mkdir -p /tmp/afepack-unit-square
cp D.d /tmp/afepack-unit-square/
cd /tmp/afepack-unit-square
"${EASYMESH_BIN:-$HOME/bin/easymesh}" D.d -g -m
```

With EasyMesh 1.4, this `D.d` reproduces the checked-in mesh with 506 nodes,
930 triangles, and 1,435 sides.

These checked-in files are inputs and should not be overwritten by example
outputs. See [`../../README.md`](../../README.md) for the full EasyMesh
tutorial.
