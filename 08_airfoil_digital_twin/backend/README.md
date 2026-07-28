# Digital-twin C++ backend

This directory builds two internal executables:

- `generate_airfoil_geometry`: fit `dat` data and write EasyMesh geometry;
- `move_and_smooth`: move/validate/smooth an AFEPack mesh.

They are normally built automatically by `../prepare_case.py`. To build
manually:

```bash
make
```

Define `CXX`, `AFEPACK_PREFIX`, `OPENBLAS_PREFIX`, and `BOOST_INCLUDE`
when the defaults do not match the current system. `OPENBLAS_PREFIX` and
`BOOST_INCLUDE` default to MacPorts paths under `/opt/local`.
Object files and executables are build products; experiment results go to
`../output/`.
