# Barycentric-motion C++ backend

The backend builds `generate_airfoil_geometry` and `move_and_smooth`, used
internally by `../barycentric_motion.py`.

```bash
make
```

Define `CXX`, `AFEPACK_PREFIX`, `OPENBLAS_PREFIX`, and `BOOST_INCLUDE` in the
top-level optional `course_config.local` when the defaults do not match the
current system. Object files and executables are build products; continuation
results are written to `../output/`.
