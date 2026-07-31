# ML-iteration tests

Run all tests from the parent example directory:

```bash
python3 -m unittest discover -s tests -v
```

The geometry and point-cloud tests use checked-in/in-memory data. The model
test uses an automatically removed temporary directory. No persistent output
is written to the project.
