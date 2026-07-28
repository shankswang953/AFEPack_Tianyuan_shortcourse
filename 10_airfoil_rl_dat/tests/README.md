# Pure-`dat` RL tests

Run from the parent example directory:

```bash
python3 -m unittest discover -s tests -v
```

The tests read the checked-in target `dat` file, use in-memory environments,
and create no persistent output.
