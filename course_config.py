"""Load the optional repository-wide machine configuration for Python tools."""

from __future__ import annotations

import os
from pathlib import Path


def load_course_config(repository_root: Path) -> Path | None:
    """Load ``course_config.local`` into the process environment if it exists.

    The file deliberately uses simple ``NAME=value`` lines so the same local
    configuration can be read by POSIX shell, Make, and Python.
    """

    config_file = repository_root.resolve() / "course_config.local"
    if not config_file.is_file():
        return None

    for line_number, raw_line in enumerate(
        config_file.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(
                f"{config_file}:{line_number}: expected NAME=value"
            )
        name, value = line.split("=", 1)
        name = name.strip()
        value = value.strip()
        if (
            not name
            or not name.replace("_", "A").isalnum()
            or name[0].isdigit()
        ):
            raise ValueError(
                f"{config_file}:{line_number}: invalid variable name {name!r}"
            )
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
            value = value[1:-1]
        os.environ[name] = os.path.expandvars(os.path.expanduser(value))

    return config_file
