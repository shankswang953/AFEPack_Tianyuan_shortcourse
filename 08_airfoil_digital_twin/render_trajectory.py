#!/usr/bin/env python3

"""Regenerate the saved policy-rollout animation from its CSV history.

Command-line usage:
    python3 render_trajectory.py [--fps N]

Reads `output/policy_rollout/data_history.csv` and the target airfoil data.
Writes `shape_evolution.gif` and `shape_final.png` in the same rollout
directory. Run a controller or policy rollout before using this helper.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from dt_airfoil.trajectory import render_trajectory


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fps", type=int, default=2)
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parent
    rollout = root / "output" / "policy_rollout"
    history = rollout / "data_history.csv"
    if not history.exists():
        raise SystemExit("run run_policy.py before rendering a trajectory")
    animation = rollout / "shape_evolution.gif"
    final_figure = rollout / "shape_final.png"
    render_trajectory(
        history,
        root / "data" / "target_naca0012.dat",
        animation,
        final_figure,
        frames_per_second=arguments.fps,
    )
    print(f"history:   {history}")
    print(f"animation: {animation}")
    print(f"final PNG: {final_figure}")


if __name__ == "__main__":
    main()

