#!/usr/bin/env python3

"""Train the reward-prediction digital twin from recorded real transitions.

Command-line usage:
    python3 train_twin.py [--epochs N] [--seed N]

Reads `output/replay.jsonl` and writes
`output/models/twin.pt`. At least four real transition records are required.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from dt_airfoil.learning import train_twin
from dt_airfoil.replay import load_replay


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=600)
    parser.add_argument("--seed", type=int, default=7)
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parent
    replay = load_replay(root / "output" / "replay.jsonl")
    checkpoint = root / "output" / "models" / "twin.pt"
    report = train_twin(
        replay,
        checkpoint,
        epochs=arguments.epochs,
        seed=arguments.seed,
    )
    print(f"real transitions: {report.samples}")
    print(f"normalized training MSE: {report.training_mse:.6g}")
    if report.validation_mse is not None:
        print(f"normalized validation MSE: {report.validation_mse:.6g}")
    print(f"digital twin: {checkpoint}")


if __name__ == "__main__":
    main()

