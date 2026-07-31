#!/usr/bin/env python3

"""Train the one-step reward model from recorded real transitions.

Command-line usage:
    python3 train_reward_model.py [--epochs N] [--seed N]

Reads `output/replay.jsonl` and writes
`output/models/reward_model.pt`. At least four real transition records are
required.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ml_airfoil.learning import train_reward_model
from ml_airfoil.replay import load_replay


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=600)
    parser.add_argument("--seed", type=int, default=7)
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parent
    replay = load_replay(root / "output" / "replay.jsonl")
    checkpoint = root / "output" / "models" / "reward_model.pt"
    report = train_reward_model(
        replay,
        checkpoint,
        epochs=arguments.epochs,
        seed=arguments.seed,
    )
    print(f"real transitions: {report.samples}")
    print(f"normalized training MSE: {report.training_mse:.6g}")
    if report.validation_mse is not None:
        print(f"normalized validation MSE: {report.validation_mse:.6g}")
    print(f"reward model: {checkpoint}")


if __name__ == "__main__":
    main()
