#!/usr/bin/env python3

"""Collect real-mesh transitions for initial digital-twin training.

Command-line usage:
    python3 collect_random.py [--episodes N] [--steps N] [--seed N] [--fresh]

The script prepares/resets the mesh environment, executes random actions, and
appends JSON records to `output/replay.jsonl`; `--fresh` removes the prior
replay file first. EasyMesh paths are selected through EASYMESH_BIN and
EASYMESH2MESH_BIN.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from dt_airfoil.environment import AirfoilMeshEnvironment
from dt_airfoil.learning import sample_actions
from dt_airfoil.replay import (
    TransitionRecord,
    append_replay,
    load_replay,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--episodes", type=int, default=6)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="discard the existing replay file before collection",
    )
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parent
    replay_file = root / "output" / "replay.jsonl"
    if arguments.fresh and replay_file.exists():
        replay_file.unlink()
    existing = load_replay(replay_file)
    episode_offset = (
        1 + max((record.episode for record in existing), default=-1)
    )

    rng = np.random.default_rng(arguments.seed)
    environment = AirfoilMeshEnvironment(root)
    environment.prepare()
    for local_episode in range(arguments.episodes):
        episode = episode_offset + local_episode
        environment.reset()
        print(
            f"episode {episode}: initial loss = "
            f"{environment.current_loss:.8e}"
        )
        for step in range(arguments.steps):
            action = sample_actions(rng, 1)[0]
            result = environment.step(action)
            record = TransitionRecord.from_result(
                result,
                phase="warmup",
                episode=episode,
                step=step,
            )
            append_replay(replay_file, record)
            status = (
                "valid"
                if result.valid_mesh
                else "invalid and rolled back"
            )
            print(
                f"  step {step}: {action.surface} "
                f"x={action.center:.3f} dy={action.shift:+.4f} "
                f"reward={result.reward:+.4e} {status}"
            )
    print(f"replay file: {replay_file}")


if __name__ == "__main__":
    main()

