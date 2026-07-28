#!/usr/bin/env python3

"""Rank candidate actions with the twin and verify them on the real mesh.

Command-line usage:
    python3 optimize_with_twin.py [--episodes N] [--steps N]
        [--candidates N] [--seed N] [--retrain-epochs N]

Requires `output/models/twin.pt`. It updates `output/replay.jsonl`,
`output/policy_examples.jsonl`, and the twin/policy checkpoints below
`output/models/`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from dt_airfoil.environment import AirfoilMeshEnvironment
from dt_airfoil.learning import (
    PolicyPredictor,
    TwinPredictor,
    sample_actions,
    train_policy,
    train_twin,
)
from dt_airfoil.replay import (
    PolicyExample,
    TransitionRecord,
    append_policy_example,
    append_replay,
    load_policy_examples,
    load_replay,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--episodes", type=int, default=4)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--candidates", type=int, default=256)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--retrain-epochs", type=int, default=250)
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parent
    output = root / "output"
    replay_file = output / "replay.jsonl"
    policy_data_file = output / "policy_examples.jsonl"
    twin_checkpoint = output / "models" / "twin.pt"
    policy_checkpoint = output / "models" / "policy.pt"
    if not twin_checkpoint.exists():
        raise SystemExit("run train_twin.py before twin optimization")

    rng = np.random.default_rng(arguments.seed)
    environment = AirfoilMeshEnvironment(root)
    environment.prepare()
    replay = load_replay(replay_file)
    episode_offset = (
        1 + max((record.episode for record in replay), default=-1)
    )

    for local_episode in range(arguments.episodes):
        episode = episode_offset + local_episode
        state = environment.reset()
        twin = TwinPredictor(twin_checkpoint)
        policy = (
            PolicyPredictor(policy_checkpoint)
            if policy_checkpoint.exists()
            else None
        )
        print(
            f"twin episode {episode}: initial loss = "
            f"{environment.current_loss:.8e}"
        )
        for step in range(arguments.steps):
            candidates = sample_actions(
                rng,
                arguments.candidates,
                good_direction_probability=0.50,
            )
            if policy is not None:
                candidates.append(policy.predict(state))
            predicted_rewards = twin.predict(state, candidates)
            selected = candidates[int(np.argmax(predicted_rewards))]
            result = environment.step(
                selected,
                accept_only_if_improves=True,
            )
            record = TransitionRecord.from_result(
                result,
                phase="twin",
                episode=episode,
                step=step,
            )
            append_replay(replay_file, record)
            replay.append(record)
            if result.accepted:
                example = PolicyExample.create(state, result.effective_action)
                append_policy_example(policy_data_file, example)
                state = result.state_after_trial
            print(
                f"  step {step}: predicted={np.max(predicted_rewards):+.4e} "
                f"real={result.reward:+.4e} "
                f"{'accepted' if result.accepted else 'rolled back'}"
            )

        train_twin(
            replay,
            twin_checkpoint,
            epochs=arguments.retrain_epochs,
            seed=arguments.seed + local_episode + 1,
        )
        policy_examples = load_policy_examples(policy_data_file)
        if len(policy_examples) >= 4:
            train_policy(
                policy_examples,
                policy_checkpoint,
                epochs=400,
                seed=arguments.seed + local_episode + 101,
            )
        print(f"  episode final loss = {environment.current_loss:.8e}")

    print(f"updated twin:   {twin_checkpoint}")
    if policy_checkpoint.exists():
        print(f"direct policy:  {policy_checkpoint}")


if __name__ == "__main__":
    main()

