#!/usr/bin/env python3

"""Run the optional continuous policy with online reward-model correction.

Command-line usage:
    python3 run_policy.py [--steps N] [--candidates N] [--max-retries N]
        [--model-epochs N] [--policy-epochs N] [--min-improvement FLOAT]
        [--seed N] [--animation-fps N] [--no-animation]

Requires `output/models/policy.pt` and `reward_model.pt`. It updates replay/model
data and writes trajectory CSV, GIF, and final PNG files below
`output/policy_rollout/`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ml_airfoil.environment import AirfoilMeshEnvironment
from ml_airfoil.learning import (
    PolicyPredictor,
    RewardPredictor,
    sample_actions,
    train_policy,
    train_reward_model,
)
from ml_airfoil.replay import (
    PolicyExample,
    TransitionRecord,
    append_policy_example,
    append_replay,
    load_policy_examples,
    load_replay,
)
from ml_airfoil.trajectory import TrajectoryRecorder, render_trajectory


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--candidates", type=int, default=256)
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--model-epochs", type=int, default=150)
    parser.add_argument("--policy-epochs", type=int, default=250)
    parser.add_argument(
        "--min-improvement",
        type=float,
        default=1.0e-4,
        help=(
            "reject smaller loss reductions and ask the corrected reward model "
            "for a stronger recovery action"
        ),
    )
    parser.add_argument("--seed", type=int, default=37)
    parser.add_argument("--animation-fps", type=int, default=2)
    parser.add_argument(
        "--no-animation",
        action="store_true",
        help="record the CSV but skip GIF and final PNG generation",
    )
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parent
    output = root / "output"
    policy_checkpoint = output / "models" / "policy.pt"
    reward_model_checkpoint = output / "models" / "reward_model.pt"
    replay_file = output / "replay.jsonl"
    policy_data_file = output / "policy_examples.jsonl"
    rollout_dir = output / "policy_rollout"
    trajectory_file = rollout_dir / "data_history.csv"
    animation_file = rollout_dir / "shape_evolution.gif"
    final_figure = rollout_dir / "shape_final.png"
    if not policy_checkpoint.exists():
        raise SystemExit("run optimize_with_reward_model.py before run_policy.py")
    if not reward_model_checkpoint.exists():
        raise SystemExit("run train_reward_model.py before run_policy.py")
    if arguments.candidates < 2:
        raise SystemExit("--candidates must be at least 2")
    if arguments.max_retries < 0:
        raise SystemExit("--max-retries must be nonnegative")
    if arguments.min_improvement < 0.0:
        raise SystemExit("--min-improvement must be nonnegative")

    rng = np.random.default_rng(arguments.seed)
    environment = AirfoilMeshEnvironment(root)
    environment.prepare()
    state = environment.reset()
    recorder = TrajectoryRecorder(trajectory_file)
    assert environment.current_loss is not None
    recorder.start(
        environment.working_dat,
        loss=environment.current_loss,
    )
    replay = load_replay(replay_file)
    policy_examples = load_policy_examples(policy_data_file)
    episode = 1 + max(
        (record.episode for record in replay),
        default=-1,
    )
    transition_index = 0
    print(f"initial loss = {environment.current_loss:.8e}")
    print(
        "minimum accepted improvement = "
        f"{arguments.min_improvement:.3e}"
    )
    for step in range(arguments.steps):
        policy = PolicyPredictor(policy_checkpoint)
        action = policy.predict(state)
        accepted = False

        for attempt in range(arguments.max_retries + 1):
            source = "policy" if attempt == 0 else "reward-model recovery"
            result = environment.step(
                action,
                accept_only_if_improves=True,
                minimum_improvement=arguments.min_improvement,
            )
            record = TransitionRecord.from_result(
                result,
                phase="policy",
                episode=episode,
                step=transition_index,
            )
            transition_index += 1
            append_replay(replay_file, record)
            replay.append(record)

            # Every real outcome, including a rejected or invalid action,
            # immediately corrects the reward model.
            report = train_reward_model(
                replay,
                reward_model_checkpoint,
                epochs=arguments.model_epochs,
                seed=arguments.seed + transition_index,
            )
            print(
                f"step {step}, attempt {attempt} ({source}): "
                f"{action.surface} x={action.center:.3f} "
                f"dy={action.shift:+.4f} "
                f"reward={result.reward:+.4e} "
                f"{'accepted' if result.accepted else 'rolled back'} "
                f"[model samples={report.samples}]"
            )
            if not result.valid_mesh:
                detail_lines = [
                    line.strip()
                    for line in result.message.splitlines()
                    if line.strip()
                ]
                detail = (
                    detail_lines[-1]
                    if detail_lines
                    else "mesh transition failed"
                )
                print(
                    "  invalid mesh: "
                    f"minimum quality={result.minimum_quality}; {detail}"
                )

            if result.accepted:
                accepted = True
                state = result.state_after_trial
                assert environment.current_loss is not None
                recorder.append(
                    environment.working_dat,
                    accepted_step=step + 1,
                    loss=environment.current_loss,
                    reward=result.reward,
                    source=source,
                    action=result.effective_action,
                )
                example = PolicyExample.create(
                    result.state_before,
                    result.effective_action,
                )
                append_policy_example(policy_data_file, example)
                policy_examples.append(example)
                if len(policy_examples) >= 4:
                    train_policy(
                        policy_examples,
                        policy_checkpoint,
                        epochs=arguments.policy_epochs,
                        seed=arguments.seed + 1000 + transition_index,
                    )
                break

            if attempt == arguments.max_retries:
                break

            # The state was rolled back and is therefore unchanged.  Use the
            # corrected reward model to search for a replacement action at that state.
            reward_model = RewardPredictor(reward_model_checkpoint)
            candidates = sample_actions(
                rng,
                arguments.candidates,
                # For this specific circle-to-NACA problem, recovery stays
                # inside the known thickness-reducing direction: upper down
                # and lower up.  The reward model still ranks locations/magnitudes,
                # and the real mesh still accepts or rejects the result.
                good_direction_probability=1.0,
            )
            predicted_rewards = reward_model.predict(state, candidates)
            action = candidates[int(np.argmax(predicted_rewards))]

        if not accepted:
            print(
                f"step {step}: no improving action after "
                f"{arguments.max_retries + 1} real evaluations; stopping"
            )
            break
    print(f"final loss = {environment.current_loss:.8e}")
    print(f"updated reward model:  {reward_model_checkpoint}")
    print(f"updated policy: {policy_checkpoint}")
    print(f"data history:   {trajectory_file}")
    if not arguments.no_animation:
        render_trajectory(
            trajectory_file,
            environment.target_dat,
            animation_file,
            final_figure,
            frames_per_second=arguments.animation_fps,
        )
        print(f"animation:      {animation_file}")
        print(f"final PNG:      {final_figure}")


if __name__ == "__main__":
    main()
