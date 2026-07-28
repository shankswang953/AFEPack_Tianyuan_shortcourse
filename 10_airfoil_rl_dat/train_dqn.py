#!/usr/bin/env python3

"""Train and evaluate Double DQN entirely in airfoil `dat` space.

Command-line usage:
    python3 train_dqn.py [--episodes N] [--max-steps N] [--seed N]
        [--batch-size N] [--warmup N] [--train-every N]
        [--terminal-loss FLOAT] [--fps N] [--evaluate-only]

Reads `data/target_naca0012.dat`. Training/evaluation artifacts are written
below `output/`, including the best checkpoint, configuration, training CSV
and plot, rollout CSV/`dat` history, final shape, loss plot, and animation.
`--evaluate-only` requires `output/checkpoints/dqn_best.pt`.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
from pathlib import Path

import numpy as np

from airfoil_rl.dqn import DoubleDQNAgent, ReplayBuffer
from airfoil_rl.environment import DatAirfoilEnvironment
from airfoil_rl.geometry import (
    Airfoil,
    make_circle_from_x_grid,
    read_airfoil,
    write_airfoil,
)
from airfoil_rl.plotting import (
    render_rollout,
    render_training_history,
    write_rollout_csv,
)


ROOT = Path(__file__).resolve().parent


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train and evaluate the pure-dat airfoil DQN example."
    )
    parser.add_argument("--episodes", type=int, default=400)
    parser.add_argument("--max-steps", type=int, default=140)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--warmup", type=int, default=1200)
    parser.add_argument("--train-every", type=int, default=2)
    parser.add_argument("--terminal-loss", type=float, default=2.0e-6)
    parser.add_argument("--fps", type=int, default=3)
    parser.add_argument(
        "--evaluate-only",
        action="store_true",
        help="load output/checkpoints/dqn_best.pt and render one rollout",
    )
    arguments = parser.parse_args()
    if arguments.episodes < 1:
        raise SystemExit("--episodes must be positive")
    if arguments.max_steps < 1:
        raise SystemExit("--max-steps must be positive")
    return arguments


def epsilon_at(episode: int, total_episodes: int) -> float:
    decay_end = max(1, int(0.75 * total_episodes))
    fraction = min(1.0, episode / decay_end)
    return 1.0 + fraction * (0.05 - 1.0)


def evaluate(
    agent: DoubleDQNAgent,
    environment: DatAirfoilEnvironment,
) -> tuple[list[Airfoil], list, list[float], list[float]]:
    state = environment.reset(evaluation=True)
    states = [environment.current.copy()]
    actions = []
    losses = [environment.loss]
    rewards = []
    for _ in range(environment.max_steps):
        action_index = agent.greedy_action(
            state,
            environment.valid_action_mask(),
        )
        result = environment.step(action_index)
        state = result.state
        states.append(environment.current.copy())
        actions.append(result.action)
        losses.append(result.loss)
        rewards.append(result.reward)
        if result.terminated or result.truncated:
            break
    return states, actions, losses, rewards


def write_training_csv(
    filename: Path,
    rows: list[dict[str, float | int]],
) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(rows[0]),
        )
        writer.writeheader()
        writer.writerows(rows)


def train(
    *,
    environment: DatAirfoilEnvironment,
    episodes: int,
    seed: int,
    batch_size: int,
    warmup: int,
    train_every: int,
    checkpoint: Path,
) -> tuple[DoubleDQNAgent, list[dict[str, float | int]]]:
    agent = DoubleDQNAgent(
        environment.state_dimension,
        environment.action_count,
        seed=seed,
    )
    replay = ReplayBuffer(capacity=100_000, seed=seed + 10)
    history: list[dict[str, float | int]] = []
    best_loss = math.inf
    best_state = copy.deepcopy(agent.online.state_dict())
    total_steps = 0

    for episode in range(episodes):
        state = environment.reset()
        epsilon = epsilon_at(episode, episodes)
        episode_return = 0.0
        losses = []
        for _ in range(environment.max_steps):
            action_index = agent.select_action(
                state,
                epsilon,
                environment.valid_action_mask(),
            )
            result = environment.step(action_index)
            done = result.terminated or result.truncated
            replay.append(
                state,
                action_index,
                result.reward,
                result.state,
                done,
                environment.valid_action_mask(),
            )
            state = result.state
            episode_return += result.reward
            total_steps += 1
            if (
                len(replay) >= max(warmup, batch_size)
                and total_steps % train_every == 0
            ):
                losses.append(agent.optimize(replay.sample(batch_size)))
            if done:
                break

        row: dict[str, float | int] = {
            "episode": episode + 1,
            "steps": environment.step_count,
            "return": episode_return,
            "final_loss": environment.loss,
            "epsilon": epsilon,
            "mean_td_loss": (
                float(np.mean(losses)) if losses else float("nan")
            ),
        }
        history.append(row)

        if (episode + 1) % 10 == 0 or episode == episodes - 1:
            _, _, evaluation_losses, _ = evaluate(agent, environment)
            evaluation_loss = evaluation_losses[-1]
            if evaluation_loss < best_loss:
                best_loss = evaluation_loss
                best_state = copy.deepcopy(agent.online.state_dict())
            print(
                f"episode {episode + 1:4d}/{episodes}: "
                f"train loss={float(row['final_loss']):.3e}, "
                f"eval loss={evaluation_loss:.3e}, "
                f"epsilon={epsilon:.3f}"
            )

    agent.online.load_state_dict(best_state)
    agent.target.load_state_dict(best_state)
    agent.save(
        checkpoint,
        metadata={
            "seed": seed,
            "episodes": episodes,
            "best_evaluation_loss": best_loss,
        },
    )
    return agent, history


def main() -> None:
    arguments = parse_arguments()
    data = ROOT / "data"
    output = ROOT / "output"
    target = read_airfoil(data / "target_naca0012.dat")
    circle = make_circle_from_x_grid(target)
    write_airfoil(data / "initial_circle.dat", circle)
    environment = DatAirfoilEnvironment(
        target,
        circle,
        max_steps=arguments.max_steps,
        terminal_loss=arguments.terminal_loss,
        seed=arguments.seed,
    )
    checkpoint = output / "checkpoints" / "dqn_best.pt"

    history: list[dict[str, float | int]] = []
    if arguments.evaluate_only:
        if not checkpoint.exists():
            raise SystemExit(
                f"checkpoint does not exist: {checkpoint}; train first"
            )
        agent = DoubleDQNAgent.load(checkpoint)
    else:
        agent, history = train(
            environment=environment,
            episodes=arguments.episodes,
            seed=arguments.seed,
            batch_size=arguments.batch_size,
            warmup=arguments.warmup,
            train_every=arguments.train_every,
            checkpoint=checkpoint,
        )
        write_training_csv(output / "training.csv", history)
        render_training_history(
            history,
            output / "training_history.png",
        )

    states, actions, losses, rewards = evaluate(agent, environment)
    evaluation = output / "evaluation"
    dat_history = evaluation / "dat_history"
    if dat_history.exists():
        for old_state in dat_history.glob("step_*.dat"):
            old_state.unlink()
    for step, state in enumerate(states):
        write_airfoil(dat_history / f"step_{step:04d}.dat", state)
    write_airfoil(evaluation / "final.dat", states[-1])
    write_rollout_csv(
        evaluation / "rollout.csv",
        actions,
        losses,
        rewards,
    )
    render_rollout(
        states=states,
        actions=actions,
        losses=losses,
        target=target,
        animation_file=evaluation / "shape_evolution.gif",
        final_figure=evaluation / "final_shape.png",
        loss_figure=evaluation / "loss_history.png",
        fps=arguments.fps,
    )
    configuration = {
        "seed": arguments.seed,
        "episodes": arguments.episodes,
        "max_steps": arguments.max_steps,
        "terminal_loss": arguments.terminal_loss,
        "state_dimension": environment.state_dimension,
        "action_count": environment.action_count,
        "initial_loss": losses[0],
        "final_loss": losses[-1],
        "evaluation_steps": len(actions),
        "terminated": losses[-1] <= arguments.terminal_loss,
        "mesh_used": False,
    }
    (output / "configuration.json").write_text(
        json.dumps(configuration, indent=2)
    )
    print(f"evaluation steps: {len(actions)}")
    print(f"initial loss: {losses[0]:.8e}")
    print(f"final loss:   {losses[-1]:.8e}")
    print(f"checkpoint:   {checkpoint}")
    print(f"animation:    {evaluation / 'shape_evolution.gif'}")


if __name__ == "__main__":
    main()
