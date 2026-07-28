#!/usr/bin/env python3

"""Run the robust model-based controller and verify actions on the real mesh.

Command-line usage:
    python3 run_twin_controller.py [--steps N] [--center-count N]
        [--max-real-trials N] [--twin-epochs N] [--seed N] [OTHER OPTIONS]

Run with `--help` for improvement, quality, exploration, safeguard,
backtracking, remeshing, and animation controls. Requires
`output/models/twin.pt`; updates replay/model state and writes the accepted
trajectory and figures below `output/policy_rollout/`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from dt_airfoil.environment import AirfoilMeshEnvironment
from dt_airfoil.geometry import (
    Action,
    action_points_toward_target,
    apply_action,
    data_shape_mse,
    read_airfoil,
)
from dt_airfoil.learning import (
    TwinPredictor,
    discrete_actions,
    train_twin,
)
from dt_airfoil.replay import (
    TransitionRecord,
    append_replay,
    load_replay,
)
from dt_airfoil.trajectory import TrajectoryRecorder, render_trajectory


def action_key(action: Action) -> tuple[str, float, float]:
    return (
        action.surface,
        round(action.center, 10),
        round(action.shift, 10),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=80)
    parser.add_argument("--center-count", type=int, default=21)
    parser.add_argument("--max-real-trials", type=int, default=10)
    parser.add_argument("--twin-epochs", type=int, default=150)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--min-improvement", type=float, default=1.0e-8)
    parser.add_argument("--minimum-mesh-quality", type=float, default=0.40)
    parser.add_argument("--exploration-bonus", type=float, default=1.0e-3)
    parser.add_argument("--twin-safeguard-ratio", type=float, default=0.25)
    parser.add_argument("--fine-action-threshold", type=float, default=0.01)
    parser.add_argument("--min-backtrack-shift", type=float, default=0.003)
    parser.add_argument(
        "--no-remesh",
        action="store_true",
        help="disable the EasyMesh fallback and keep topology fixed",
    )
    parser.add_argument("--animation-fps", type=int, default=2)
    parser.add_argument("--no-animation", action="store_true")
    arguments = parser.parse_args()

    if arguments.steps < 1:
        raise SystemExit("--steps must be positive")
    if arguments.center_count < 2:
        raise SystemExit("--center-count must be at least two")
    if arguments.max_real_trials < 1:
        raise SystemExit("--max-real-trials must be positive")
    if arguments.min_improvement < 0.0:
        raise SystemExit("--min-improvement must be nonnegative")
    if arguments.minimum_mesh_quality < 0.0:
        raise SystemExit("--minimum-mesh-quality must be nonnegative")
    if arguments.exploration_bonus < 0.0:
        raise SystemExit("--exploration-bonus must be nonnegative")
    if not 0.0 <= arguments.twin_safeguard_ratio <= 1.0:
        raise SystemExit("--twin-safeguard-ratio must lie in [0, 1]")
    if arguments.fine_action_threshold < 0.0:
        raise SystemExit("--fine-action-threshold must be nonnegative")
    if arguments.min_backtrack_shift <= 0.0:
        raise SystemExit("--min-backtrack-shift must be positive")

    root = Path(__file__).resolve().parent
    output = root / "output"
    replay_file = output / "replay.jsonl"
    twin_checkpoint = output / "models" / "twin.pt"
    rollout_dir = output / "policy_rollout"
    history_file = rollout_dir / "data_history.csv"
    animation_file = rollout_dir / "shape_evolution.gif"
    final_figure = rollout_dir / "shape_final.png"
    if not twin_checkpoint.exists():
        raise SystemExit("run train_twin.py before the twin controller")

    environment = AirfoilMeshEnvironment(root)
    environment.prepare()
    state = environment.reset()
    assert environment.current_loss is not None
    target_airfoil = read_airfoil(environment.target_dat)

    recorder = TrajectoryRecorder(history_file)
    recorder.start(environment.working_dat, loss=environment.current_loss)
    replay = load_replay(replay_file)
    episode = 1 + max(
        (record.episode for record in replay),
        default=-1,
    )
    transition_index = 0
    base_actions = discrete_actions(
        center_count=arguments.center_count,
    )
    # An invalid action lowers the allowed magnitude at that surface/center
    # for the remainder of this rollout.
    local_shift_limits: dict[tuple[str, float, int], float] = {}
    center_visits: dict[tuple[str, float, int], int] = {}

    print(f"initial loss = {environment.current_loss:.8e}")
    print(
        f"action bank = {len(base_actions)} candidates over "
        f"{arguments.center_count} full-chord centers"
    )
    print(
        "minimum accepted improvement = "
        f"{arguments.min_improvement:.3e}"
    )
    print(
        "EasyMesh fallback = "
        f"{'disabled' if arguments.no_remesh else 'enabled'}"
    )
    print(f"controller seed = {arguments.seed}")

    for step in range(arguments.steps):
        tried: set[tuple[str, float, float]] = set()
        accepted = False
        objective_converged = False
        next_action: Action | None = None
        source = "twin grid"

        for trial in range(arguments.max_real_trials):
            if next_action is None:
                allowed_surfaces = {"U", "L"}
                candidates = []
                candidate_keys: set[tuple[str, float, float]] = set()
                current_airfoil = read_airfoil(environment.working_dat)
                for candidate in base_actions:
                    if candidate.surface not in allowed_surfaces:
                        continue
                    if (
                        environment.current_loss
                        > arguments.fine_action_threshold
                        and abs(candidate.shift) < 0.02
                    ):
                        continue
                    if not action_points_toward_target(
                        current_airfoil,
                        target_airfoil,
                        candidate,
                    ):
                        continue
                    _, effective_candidate = apply_action(
                        current_airfoil,
                        candidate,
                        width=environment.gaussian_width,
                        shift_max=environment.shift_max,
                        thickness_reference=target_airfoil,
                    )
                    if abs(effective_candidate.shift) < 1.0e-12:
                        continue
                    effective_key = action_key(effective_candidate)
                    if effective_key in candidate_keys:
                        continue
                    direction = (
                        1 if effective_candidate.shift > 0.0 else -1
                    )
                    local_key = (
                        effective_candidate.surface,
                        round(effective_candidate.center, 10),
                        direction,
                    )
                    local_limit = local_shift_limits.get(
                        local_key,
                        float("inf"),
                    )
                    if abs(effective_candidate.shift) > local_limit:
                        continue
                    if effective_key in tried:
                        continue
                    candidates.append(effective_candidate)
                    candidate_keys.add(effective_key)
                if not candidates:
                    break
                # The dat objective is intentionally cheap in this teaching
                # problem.  Evaluate it for every candidate as a safety oracle
                # before spending a real AFEPack mesh transition.
                exact_rewards = []
                for candidate in candidates:
                    moved, _ = apply_action(
                        current_airfoil,
                        candidate,
                        width=environment.gaussian_width,
                        shift_max=environment.shift_max,
                        thickness_reference=target_airfoil,
                    )
                    exact_rewards.append(
                        environment.current_loss
                        - data_shape_mse(moved, target_airfoil)
                    )
                best_index = int(np.argmax(exact_rewards))
                best_exact_reward = exact_rewards[best_index]
                if best_exact_reward <= arguments.min_improvement:
                    objective_converged = True
                    break
                if trial == 0:
                    # The learned digital twin always makes the first
                    # proposal.  Coverage discourages repeatedly selecting
                    # the same center while the model is still uncertain.
                    twin = TwinPredictor(twin_checkpoint)
                    predicted_rewards = twin.predict(state, candidates)
                    coverage_bonus = np.asarray(
                        [
                            arguments.exploration_bonus
                            / np.sqrt(
                                1.0
                                + center_visits.get(
                                    (
                                        candidate.surface,
                                        round(candidate.center, 10),
                                        (
                                            1
                                            if candidate.shift > 0.0
                                            else -1
                                        ),
                                    ),
                                    0,
                                )
                            )
                            for candidate in candidates
                        ]
                    )
                    acquisition_score = predicted_rewards + coverage_bonus
                    twin_index = int(np.argmax(acquisition_score))
                    if (
                        exact_rewards[twin_index]
                        >= arguments.twin_safeguard_ratio
                        * best_exact_reward
                    ):
                        next_action = candidates[twin_index]
                        source = "twin grid"
                    else:
                        next_action = candidates[best_index]
                        source = "data-MSE safeguard"
                else:
                    # A rejected NN proposal is useful training data, but
                    # retraining need not fix the very next ranking.  The dat
                    # objective is cheap to evaluate before moving the mesh,
                    # so use it to rank deterministic recovery candidates.
                    next_action = candidates[best_index]
                    source = "data-MSE recovery"

            action = next_action
            next_action = None
            result = environment.step(
                action,
                accept_only_if_improves=True,
                minimum_improvement=arguments.min_improvement,
                minimum_mesh_quality=arguments.minimum_mesh_quality,
                remesh_if_needed=not arguments.no_remesh,
            )
            tried.add(action_key(action))
            tried.add(action_key(result.effective_action))
            record = TransitionRecord.from_result(
                result,
                phase="twin_controller",
                episode=episode,
                step=transition_index,
            )
            transition_index += 1
            append_replay(replay_file, record)
            replay.append(record)
            report = train_twin(
                replay,
                twin_checkpoint,
                epochs=arguments.twin_epochs,
                seed=arguments.seed + transition_index,
            )

            print(
                f"step {step}, trial {trial} ({source}): "
                f"{action.surface} x={action.center:.3f} "
                f"dy={action.shift:+.4f} "
                f"reward={result.reward:+.4e} "
                f"{'accepted' if result.accepted else 'rolled back'} "
                f"{'[remeshed] ' if result.remeshed else ''}"
                f"[twin samples={report.samples}]"
            )

            if result.accepted:
                accepted = True
                state = result.state_after_trial
                visit_key = (
                    result.effective_action.surface,
                    round(result.effective_action.center, 10),
                    1 if result.effective_action.shift > 0.0 else -1,
                )
                center_visits[visit_key] = (
                    center_visits.get(visit_key, 0) + 1
                )
                if result.remeshed:
                    # Limits inferred from the old topology no longer apply.
                    local_shift_limits.clear()
                assert environment.current_loss is not None
                recorder.append(
                    environment.working_dat,
                    accepted_step=step + 1,
                    loss=environment.current_loss,
                    reward=result.reward,
                    source=(
                        f"{source} + EasyMesh remesh"
                        if result.remeshed
                        else source
                    ),
                    action=result.effective_action,
                )
                break

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
                reduced_shift = 0.5 * result.effective_action.shift
                local_key = (
                    result.effective_action.surface,
                    round(result.effective_action.center, 10),
                    1 if result.effective_action.shift > 0.0 else -1,
                )
                local_shift_limits[local_key] = min(
                    local_shift_limits.get(local_key, float("inf")),
                    abs(reduced_shift),
                )
                if abs(reduced_shift) >= arguments.min_backtrack_shift:
                    next_action = Action(
                        result.effective_action.surface,
                        result.effective_action.center,
                        reduced_shift,
                    )
                    source = "mesh backtracking"

        if not accepted:
            assert environment.current_loss is not None
            if (
                objective_converged
                or environment.current_loss <= arguments.min_improvement
            ):
                print(
                    f"step {step}: converged at the current action "
                    f"resolution (loss={environment.current_loss:.3e})"
                )
            else:
                print(
                    f"step {step}: no acceptable action after "
                    f"{arguments.max_real_trials} real trials; stopping"
                )
            break

    assert environment.current_loss is not None
    print(f"final loss = {environment.current_loss:.8e}")
    print(f"updated twin: {twin_checkpoint}")
    print(f"data history: {history_file}")
    if not arguments.no_animation:
        render_trajectory(
            history_file,
            environment.target_dat,
            animation_file,
            final_figure,
            frames_per_second=arguments.animation_fps,
        )
        print(f"animation:    {animation_file}")
        print(f"final PNG:    {final_figure}")


if __name__ == "__main__":
    main()
