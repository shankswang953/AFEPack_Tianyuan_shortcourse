#!/usr/bin/env python3

"""Run the complete five-stage airfoil digital-twin experiment.

Command-line usage:
    python3 run_experiment.py [--reset] [--seed N] [STAGE OPTIONS]

Use `python3 run_experiment.py --help` for warm-up, training, optimization,
controller, mesh-quality, and animation parameters. Results are kept below
`output/`, including `run_configuration.json`, replay data, models, mesh
history, and the policy rollout. `--clean-build` is valid only with
`--reset`. External mesh tools are configured with EASYMESH_BIN and
EASYMESH2MESH_BIN.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import torch

from reset_project import reset_generated_state


def run_stage(root: Path, title: str, arguments: list[str]) -> None:
    command = [sys.executable, *arguments]
    print()
    print("=" * 72)
    print(title)
    print("$ " + " ".join(command))
    print("=" * 72, flush=True)
    start = time.monotonic()
    subprocess.run(command, cwd=root, check=True)
    elapsed = time.monotonic() - start
    print(f"{title} completed in {elapsed:.1f} s", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reset",
        action="store_true",
        help="discard all previous meshes, replay data, models, and histories",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2026,
        help="master seed used to derive every stochastic stage seed",
    )
    parser.add_argument(
        "--clean-build",
        action="store_true",
        help="with --reset, also rebuild the C++ executables",
    )
    parser.add_argument("--warmup-episodes", type=int, default=6)
    parser.add_argument("--warmup-steps", type=int, default=5)
    parser.add_argument("--twin-epochs", type=int, default=600)
    parser.add_argument(
        "--controller-twin-epochs",
        type=int,
        default=30,
        help="online twin-training epochs after each real controller trial",
    )
    parser.add_argument("--optimization-episodes", type=int, default=4)
    parser.add_argument("--optimization-steps", type=int, default=5)
    parser.add_argument(
        "--controller-steps",
        "--policy-steps",
        dest="controller_steps",
        type=int,
        default=20,
    )
    parser.add_argument("--candidates", type=int, default=256)
    parser.add_argument(
        "--max-real-trials",
        "--max-retries",
        dest="max_real_trials",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--min-improvement",
        type=float,
        default=1.0e-8,
    )
    parser.add_argument(
        "--minimum-mesh-quality",
        type=float,
        default=0.40,
    )
    parser.add_argument(
        "--twin-safeguard-ratio",
        type=float,
        default=0.25,
    )
    parser.add_argument(
        "--exploration-bonus",
        type=float,
        default=1.0e-3,
    )
    parser.add_argument(
        "--fine-action-threshold",
        type=float,
        default=0.01,
    )
    parser.add_argument("--animation-fps", type=int, default=2)
    parser.add_argument(
        "--no-animation",
        action="store_true",
        help="record the final rollout but do not render its GIF",
    )
    arguments = parser.parse_args()

    if arguments.clean_build and not arguments.reset:
        parser.error("--clean-build requires --reset")

    root = Path(__file__).resolve().parent
    total_start = time.monotonic()
    stage_seeds = {
        "warmup": arguments.seed + 101,
        "initial_twin": arguments.seed + 201,
        "twin_optimization": arguments.seed + 301,
        "controller": arguments.seed + 401,
    }

    if arguments.reset:
        print("Resetting all reproducible experiment state ...")
        reset_generated_state(
            root,
            clean_build=arguments.clean_build,
        )
    else:
        print(
            "Continuing existing replay/model state: the seed is fixed, "
            "but exact clean-run reproduction requires --reset."
        )

    output = root / "output"
    output.mkdir(parents=True, exist_ok=True)
    run_configuration = {
        "arguments": vars(arguments),
        "stage_seeds": stage_seeds,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "torch": torch.__version__,
        },
        "reproducibility_scope": (
            "Exact action sequence on the same machine and software "
            "environment when run with --reset."
        ),
    }
    configuration_file = output / "run_configuration.json"
    configuration_file.write_text(
        json.dumps(run_configuration, indent=2, sort_keys=True) + "\n"
    )
    print(f"master seed = {arguments.seed}")
    print(f"run configuration: {configuration_file}")

    run_stage(
        root,
        "1/5 Prepare AFEPack and EasyMesh reference meshes",
        ["prepare_case.py"],
    )
    warmup_arguments = [
        "collect_random.py",
        "--episodes",
        str(arguments.warmup_episodes),
        "--steps",
        str(arguments.warmup_steps),
        "--seed",
        str(stage_seeds["warmup"]),
    ]
    if arguments.reset:
        warmup_arguments.append("--fresh")
    run_stage(
        root,
        "2/5 Collect real random mesh transitions",
        warmup_arguments,
    )
    run_stage(
        root,
        "3/5 Train the reward-predicting digital twin",
        [
            "train_twin.py",
            "--epochs",
            str(arguments.twin_epochs),
            "--seed",
            str(stage_seeds["initial_twin"]),
        ],
    )
    run_stage(
        root,
        "4/5 Optimize with twin-ranked, real-verified actions",
        [
            "optimize_with_twin.py",
            "--episodes",
            str(arguments.optimization_episodes),
            "--steps",
            str(arguments.optimization_steps),
            "--candidates",
            str(arguments.candidates),
            "--seed",
            str(stage_seeds["twin_optimization"]),
        ],
    )
    policy_arguments = [
        "run_twin_controller.py",
        "--steps",
        str(arguments.controller_steps),
        "--center-count",
        "21",
        "--max-real-trials",
        str(arguments.max_real_trials),
        "--twin-epochs",
        str(arguments.controller_twin_epochs),
        "--seed",
        str(stage_seeds["controller"]),
        "--min-improvement",
        str(arguments.min_improvement),
        "--minimum-mesh-quality",
        str(arguments.minimum_mesh_quality),
        "--exploration-bonus",
        str(arguments.exploration_bonus),
        "--twin-safeguard-ratio",
        str(arguments.twin_safeguard_ratio),
        "--fine-action-threshold",
        str(arguments.fine_action_threshold),
        "--animation-fps",
        str(arguments.animation_fps),
    ]
    if arguments.no_animation:
        policy_arguments.append("--no-animation")
    run_stage(
        root,
        "5/5 Run the model-based twin controller and record shape history",
        policy_arguments,
    )

    elapsed = time.monotonic() - total_start
    rollout = root / "output" / "policy_rollout"
    print()
    print(f"Complete experiment finished in {elapsed:.1f} s")
    print(f"trajectory data: {rollout / 'data_history.csv'}")
    if not arguments.no_animation:
        print(f"animation:       {rollout / 'shape_evolution.gif'}")
        print(f"final figure:    {rollout / 'shape_final.png'}")


if __name__ == "__main__":
    main()
