"""Provide deterministic training and inference helpers for twin/policy models.

This import-only library has no command-line interface. Training functions
write checkpoints to caller-supplied paths (normally `output/models/*.pt`);
prediction helpers load those checkpoints. NumPy and PyTorch are required.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch import nn

from .geometry import Action, action_vector
from .models import PolicyNetwork, TwinNetwork
from .replay import PolicyExample, TransitionRecord


@dataclass(frozen=True)
class TrainingReport:
    samples: int
    training_mse: float
    validation_mse: float | None


def configure_reproducible_torch(seed: int) -> None:
    """Make the small CPU networks deterministic for a fixed environment."""

    torch.manual_seed(seed)
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)


def _state_statistics(states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.mean(states, axis=0)
    standard_deviation = np.std(states, axis=0)
    standard_deviation[standard_deviation < 1.0e-6] = 1.0
    return mean.astype(np.float32), standard_deviation.astype(np.float32)


def train_twin(
    records: list[TransitionRecord],
    checkpoint: Path,
    *,
    epochs: int = 600,
    seed: int = 7,
) -> TrainingReport:
    if len(records) < 4:
        raise ValueError("at least four real transitions are needed")
    configure_reproducible_torch(seed)
    rng = np.random.default_rng(seed)

    states = np.stack([record.state for record in records])
    actions = np.stack(
        [action_vector(record.action) for record in records]
    )
    rewards = np.asarray(
        [record.reward for record in records],
        dtype=np.float32,
    )
    state_mean, state_std = _state_statistics(states)
    reward_mean = float(np.mean(rewards))
    reward_std = float(np.std(rewards))
    if reward_std < 1.0e-7:
        reward_std = 1.0

    normalized_states = (states - state_mean) / state_std
    normalized_rewards = (rewards - reward_mean) / reward_std
    permutation = rng.permutation(len(records))
    validation_count = (
        max(1, len(records) // 5)
        if len(records) >= 10
        else 0
    )
    validation_indices = permutation[:validation_count]
    training_indices = permutation[validation_count:]

    state_tensor = torch.from_numpy(normalized_states)
    action_tensor = torch.from_numpy(actions)
    reward_tensor = torch.from_numpy(normalized_rewards)
    model = TwinNetwork(states.shape[1])
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=2.0e-3,
        weight_decay=1.0e-4,
    )
    loss_function = nn.MSELoss()

    train_index_tensor = torch.from_numpy(training_indices)
    for _ in range(epochs):
        model.train()
        prediction = model(
            state_tensor[train_index_tensor],
            action_tensor[train_index_tensor],
        )
        loss = loss_function(
            prediction,
            reward_tensor[train_index_tensor],
        )
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    model.eval()
    with torch.no_grad():
        training_prediction = model(
            state_tensor[train_index_tensor],
            action_tensor[train_index_tensor],
        )
        training_mse = float(
            loss_function(
                training_prediction,
                reward_tensor[train_index_tensor],
            )
        )
        validation_mse = None
        if validation_count:
            validation_tensor = torch.from_numpy(validation_indices)
            validation_prediction = model(
                state_tensor[validation_tensor],
                action_tensor[validation_tensor],
            )
            validation_mse = float(
                loss_function(
                    validation_prediction,
                    reward_tensor[validation_tensor],
                )
            )

    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state_dimension": states.shape[1],
            "model": model.state_dict(),
            "state_mean": state_mean,
            "state_std": state_std,
            "reward_mean": reward_mean,
            "reward_std": reward_std,
        },
        checkpoint,
    )
    return TrainingReport(
        samples=len(records),
        training_mse=training_mse,
        validation_mse=validation_mse,
    )


class TwinPredictor:
    def __init__(self, checkpoint: Path) -> None:
        bundle = torch.load(
            checkpoint,
            map_location="cpu",
            weights_only=False,
        )
        self.state_mean = np.asarray(
            bundle["state_mean"],
            dtype=np.float32,
        )
        self.state_std = np.asarray(
            bundle["state_std"],
            dtype=np.float32,
        )
        self.reward_mean = float(bundle["reward_mean"])
        self.reward_std = float(bundle["reward_std"])
        self.model = TwinNetwork(int(bundle["state_dimension"]))
        self.model.load_state_dict(bundle["model"])
        self.model.eval()

    def predict(
        self,
        state: np.ndarray,
        actions: list[Action],
    ) -> np.ndarray:
        repeated = np.repeat(
            ((state - self.state_mean) / self.state_std)[None, :],
            len(actions),
            axis=0,
        ).astype(np.float32)
        action_array = np.stack([action_vector(action) for action in actions])
        with torch.no_grad():
            normalized = self.model(
                torch.from_numpy(repeated),
                torch.from_numpy(action_array),
            ).numpy()
        return normalized * self.reward_std + self.reward_mean


def train_policy(
    examples: list[PolicyExample],
    checkpoint: Path,
    *,
    epochs: int = 800,
    seed: int = 11,
) -> TrainingReport:
    if len(examples) < 4:
        raise ValueError("at least four accepted twin actions are needed")
    configure_reproducible_torch(seed)
    states = np.asarray(
        [example.state for example in examples],
        dtype=np.float32,
    )
    state_mean, state_std = _state_statistics(states)
    normalized_states = (states - state_mean) / state_std
    surfaces = np.asarray(
        [
            0 if str(example.action["surface"]) == "U" else 1
            for example in examples
        ],
        dtype=np.int64,
    )
    centers = np.asarray(
        [float(example.action["center"]) for example in examples],
        dtype=np.float32,
    )
    shifts = np.asarray(
        [float(example.action["shift"]) for example in examples],
        dtype=np.float32,
    )

    model = PolicyNetwork(states.shape[1])
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=2.0e-3,
        weight_decay=1.0e-4,
    )
    state_tensor = torch.from_numpy(normalized_states)
    surface_tensor = torch.from_numpy(surfaces)
    center_tensor = torch.from_numpy(centers)
    shift_tensor = torch.from_numpy(shifts)

    for _ in range(epochs):
        logits, center_prediction, shift_prediction = model(state_tensor)
        loss = (
            nn.functional.cross_entropy(logits, surface_tensor)
            + 4.0 * nn.functional.mse_loss(
                center_prediction,
                center_tensor,
            )
            + 2.0 * nn.functional.mse_loss(
                shift_prediction / 0.08,
                shift_tensor / 0.08,
            )
        )
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    model.eval()
    with torch.no_grad():
        logits, center_prediction, shift_prediction = model(state_tensor)
        surface_error = (
            torch.argmax(logits, dim=1) != surface_tensor
        ).float()
        regression_error = (
            (center_prediction - center_tensor) ** 2
            + (shift_prediction - shift_tensor) ** 2
        )
        training_mse = float(
            torch.mean(surface_error + regression_error)
        )

    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state_dimension": states.shape[1],
            "model": model.state_dict(),
            "state_mean": state_mean,
            "state_std": state_std,
        },
        checkpoint,
    )
    return TrainingReport(
        samples=len(examples),
        training_mse=training_mse,
        validation_mse=None,
    )


class PolicyPredictor:
    def __init__(self, checkpoint: Path) -> None:
        bundle = torch.load(
            checkpoint,
            map_location="cpu",
            weights_only=False,
        )
        self.state_mean = np.asarray(
            bundle["state_mean"],
            dtype=np.float32,
        )
        self.state_std = np.asarray(
            bundle["state_std"],
            dtype=np.float32,
        )
        self.model = PolicyNetwork(int(bundle["state_dimension"]))
        self.model.load_state_dict(bundle["model"])
        self.model.eval()

    def predict(self, state: np.ndarray) -> Action:
        normalized = ((state - self.state_mean) / self.state_std)[None, :]
        with torch.no_grad():
            logits, center, shift = self.model(
                torch.from_numpy(normalized.astype(np.float32))
            )
        surface = "U" if int(torch.argmax(logits, dim=1)[0]) == 0 else "L"
        return Action(
            surface=surface,
            center=float(center[0]),
            shift=float(shift[0]),
        )


def sample_actions(
    rng: np.random.Generator,
    count: int,
    *,
    good_direction_probability: float = 0.70,
) -> list[Action]:
    actions: list[Action] = []
    for _ in range(count):
        surface = "U" if rng.random() < 0.5 else "L"
        center = float(rng.uniform(0.0, 1.0))
        magnitude = float(rng.uniform(0.015, 0.040))
        good_direction = rng.random() < good_direction_probability
        if surface == "U":
            shift = -magnitude if good_direction else magnitude
        else:
            shift = magnitude if good_direction else -magnitude
        actions.append(Action(surface, center, shift))
    return actions


def discrete_actions(
    *,
    center_count: int = 21,
    magnitudes: tuple[float, ...] = (0.005, 0.01, 0.02, 0.03, 0.04),
    include_reverse: bool = True,
) -> list[Action]:
    """Build a transparent full-chord action bank for this teaching case."""

    if center_count < 2:
        raise ValueError("center_count must be at least two")
    if not magnitudes or any(
        magnitude <= 0.0 or magnitude > 0.04
        for magnitude in magnitudes
    ):
        raise ValueError("magnitudes must lie in (0, 0.04]")

    actions: list[Action] = []
    for center in np.linspace(0.0, 1.0, center_count):
        for magnitude in magnitudes:
            actions.append(Action("U", float(center), -magnitude))
            actions.append(Action("L", float(center), magnitude))
            if include_reverse:
                actions.append(Action("U", float(center), magnitude))
                actions.append(Action("L", float(center), -magnitude))
    return actions


def discrete_thickness_reducing_actions(
    *,
    center_count: int = 21,
    magnitudes: tuple[float, ...] = (0.005, 0.01, 0.02, 0.03, 0.04),
) -> list[Action]:
    """Compatibility wrapper for the earlier one-direction controller."""

    return discrete_actions(
        center_count=center_count,
        magnitudes=magnitudes,
        include_reverse=False,
    )
