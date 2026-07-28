"""Define the deterministic pure-`dat` airfoil RL environment.

This import-only module performs shape transitions in memory and has no direct
file output. Construct `DatAirfoilEnvironment` from caller-provided target
and initial shapes; use `../train_dqn.py` for the complete workflow.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import (
    Action,
    Airfoil,
    action_grid,
    apply_action,
    data_mse,
    error_state,
    gaussian_weights,
    interpolate_shapes,
    smooth_random_perturbation,
)


@dataclass(frozen=True)
class StepResult:
    state: np.ndarray
    reward: float
    terminated: bool
    truncated: bool
    loss: float
    improvement: float
    realized_shift: float
    action: Action


class DatAirfoilEnvironment:
    """Shape matching without a mesh or CFD solver."""

    def __init__(
        self,
        target: Airfoil,
        initial: Airfoil,
        *,
        max_steps: int = 140,
        terminal_loss: float = 2.0e-6,
        seed: int = 2026,
    ) -> None:
        self.target = target.copy()
        self.initial = initial.copy()
        self.actions = action_grid()
        self.max_steps = max_steps
        self.terminal_loss = terminal_loss
        self.rng = np.random.default_rng(seed)
        self.current = initial.copy()
        self.episode_initial_loss = data_mse(self.current, self.target)
        self.reference_loss = self.episode_initial_loss
        self.step_count = 0
        self._action_weights = np.stack(
            [
                gaussian_weights(self.target.x, action.center)
                for action in self.actions
            ]
        )

    @property
    def state_dimension(self) -> int:
        return 2 * (self.target.x.size - 2)

    @property
    def action_count(self) -> int:
        return len(self.actions)

    @property
    def loss(self) -> float:
        return data_mse(self.current, self.target)

    def state(self) -> np.ndarray:
        return error_state(self.current, self.target)

    def valid_action_mask(self) -> np.ndarray:
        upper_remaining = np.maximum(
            self.current.upper_y - self.target.upper_y,
            0.0,
        )
        lower_remaining = np.maximum(
            self.target.lower_y - self.current.lower_y,
            0.0,
        )
        mask = np.asarray(
            [
                float(
                    np.dot(
                        weights,
                        upper_remaining
                        if action.surface == "U"
                        else lower_remaining,
                    )
                )
                > 1.0e-9
                for action, weights in zip(
                    self.actions,
                    self._action_weights,
                )
            ],
            dtype=bool,
        )
        if not np.any(mask):
            mask[:] = True
        return mask

    def reset(
        self,
        *,
        evaluation: bool = False,
        alpha: float | None = None,
    ) -> np.ndarray:
        if evaluation:
            self.current = self.initial.copy(title="RL EVALUATION")
        else:
            if alpha is None:
                alpha = float(self.rng.uniform(0.30, 1.00))
            base = interpolate_shapes(self.target, self.initial, alpha)
            amplitude = 0.008 * float(alpha)
            self.current = smooth_random_perturbation(
                base,
                self.rng,
                amplitude=amplitude,
            )
        self.step_count = 0
        self.episode_initial_loss = max(self.loss, 1.0e-10)
        return self.state()

    def step(self, action_index: int) -> StepResult:
        if not 0 <= action_index < self.action_count:
            raise IndexError("action index is outside the action grid")
        action = self.actions[action_index]
        previous_loss = self.loss
        moved, realized_shift = apply_action(
            self.current,
            action,
            target_reference=self.target,
        )
        self.current = moved
        self.step_count += 1
        loss = self.loss
        improvement = previous_loss - loss

        normalized_improvement = improvement / self.reference_loss
        reward = 20.0 * normalized_improvement - 2.0e-3
        if realized_shift < 1.0e-10:
            reward -= 2.0e-2
        terminated = loss <= self.terminal_loss
        truncated = self.step_count >= self.max_steps and not terminated
        if terminated:
            reward += 1.0
        return StepResult(
            state=self.state(),
            reward=float(reward),
            terminated=terminated,
            truncated=truncated,
            loss=loss,
            improvement=improvement,
            realized_shift=realized_shift,
            action=action,
        )
