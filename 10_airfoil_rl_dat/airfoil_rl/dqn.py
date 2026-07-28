"""Implement deterministic Double DQN, replay storage, and checkpoint I/O.

This import-only module has no command-line interface. Checkpoints are read
from or written to caller-supplied paths, normally
`output/checkpoints/dqn_best.pt`. NumPy and PyTorch are required.
"""

from __future__ import annotations

import random
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch import nn


def configure_reproducibility(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)


class QNetwork(nn.Module):
    def __init__(self, state_dimension: int, action_count: int) -> None:
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(state_dimension, 192),
            nn.SiLU(),
            nn.Linear(192, 192),
            nn.SiLU(),
            nn.Linear(192, action_count),
        )

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        return self.network(state)


@dataclass(frozen=True)
class Batch:
    states: np.ndarray
    actions: np.ndarray
    rewards: np.ndarray
    next_states: np.ndarray
    dones: np.ndarray
    next_masks: np.ndarray


class ReplayBuffer:
    def __init__(self, capacity: int, seed: int) -> None:
        self.storage: deque[
            tuple[np.ndarray, int, float, np.ndarray, bool, np.ndarray]
        ] = deque(maxlen=capacity)
        self.rng = np.random.default_rng(seed)

    def __len__(self) -> int:
        return len(self.storage)

    def append(
        self,
        state: np.ndarray,
        action: int,
        reward: float,
        next_state: np.ndarray,
        done: bool,
        next_mask: np.ndarray,
    ) -> None:
        self.storage.append(
            (
                state.astype(np.float32, copy=True),
                int(action),
                float(reward),
                next_state.astype(np.float32, copy=True),
                bool(done),
                next_mask.astype(bool, copy=True),
            )
        )

    def sample(self, batch_size: int) -> Batch:
        indices = self.rng.choice(
            len(self.storage),
            size=batch_size,
            replace=False,
        )
        values = [self.storage[int(index)] for index in indices]
        return Batch(
            states=np.stack([value[0] for value in values]),
            actions=np.asarray([value[1] for value in values], dtype=np.int64),
            rewards=np.asarray([value[2] for value in values], dtype=np.float32),
            next_states=np.stack([value[3] for value in values]),
            dones=np.asarray([value[4] for value in values], dtype=np.float32),
            next_masks=np.stack([value[5] for value in values]),
        )


class DoubleDQNAgent:
    def __init__(
        self,
        state_dimension: int,
        action_count: int,
        *,
        seed: int,
        learning_rate: float = 5.0e-4,
        gamma: float = 0.98,
        target_update_interval: int = 300,
    ) -> None:
        configure_reproducibility(seed)
        self.state_dimension = state_dimension
        self.action_count = action_count
        self.gamma = gamma
        self.target_update_interval = target_update_interval
        self.online = QNetwork(state_dimension, action_count)
        self.target = QNetwork(state_dimension, action_count)
        self.target.load_state_dict(self.online.state_dict())
        self.target.eval()
        self.optimizer = torch.optim.AdamW(
            self.online.parameters(),
            lr=learning_rate,
            weight_decay=1.0e-5,
        )
        self.update_count = 0
        self.rng = np.random.default_rng(seed + 1)

    def select_action(
        self,
        state: np.ndarray,
        epsilon: float,
        valid_mask: np.ndarray | None = None,
    ) -> int:
        valid_indices = (
            np.arange(self.action_count)
            if valid_mask is None
            else np.flatnonzero(valid_mask)
        )
        if valid_indices.size == 0:
            raise ValueError("valid action mask contains no action")
        if float(self.rng.random()) < epsilon:
            return int(self.rng.choice(valid_indices))
        with torch.no_grad():
            values = self.online(
                torch.from_numpy(state[None, :]).float()
            ).squeeze(0)
            if valid_mask is not None:
                values = values.masked_fill(
                    ~torch.from_numpy(valid_mask),
                    -torch.inf,
                )
        return int(torch.argmax(values).item())

    def greedy_action(
        self,
        state: np.ndarray,
        valid_mask: np.ndarray | None = None,
    ) -> int:
        return self.select_action(
            state,
            epsilon=0.0,
            valid_mask=valid_mask,
        )

    def optimize(self, batch: Batch) -> float:
        states = torch.from_numpy(batch.states)
        actions = torch.from_numpy(batch.actions)
        rewards = torch.from_numpy(batch.rewards)
        next_states = torch.from_numpy(batch.next_states)
        dones = torch.from_numpy(batch.dones)
        next_masks = torch.from_numpy(batch.next_masks)

        q_values = self.online(states).gather(
            1,
            actions[:, None],
        ).squeeze(1)
        with torch.no_grad():
            online_next = self.online(next_states).masked_fill(
                ~next_masks,
                -torch.inf,
            )
            next_actions = torch.argmax(online_next, dim=1)
            next_values = self.target(next_states).gather(
                1,
                next_actions[:, None],
            ).squeeze(1)
            targets = rewards + self.gamma * (1.0 - dones) * next_values
        loss = nn.functional.smooth_l1_loss(q_values, targets)
        self.optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(self.online.parameters(), 5.0)
        self.optimizer.step()
        self.update_count += 1
        if self.update_count % self.target_update_interval == 0:
            self.target.load_state_dict(self.online.state_dict())
        return float(loss.detach())

    def save(self, filename: Path, *, metadata: dict) -> None:
        filename.parent.mkdir(parents=True, exist_ok=True)
        torch.save(
            {
                "state_dimension": self.state_dimension,
                "action_count": self.action_count,
                "online": self.online.state_dict(),
                "metadata": metadata,
            },
            filename,
        )

    @classmethod
    def load(cls, filename: Path) -> "DoubleDQNAgent":
        bundle = torch.load(
            filename,
            map_location="cpu",
            weights_only=False,
        )
        metadata = bundle.get("metadata", {})
        agent = cls(
            int(bundle["state_dimension"]),
            int(bundle["action_count"]),
            seed=int(metadata.get("seed", 2026)),
        )
        agent.online.load_state_dict(bundle["online"])
        agent.target.load_state_dict(bundle["online"])
        agent.online.eval()
        agent.target.eval()
        return agent
