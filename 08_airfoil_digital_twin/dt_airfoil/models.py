"""Define the PyTorch twin and direct-policy neural networks.

This import-only module has no command-line interface and performs no file I/O.
Model training and checkpoint output are handled by `learning.py` and the
parent-directory entry scripts.
"""

from __future__ import annotations

import torch
from torch import nn


class TwinNetwork(nn.Module):
    """Predict normalized real-mesh reward from state and action."""

    def __init__(self, state_dimension: int) -> None:
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(state_dimension + 3, 256),
            nn.SiLU(),
            nn.Linear(256, 128),
            nn.SiLU(),
            nn.Linear(128, 1),
        )

    def forward(
        self,
        state: torch.Tensor,
        action: torch.Tensor,
    ) -> torch.Tensor:
        return self.network(torch.cat((state, action), dim=-1)).squeeze(-1)


class PolicyNetwork(nn.Module):
    """Map the current ``dat`` state to one action inside its training domain."""

    def __init__(self, state_dimension: int) -> None:
        super().__init__()
        self.trunk = nn.Sequential(
            nn.Linear(state_dimension, 256),
            nn.SiLU(),
            nn.Linear(256, 128),
            nn.SiLU(),
        )
        self.surface_head = nn.Linear(128, 2)
        self.center_head = nn.Linear(128, 1)
        self.shift_head = nn.Linear(128, 1)

    def forward(
        self,
        state: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        hidden = self.trunk(state)
        surface_logits = self.surface_head(hidden)
        # The center covers the full chord.  geometry.py fixes only the two
        # endpoint coordinates while allowing their neighboring points to move.
        center = torch.sigmoid(self.center_head(hidden).squeeze(-1))
        # Large deformations accumulate through multiple safe local updates.
        shift = 0.04 * torch.tanh(
            self.shift_head(hidden).squeeze(-1)
        )
        return surface_logits, center, shift
