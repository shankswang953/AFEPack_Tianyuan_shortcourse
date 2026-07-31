"""Serialize transition and policy examples used by the reward model.

This import-only module reads, appends, or rewrites caller-supplied JSONL
paths, normally `output/replay.jsonl` and `output/policy_examples.jsonl`.
It has no standalone command-line interface.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .environment import StepResult
from .geometry import Action


@dataclass
class TransitionRecord:
    phase: str
    episode: int
    step: int
    state_before: list[float]
    state_after_trial: list[float]
    requested_action: dict[str, float | str]
    effective_action: dict[str, float | str]
    loss_before: float
    loss_after_trial: float
    reward: float
    valid_mesh: bool
    accepted: bool
    minimum_quality: float | None
    remeshed: bool = False

    @classmethod
    def from_result(
        cls,
        result: StepResult,
        *,
        phase: str,
        episode: int,
        step: int,
    ) -> "TransitionRecord":
        return cls(
            phase=phase,
            episode=episode,
            step=step,
            state_before=result.state_before.astype(float).tolist(),
            state_after_trial=(
                result.state_after_trial.astype(float).tolist()
            ),
            requested_action=asdict(result.requested_action),
            effective_action=asdict(result.effective_action),
            loss_before=float(result.loss_before),
            loss_after_trial=float(result.loss_after_trial),
            reward=float(result.reward),
            valid_mesh=bool(result.valid_mesh),
            accepted=bool(result.accepted),
            minimum_quality=(
                None
                if result.minimum_quality is None
                else float(result.minimum_quality)
            ),
            remeshed=bool(result.remeshed),
        )

    @property
    def state(self) -> np.ndarray:
        return np.asarray(self.state_before, dtype=np.float32)

    @property
    def action(self) -> Action:
        return Action(
            surface=str(self.effective_action["surface"]),
            center=float(self.effective_action["center"]),
            shift=float(self.effective_action["shift"]),
        )


def load_replay(filename: Path) -> list[TransitionRecord]:
    if not filename.exists():
        return []
    records: list[TransitionRecord] = []
    with filename.open() as stream:
        for line in stream:
            if line.strip():
                records.append(TransitionRecord(**json.loads(line)))
    return records


def append_replay(filename: Path, record: TransitionRecord) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("a") as stream:
        stream.write(json.dumps(asdict(record), separators=(",", ":")))
        stream.write("\n")


@dataclass
class PolicyExample:
    state: list[float]
    action: dict[str, float | str]

    @classmethod
    def create(cls, state: np.ndarray, action: Action) -> "PolicyExample":
        return cls(
            state=state.astype(float).tolist(),
            action=asdict(action),
        )


def load_policy_examples(filename: Path) -> list[PolicyExample]:
    if not filename.exists():
        return []
    with filename.open() as stream:
        return [
            PolicyExample(**json.loads(line))
            for line in stream
            if line.strip()
        ]


def append_policy_example(
    filename: Path,
    example: PolicyExample,
) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("a") as stream:
        stream.write(json.dumps(asdict(example), separators=(",", ":")))
        stream.write("\n")
