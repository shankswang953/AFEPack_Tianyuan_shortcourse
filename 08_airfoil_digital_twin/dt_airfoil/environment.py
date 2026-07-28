"""Manage real AFEPack mesh transitions, validation, rollback, and remeshing.

This is an internal library module, normally imported by the parent-directory
entry scripts. `AirfoilMeshEnvironment` reads `data/*.dat`, builds the
`backend/` programs, and maintains generated state under `output/current/`,
`output/reference/`, and `output/history/`. Set EASYMESH_BIN and
EASYMESH2MESH_BIN when the executables are not under `$HOME/bin`.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .geometry import (
    Action,
    Airfoil,
    apply_action,
    data_shape_mse,
    make_circle_from_x_grid,
    read_airfoil,
    state_vector,
    write_airfoil,
)


@dataclass
class StepResult:
    requested_action: Action
    effective_action: Action
    state_before: np.ndarray
    state_after_trial: np.ndarray
    loss_before: float
    loss_after_trial: float
    reward: float
    valid_mesh: bool
    accepted: bool
    minimum_quality: float | None
    remeshed: bool = False
    message: str = ""


class AirfoilMeshEnvironment:
    """One real transition is boundary motion plus 200 smoothing sweeps."""

    boundary_points = 96
    smoothing_iterations = 200
    relaxation = 0.45
    gaussian_width = 0.12
    shift_max = 0.08
    minimum_mesh_quality = 0.40

    def __init__(
        self,
        root: Path | None = None,
        *,
        easymesh: Path | None = None,
        easymesh2mesh: Path | None = None,
    ) -> None:
        self.root = (
            Path(root).resolve()
            if root is not None
            else Path(__file__).resolve().parents[1]
        )
        self.backend = self.root / "backend"
        self.data_dir = self.root / "data"
        self.output_dir = self.root / "output"
        self.reference_dir = self.output_dir / "reference"
        self.current_dir = self.output_dir / "current"
        self.trial_dir = self.output_dir / "trial"
        self.remesh_trial_dir = self.output_dir / "remesh_trial"
        self.history_dir = self.output_dir / "history"

        self.target_dat = self.data_dir / "target_naca0012.dat"
        self.initial_dat = self.data_dir / "initial_circle.dat"
        self.working_dat = self.data_dir / "working.dat"

        self.generate_geometry = (
            self.backend / "generate_airfoil_geometry"
        )
        self.move_and_smooth = self.backend / "move_and_smooth"
        self.easymesh = (
            Path(easymesh)
            if easymesh is not None
            else Path(
                os.environ.get(
                    "EASYMESH_BIN",
                    str(Path.home() / "bin" / "easymesh"),
                )
            )
        )
        self.easymesh2mesh = (
            Path(easymesh2mesh)
            if easymesh2mesh is not None
            else Path(
                os.environ.get(
                    "EASYMESH2MESH_BIN",
                    str(Path.home() / "bin" / "easymesh2mesh"),
                )
            )
        )

        self.target_airfoil: Airfoil | None = None
        self.current_loss: float | None = None
        self.accepted_steps = 0

    @property
    def initial_reference(self) -> Path:
        return self.reference_dir / "initial"

    @property
    def target_reference(self) -> Path:
        return self.reference_dir / "target"

    @property
    def current_mesh(self) -> Path:
        return self.current_dir / "mesh_current.mesh"

    @property
    def current_boundary(self) -> Path:
        return self.current_dir / "boundary_current.dat"

    def _run(
        self,
        command: list[str],
        *,
        cwd: Path | None = None,
        allow_failure: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            command,
            cwd=cwd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if result.returncode != 0 and not allow_failure:
            raise RuntimeError(
                f"command failed ({result.returncode}): "
                f"{' '.join(command)}\n{result.stdout}"
            )
        return result

    def build_backend(self) -> None:
        self._run(["make", "-C", str(self.backend)])
        for executable in (self.generate_geometry, self.move_and_smooth):
            if not executable.exists():
                raise RuntimeError(f"backend executable is missing: {executable}")

    def _generate(
        self,
        initial_dat: Path,
        moved_dat: Path,
        destination: Path,
    ) -> None:
        destination.mkdir(parents=True, exist_ok=True)
        self._run(
            [
                str(self.generate_geometry),
                str(initial_dat),
                str(destination),
                str(self.boundary_points),
                "0.35",
                str(self.gaussian_width),
                "0.0",
                "0.0",
                str(moved_dat),
            ]
        )

    def _easy_mesh(self, destination: Path) -> None:
        result = self._run(
            [str(self.easymesh), "airfoil.d"],
            cwd=destination,
            allow_failure=True,
        )
        expected = [
            destination / "airfoil.n",
            destination / "airfoil.e",
            destination / "airfoil.s",
        ]
        if result.returncode != 0 and not all(path.exists() for path in expected):
            raise RuntimeError(
                "EasyMesh failed and did not create its mesh files:\n"
                + result.stdout
            )
        self._run(
            [str(self.easymesh2mesh), "airfoil", "airfoil.mesh"],
            cwd=destination,
        )

    def _smooth_transition(
        self,
        mesh: Path,
        initial_boundary: Path,
        moved_boundary: Path,
        destination: Path,
    ) -> subprocess.CompletedProcess[str]:
        return self._run(
            [
                str(self.move_and_smooth),
                str(mesh),
                str(initial_boundary),
                str(moved_boundary),
                str(destination),
                str(self.smoothing_iterations),
                str(self.relaxation),
            ],
            allow_failure=True,
        )

    def _build_reference(self, data_file: Path, destination: Path) -> None:
        if destination.exists():
            shutil.rmtree(destination)
        destination.mkdir(parents=True)
        self._generate(data_file, data_file, destination)
        self._easy_mesh(destination)
        result = self._smooth_transition(
            destination / "airfoil.mesh",
            destination / "boundary_initial.dat",
            destination / "boundary_moved.dat",
            destination,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"reference mesh is invalid for {data_file}:\n{result.stdout}"
            )

    def _remesh_candidate(
        self,
        candidate_file: Path,
    ) -> subprocess.CompletedProcess[str]:
        if self.remesh_trial_dir.exists():
            shutil.rmtree(self.remesh_trial_dir)
        self.remesh_trial_dir.mkdir(parents=True)
        self._generate(
            candidate_file,
            candidate_file,
            self.remesh_trial_dir,
        )
        self._easy_mesh(self.remesh_trial_dir)
        return self._smooth_transition(
            self.remesh_trial_dir / "airfoil.mesh",
            self.remesh_trial_dir / "boundary_initial.dat",
            self.remesh_trial_dir / "boundary_moved.dat",
            self.remesh_trial_dir,
        )

    def prepare(self, *, force: bool = False) -> None:
        if not self.target_dat.exists():
            raise FileNotFoundError(self.target_dat)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        if force or not self.initial_dat.exists():
            target = read_airfoil(self.target_dat)
            write_airfoil(
                self.initial_dat,
                make_circle_from_x_grid(target),
            )

        self.build_backend()
        required = (
            "mesh_smoothed.mesh",
            "mesh_smoothed_nodes.csv",
            "mesh_smoothed_elements.csv",
            "boundary_moved.dat",
        )
        for data_file, destination in (
            (self.initial_dat, self.initial_reference),
            (self.target_dat, self.target_reference),
        ):
            if force or not all(
                (destination / name).exists() for name in required
            ):
                self._build_reference(data_file, destination)

        self.target_airfoil = read_airfoil(self.target_dat)

    def reset(self) -> np.ndarray:
        if self.target_airfoil is None:
            self.prepare()
        if self.current_dir.exists():
            shutil.rmtree(self.current_dir)
        self.current_dir.mkdir(parents=True)
        self.history_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(
            self.initial_reference / "mesh_smoothed.mesh",
            self.current_mesh,
        )
        shutil.copy2(
            self.initial_reference / "boundary_moved.dat",
            self.current_boundary,
        )
        shutil.copy2(self.initial_dat, self.working_dat)
        current_airfoil = read_airfoil(self.working_dat)
        assert self.target_airfoil is not None
        self.current_loss = data_shape_mse(
            current_airfoil,
            self.target_airfoil,
        )
        # Accepted shapes form one global history across episode resets and
        # separate program invocations.  Never overwrite an earlier rollout.
        if self.accepted_steps == 0:
            existing_steps: list[int] = []
            for path in self.history_dir.glob("accepted_*.dat"):
                try:
                    existing_steps.append(int(path.stem.split("_")[-1]))
                except ValueError:
                    continue
            self.accepted_steps = max(existing_steps, default=0)
        return state_vector(current_airfoil)

    @staticmethod
    def _minimum_quality(directory: Path) -> float | None:
        summary = directory / "quality_summary.csv"
        if not summary.exists():
            return None
        with summary.open(newline="") as stream:
            for row in csv.DictReader(stream):
                if row["stage"] == "smoothed":
                    return float(row["min_shape_quality"])
        return None

    def current_state(self) -> np.ndarray:
        if not self.working_dat.exists():
            return self.reset()
        return state_vector(read_airfoil(self.working_dat))

    def step(
        self,
        action: Action,
        *,
        accept_only_if_improves: bool = False,
        minimum_improvement: float = 0.0,
        minimum_mesh_quality: float | None = None,
        remesh_if_needed: bool = False,
        invalid_penalty: float = 0.05,
    ) -> StepResult:
        if self.current_loss is None:
            self.reset()
        assert self.current_loss is not None
        assert self.target_airfoil is not None

        current_airfoil = read_airfoil(self.working_dat)
        state_before = state_vector(current_airfoil)
        candidate, effective = apply_action(
            current_airfoil,
            action,
            width=self.gaussian_width,
            shift_max=self.shift_max,
            thickness_reference=self.target_airfoil,
        )
        state_after = state_vector(candidate)

        if self.trial_dir.exists():
            shutil.rmtree(self.trial_dir)
        self.trial_dir.mkdir(parents=True)
        candidate_file = self.trial_dir / "candidate.dat"
        write_airfoil(candidate_file, candidate)

        try:
            self._generate(
                self.working_dat,
                candidate_file,
                self.trial_dir,
            )
            result = self._smooth_transition(
                self.current_mesh,
                self.current_boundary,
                self.trial_dir / "boundary_moved.dat",
                self.trial_dir,
            )
        except Exception as error:
            return StepResult(
                requested_action=action,
                effective_action=effective,
                state_before=state_before,
                state_after_trial=state_after,
                loss_before=self.current_loss,
                loss_after_trial=self.current_loss + invalid_penalty,
                reward=-invalid_penalty,
                valid_mesh=False,
                accepted=False,
                minimum_quality=None,
                message=str(error),
            )

        quality_floor = (
            self.minimum_mesh_quality
            if minimum_mesh_quality is None
            else minimum_mesh_quality
        )
        if quality_floor < 0.0:
            raise ValueError("minimum_mesh_quality must be nonnegative")

        evaluation_dir = self.trial_dir
        minimum_quality = self._minimum_quality(evaluation_dir)
        fixed_topology_failed = (
            result.returncode != 0
            or (
                minimum_quality is not None
                and minimum_quality < quality_floor
            )
        )
        remeshed = False
        transition_message = result.stdout

        if fixed_topology_failed and remesh_if_needed:
            try:
                remesh_result = self._remesh_candidate(candidate_file)
                remesh_quality = self._minimum_quality(
                    self.remesh_trial_dir
                )
                remesh_is_valid = (
                    remesh_result.returncode == 0
                    and (
                        remesh_quality is None
                        or remesh_quality >= quality_floor
                    )
                )
                if remesh_is_valid:
                    evaluation_dir = self.remesh_trial_dir
                    minimum_quality = remesh_quality
                    remeshed = True
                    transition_message = (
                        "fixed-topology motion failed the mesh-quality "
                        "check; EasyMesh rebuilt the candidate topology"
                    )
                else:
                    minimum_quality = remesh_quality
                    transition_message = (
                        "fixed-topology motion failed and the EasyMesh "
                        "candidate also failed the mesh-quality check\n"
                        + remesh_result.stdout
                    )
            except Exception as error:
                transition_message = (
                    "fixed-topology motion failed and EasyMesh remeshing "
                    f"failed: {error}"
                )

        if fixed_topology_failed and not remeshed:
            if (
                result.returncode == 0
                and minimum_quality is not None
                and minimum_quality < quality_floor
                and not remesh_if_needed
            ):
                transition_message = (
                    f"minimum mesh quality {minimum_quality:.6g} is below "
                    f"the safety floor {quality_floor:.6g}; "
                    "persistent mesh was not updated"
                )
            return StepResult(
                requested_action=action,
                effective_action=effective,
                state_before=state_before,
                state_after_trial=state_after,
                loss_before=self.current_loss,
                loss_after_trial=self.current_loss + invalid_penalty,
                reward=-invalid_penalty,
                valid_mesh=False,
                accepted=False,
                minimum_quality=minimum_quality,
                message=transition_message,
            )

        loss_after = data_shape_mse(candidate, self.target_airfoil)
        reward = self.current_loss - loss_after
        accepted = (
            not accept_only_if_improves
            or reward > minimum_improvement
        )
        loss_before = self.current_loss
        if accepted:
            shutil.copy2(
                evaluation_dir / "mesh_smoothed.mesh",
                self.current_mesh,
            )
            shutil.copy2(
                evaluation_dir / "boundary_moved.dat",
                self.current_boundary,
            )
            shutil.copy2(candidate_file, self.working_dat)
            self.current_loss = loss_after
            self.accepted_steps += 1
            shutil.copy2(
                self.working_dat,
                self.history_dir
                / f"accepted_{self.accepted_steps:04d}.dat",
            )
            if remeshed:
                remesh_prefix = (
                    self.history_dir
                    / f"remeshed_{self.accepted_steps:04d}"
                )
                shutil.copy2(
                    evaluation_dir / "mesh_smoothed.mesh",
                    remesh_prefix.with_suffix(".mesh"),
                )
                shutil.copy2(
                    evaluation_dir / "mesh_smoothed_nodes.csv",
                    self.history_dir
                    / f"{remesh_prefix.name}_nodes.csv",
                )
                shutil.copy2(
                    evaluation_dir / "mesh_smoothed_elements.csv",
                    self.history_dir
                    / f"{remesh_prefix.name}_elements.csv",
                )

        return StepResult(
            requested_action=action,
            effective_action=effective,
            state_before=state_before,
            state_after_trial=state_after,
            loss_before=loss_before,
            loss_after_trial=loss_after,
            reward=reward,
            valid_mesh=True,
            accepted=accepted,
            minimum_quality=minimum_quality,
            remeshed=remeshed and accepted,
            message=transition_message,
        )
