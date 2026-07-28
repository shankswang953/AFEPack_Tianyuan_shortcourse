"""Public imports for the pure-`dat` airfoil reinforcement-learning package.

This is an import-only module with no command-line interface or file output.
Use `../train_dqn.py` or `../run.sh` to train and evaluate the example.
"""

from .dqn import DoubleDQNAgent, ReplayBuffer
from .environment import DatAirfoilEnvironment
from .geometry import Action, Airfoil, read_airfoil, write_airfoil

__all__ = [
    "Action",
    "Airfoil",
    "DatAirfoilEnvironment",
    "DoubleDQNAgent",
    "ReplayBuffer",
    "read_airfoil",
    "write_airfoil",
]
