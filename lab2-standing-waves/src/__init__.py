"""Public package interface for the standing waves simulation."""

from .config import get_config
from .simulation import run_simulation
from .visualization import create_2d_animation, create_3d_animation
from .solver import relax

__all__ = [
    "get_config",
    "run_simulation",
    "create_2d_animation",
    "create_3d_animation",
    "relax"
]
