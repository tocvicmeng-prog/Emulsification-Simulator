"""Level 2: Gelation and pore formation via Cahn-Hilliard phase-field model."""

from .solver import CahnHilliardSolver, CahnHilliard2DSolver, solve_gelation
from .ternary_solver import TernaryCahnHilliard2DSolver

__all__ = [
    "CahnHilliardSolver",
    "CahnHilliard2DSolver",
    "TernaryCahnHilliard2DSolver",
    "solve_gelation",
]
