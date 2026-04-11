"""Catalytic packed-bed reactor simulation (Phase D).

Provides plug-flow reactor with axial dispersion, Michaelis-Menten kinetics,
effectiveness factor via generalized Thiele modulus, and first-order enzyme
deactivation.
"""

from .packed_bed import CatalyticResult, solve_packed_bed
from .kinetics import (
    michaelis_menten_rate,
    thiele_modulus,
    effectiveness_factor,
    generalized_thiele_modulus,
)
from .deactivation import (
    first_order_deactivation,
    half_life,
    arrhenius_deactivation_rate,
)

__all__ = [
    "CatalyticResult",
    "solve_packed_bed",
    "michaelis_menten_rate",
    "thiele_modulus",
    "effectiveness_factor",
    "generalized_thiele_modulus",
    "first_order_deactivation",
    "half_life",
    "arrhenius_deactivation_rate",
]
