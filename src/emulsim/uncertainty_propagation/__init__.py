"""Uncertainty propagation for EmulSim v6.0.

Provides Monte Carlo wrapper for M2 pipeline to produce confidence
intervals on FunctionalMediaContract outputs.
"""

from .m1_uncertainty import M1UncertaintyContract
from .monte_carlo import UncertaintyResult, run_with_uncertainty

__all__ = ["M1UncertaintyContract", "UncertaintyResult", "run_with_uncertainty"]
