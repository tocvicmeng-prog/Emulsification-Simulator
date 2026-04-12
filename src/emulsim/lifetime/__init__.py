"""Lifetime projection module for EmulSim v6.0.

Empirical exponential decay model for column resin lifetime.
Always labeled as empirical — not predictive without calibration.
"""

from .lifetime_model import LifetimeProjection, project_lifetime

__all__ = ["LifetimeProjection", "project_lifetime"]
