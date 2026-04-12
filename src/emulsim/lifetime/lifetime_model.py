"""Empirical lifetime projection for chromatography resin.

v6.0-rc: Exponential decay model for column capacity over cycles.

    capacity(n) = initial_capacity * exp(-k_deact * n)

k_deact is resin-specific and depends on CIP conditions, feed impurity
profile, storage, sanitization, and fouling. Default values are
ILLUSTRATIVE only — calibration against real cycle data is required.

This module does NOT predict:
- CIP-specific degradation mechanisms
- Fouling-dependent capacity loss
- Feed-dependent ligand leaching rates
- Storage stability

All outputs are labeled "empirical" unless user supplies calibrated k_deact.
"""

from __future__ import annotations

import math
import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class LifetimeProjection:
    """Resin lifetime projection result.

    Attributes:
        initial_capacity: Fresh resin capacity [mol/m3 or mg/mL].
        k_deactivation: First-order deactivation constant [1/cycle].
        cycles_to_80pct: Cycles until 80% of initial capacity.
        cycles_to_50pct: Cycles until 50% of initial capacity.
        capacity_at_n: Capacity at user-specified cycle number.
        n_cycles_queried: The cycle number that was queried.
        assumption_notes: CIP, storage, feed conditions assumed.
        confidence: "empirical" (always, unless calibrated).
    """
    initial_capacity: float = 0.0
    k_deactivation: float = 0.0
    cycles_to_80pct: int = 0
    cycles_to_50pct: int = 0
    capacity_at_n: float = 0.0
    n_cycles_queried: int = 0
    assumption_notes: str = ""
    confidence: str = "empirical"


def project_lifetime(
    initial_capacity: float,
    k_deactivation: float,
    n_cycles: int = 100,
    cip_description: str = "standard NaOH CIP",
) -> LifetimeProjection:
    """Project resin capacity decay over usage cycles.

    Model: capacity(n) = initial * exp(-k * n)

    Args:
        initial_capacity: Fresh resin capacity [mol/m3 or mg/mL].
        k_deactivation: First-order deactivation [1/cycle].
            Typical: Protein A ~0.002-0.01 per cycle.
            IMAC Ni: ~0.005-0.02 per cycle (metal leaching).
        n_cycles: Number of cycles to project.
        cip_description: Description of assumed CIP conditions.

    Returns:
        LifetimeProjection with cycle milestones.
    """
    if k_deactivation <= 0:
        # No deactivation: infinite lifetime
        return LifetimeProjection(
            initial_capacity=initial_capacity,
            k_deactivation=0.0,
            cycles_to_80pct=999999,
            cycles_to_50pct=999999,
            capacity_at_n=initial_capacity,
            n_cycles_queried=n_cycles,
            assumption_notes=f"k_deact=0 (no deactivation assumed). CIP: {cip_description}",
            confidence="empirical",
        )

    # Cycles to X% capacity: n = -ln(X) / k
    cycles_80 = int(math.ceil(-math.log(0.80) / k_deactivation))
    cycles_50 = int(math.ceil(-math.log(0.50) / k_deactivation))

    # Capacity at queried cycle
    cap_at_n = initial_capacity * math.exp(-k_deactivation * n_cycles)

    notes = (
        f"Exponential decay: k={k_deactivation:.4f}/cycle. "
        f"CIP: {cip_description}. "
        f"EMPIRICAL — calibrate with real cycle data."
    )

    logger.info(
        "Lifetime projection: initial=%.2f, k=%.4f, "
        "80%%@%d cycles, 50%%@%d cycles, cap@%d=%.2f",
        initial_capacity, k_deactivation,
        cycles_80, cycles_50, n_cycles, cap_at_n,
    )

    return LifetimeProjection(
        initial_capacity=initial_capacity,
        k_deactivation=k_deactivation,
        cycles_to_80pct=cycles_80,
        cycles_to_50pct=cycles_50,
        capacity_at_n=cap_at_n,
        n_cycles_queried=n_cycles,
        assumption_notes=notes,
        confidence="empirical",
    )
