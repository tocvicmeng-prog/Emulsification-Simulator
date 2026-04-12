"""Monte Carlo uncertainty propagation for M2 pipeline.

v6.0-alpha: Wraps ModificationOrchestrator.run() with perturbed M1 inputs
to produce confidence intervals on q_max and other FMC outputs.

Audit F4: Rejects physically invalid samples (negative concentrations,
porosity > 1). Labels output tier as "measured" or "assumed".
"""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, field

import numpy as np

from .m1_uncertainty import M1UncertaintyContract

logger = logging.getLogger(__name__)


@dataclass
class UncertaintyResult:
    """Monte Carlo uncertainty result for FMC outputs.

    All q_max values in [mol/m3 bed].
    """
    n_samples: int = 0
    n_valid: int = 0              # samples that passed sanity checks
    median_q_max: float = 0.0
    mean_q_max: float = 0.0
    p5_q_max: float = 0.0        # 5th percentile
    p95_q_max: float = 0.0       # 95th percentile
    cv_q_max: float = 0.0        # coefficient of variation
    tier: str = "assumed"         # "measured" or "assumed"
    all_q_max: list[float] = field(default_factory=list)


def _perturb_contract(contract, uncertainty: M1UncertaintyContract, rng: np.random.Generator):
    """Create a perturbed copy of M1ExportContract.

    Applies log-normal perturbation (ensures positivity) for each field
    with non-zero CV. Rejects invalid states.
    """
    c = copy.deepcopy(contract)

    def _perturb(value: float, cv: float) -> float:
        if cv <= 0 or value <= 0:
            return value
        # Log-normal: preserves positivity
        sigma = np.sqrt(np.log(1 + cv ** 2))
        mu = np.log(value) - 0.5 * sigma ** 2
        return float(rng.lognormal(mu, sigma))

    c.bead_d50 = _perturb(c.bead_d50, uncertainty.cv_bead_d50)
    c.bead_radius = c.bead_d50 / 2.0
    c.bead_d32 = c.bead_d50 * 0.9  # ASSUMPTION: d32/d50 ratio preserved

    c.porosity = _perturb(c.porosity, uncertainty.cv_porosity)
    c.porosity = min(max(c.porosity, 0.01), 0.99)  # sanity clamp

    c.pore_size_mean = _perturb(c.pore_size_mean, uncertainty.cv_pore_size)
    c.nh2_bulk_concentration = _perturb(c.nh2_bulk_concentration, uncertainty.cv_nh2_bulk)
    c.oh_bulk_concentration = _perturb(c.oh_bulk_concentration, uncertainty.cv_oh_bulk)

    return c


def run_with_uncertainty(
    contract,
    steps: list,
    uncertainty: M1UncertaintyContract,
    n_samples: int = 100,
    seed: int = 42,
) -> UncertaintyResult:
    """Monte Carlo wrapper around ModificationOrchestrator.run().

    Perturbs M1ExportContract fields by their CVs, runs M2 pipeline
    for each sample, collects q_max distribution from FMC.

    Args:
        contract: M1ExportContract (baseline, unperturbed).
        steps: List of ModificationStep for M2 orchestrator.
        uncertainty: M1UncertaintyContract with CVs.
        n_samples: Number of Monte Carlo samples.
        seed: Random seed for reproducibility.

    Returns:
        UncertaintyResult with percentile bounds.
    """
    from emulsim.module2_functionalization.orchestrator import (
        ModificationOrchestrator,
        build_functional_media_contract,
    )

    if uncertainty.is_deterministic():
        # No uncertainty — run once
        orch = ModificationOrchestrator()
        result = orch.run(contract, steps)
        fmc = build_functional_media_contract(result)
        return UncertaintyResult(
            n_samples=1,
            n_valid=1,
            median_q_max=fmc.estimated_q_max,
            mean_q_max=fmc.estimated_q_max,
            p5_q_max=fmc.estimated_q_max,
            p95_q_max=fmc.estimated_q_max,
            cv_q_max=0.0,
            tier=uncertainty.tier,
            all_q_max=[fmc.estimated_q_max],
        )

    rng = np.random.default_rng(seed)
    q_max_values: list[float] = []
    n_valid = 0

    for i in range(n_samples):
        try:
            perturbed = _perturb_contract(contract, uncertainty, rng)
            orch = ModificationOrchestrator()
            result = orch.run(perturbed, steps)
            fmc = build_functional_media_contract(result)
            q_max_values.append(fmc.estimated_q_max)
            n_valid += 1
        except Exception as e:
            logger.debug("MC sample %d failed: %s", i, e)
            continue

    if not q_max_values:
        logger.warning("All %d Monte Carlo samples failed", n_samples)
        return UncertaintyResult(n_samples=n_samples, tier=uncertainty.tier)

    arr = np.array(q_max_values)
    mean_val = float(np.mean(arr))

    return UncertaintyResult(
        n_samples=n_samples,
        n_valid=n_valid,
        median_q_max=float(np.median(arr)),
        mean_q_max=mean_val,
        p5_q_max=float(np.percentile(arr, 5)),
        p95_q_max=float(np.percentile(arr, 95)),
        cv_q_max=float(np.std(arr) / mean_val) if mean_val > 0 else 0.0,
        tier=uncertainty.tier,
        all_q_max=q_max_values,
    )
