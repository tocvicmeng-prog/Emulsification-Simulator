"""Isotherm adapters for M3 routing from FunctionalMediaContract.

v5.9.0: Enables M3 to consume M2's binding_model_hint and select the
appropriate isotherm with correct process-state parameters.

The EquilibriumAdapter wraps multi-parameter isotherms (SMA, IMAC) into
the single-argument equilibrium_loading(C) interface expected by the
LRM transport solver.

The IrreversibleAdsorptionIsotherm models one-way binding (streptavidin-biotin)
where desorption is effectively zero.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


class EquilibriumAdapter:
    """Adapts multi-parameter isotherms to single-C interface for LRM solver.

    SMA needs salt concentration. IMAC needs imidazole concentration.
    This adapter carries process_state and dispatches accordingly.
    """

    def __init__(self, isotherm, process_state: dict | None = None):
        self._isotherm = isotherm
        self._state = process_state or {}

    def equilibrium_loading(self, C):
        """Compute equilibrium loading, passing process state to multi-param isotherms."""
        # Check if isotherm has a multi-param signature
        # ASSUMPTION: SMA and IMAC isotherms accept a second positional arg
        _cls_name = type(self._isotherm).__name__

        if _cls_name == "SMAIsotherm":
            salt = self._state.get("salt_concentration", 0.1)  # [M] default 100 mM
            return self._isotherm.equilibrium_loading(C, salt)
        elif _cls_name == "IMACCompetitionIsotherm":
            imidazole = self._state.get("imidazole", 0.0)  # [M]
            return self._isotherm.equilibrium_loading(C, imidazole)
        else:
            # Standard single-arg isotherm (Langmuir, CompetitiveLangmuir)
            return self._isotherm.equilibrium_loading(C)

    @property
    def q_max(self):
        """Pass through q_max for compatibility."""
        return getattr(self._isotherm, 'q_max', 0.0)

    @property
    def K_L(self):
        """Pass through K_L for compatibility."""
        return getattr(self._isotherm, 'K_L', 0.0)


class IrreversibleAdsorptionIsotherm:
    """One-way binding model: dq/dt = k_ads * C * (q_max - q), no desorption.

    Used for streptavidin-biotin and similar near-irreversible systems
    where Kd << 1e-12 M. Once bound, material does not elute under
    mild conditions.

    For the LRM equilibrium interface: returns q_max when C > 0
    (all sites fill to capacity at any positive concentration).
    For kinetic use: provides kinetic_rate() for explicit ODE integration.
    """

    def __init__(self, q_max: float, k_ads: float = 100.0):
        self.q_max = q_max
        self.k_ads = k_ads
        self.K_L = 1e15  # effectively infinite affinity

    def equilibrium_loading(self, C):
        """At equilibrium, all sites are occupied (irreversible)."""
        return self.q_max if C > 0 else 0.0

    def kinetic_rate(self, C, q):
        """One-way adsorption rate: k_ads * C * (q_max - q)."""
        return self.k_ads * C * max(self.q_max - q, 0.0)


def select_isotherm_from_fmc(fmc, process_state: dict | None = None):
    """Select and configure an isotherm based on FunctionalMediaContract.

    Routes binding_model_hint to the appropriate isotherm class,
    wraps in EquilibriumAdapter if needed.

    Args:
        fmc: FunctionalMediaContract from M2.
        process_state: Optional dict with salt_concentration, imidazole, etc.

    Returns:
        Isotherm instance compatible with LRM solver's equilibrium_loading(C).
    """
    from ..isotherms.langmuir import LangmuirIsotherm

    hint = getattr(fmc, 'binding_model_hint', '')
    q_max = getattr(fmc, 'estimated_q_max', 0.0)
    state = process_state or {}

    if q_max <= 0:
        logger.warning(
            "FMC q_max is 0 for hint='%s'. Using default Langmuir(q_max=100, K_L=1000).",
            hint,
        )
        return LangmuirIsotherm(q_max=100.0, K_L=1000.0)

    if hint == "charge_exchange":
        # SMA isotherm for ion exchange
        try:
            from ..isotherms.sma import SMAIsotherm
            # ASSUMPTION: default SMA parameters for generic IEX
            Lambda = state.get("Lambda", 50.0)  # characteristic charge
            sigma = state.get("sigma", 10.0)    # steric factor
            iso = SMAIsotherm(q_max=q_max, Lambda=Lambda, sigma=sigma, K_eq=1.0)
            logger.info("M3 routing: charge_exchange → SMA (Lambda=%.0f, sigma=%.0f)", Lambda, sigma)
            return EquilibriumAdapter(iso, state)
        except ImportError:
            logger.warning("SMA isotherm not available. Falling back to Langmuir.")
            return LangmuirIsotherm(q_max=q_max, K_L=1000.0)

    elif hint == "metal_chelation":
        # IMAC competitive isotherm
        try:
            from ..isotherms.imac import IMACCompetitionIsotherm
            iso = IMACCompetitionIsotherm(q_max=q_max, K_protein=1e4, K_imidazole=10.0)
            logger.info("M3 routing: metal_chelation → IMAC competitive")
            return EquilibriumAdapter(iso, state)
        except ImportError:
            logger.warning("IMAC isotherm not available. Falling back to Langmuir.")
            return LangmuirIsotherm(q_max=q_max, K_L=1e4)

    elif hint == "near_irreversible":
        # Irreversible adsorption (streptavidin-biotin)
        logger.info("M3 routing: near_irreversible → IrreversibleAdsorption (q_max=%.2f)", q_max)
        return IrreversibleAdsorptionIsotherm(q_max=q_max, k_ads=100.0)

    elif hint == "salt_promoted":
        # HIC — Langmuir placeholder
        logger.info("M3 routing: salt_promoted → Langmuir placeholder (HIC)")
        return LangmuirIsotherm(q_max=q_max, K_L=100.0)

    elif hint in ("fc_affinity", "kappa_light_chain_affinity",
                   "gst_glutathione", "mixed_mode"):
        # Generic affinity — Langmuir
        K_L = state.get("K_affinity", 1e6)
        logger.info("M3 routing: %s → Langmuir (K_L=%.0e)", hint, K_L)
        return LangmuirIsotherm(q_max=q_max, K_L=K_L)

    elif hint.startswith("lectin_"):
        # Lectin — Langmuir placeholder, requires calibration
        logger.warning("M3 routing: %s → Langmuir placeholder (requires_user_calibration)", hint)
        return LangmuirIsotherm(q_max=q_max, K_L=1e3)

    else:
        # Default: generic Langmuir
        logger.info("M3 routing: '%s' → default Langmuir", hint)
        return LangmuirIsotherm(q_max=q_max, K_L=1000.0)
