"""Isotherm adapters for M3 routing from FunctionalMediaContract.

v5.9.0: Enables M3 to consume M2's binding_model_hint and select the
appropriate isotherm with correct process-state parameters.

v5.9.5 hardening: Fixed SMA constructor (H1), array-safe irreversible (H2),
added TODO for gradient process state (H6).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

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
        _cls_name = type(self._isotherm).__name__

        if _cls_name == "SMAIsotherm":
            # SMA needs salt concentration [mol/m3]
            salt = self._state.get("salt_concentration", 100.0)  # [mol/m3] default 100 mM
            # SMA expects C as array of component concentrations
            C_arr = np.atleast_1d(np.asarray(C, dtype=float))
            return self._isotherm.equilibrium_loading(C_arr, salt)
        elif _cls_name == "IMACCompetitionIsotherm":
            imidazole = self._state.get("imidazole", 0.0)  # [mol/m3]
            return self._isotherm.equilibrium_loading(C, imidazole)
        else:
            return self._isotherm.equilibrium_loading(C)

    @property
    def q_max(self):
        """Pass through q_max for compatibility."""
        return getattr(self._isotherm, 'q_max',
                       getattr(self._isotherm, 'Lambda', 0.0))

    @property
    def K_L(self):
        """Pass through K_L for compatibility."""
        return getattr(self._isotherm, 'K_L', 0.0)


class IrreversibleAdsorptionIsotherm:
    """One-way binding model: dq/dt = k_ads * C * (q_max - q), no desorption.

    Used for streptavidin-biotin and similar near-irreversible systems
    where Kd << 1e-12 M.

    v5.9.5 H2: Array-safe for LRM solver NumPy inputs.
    """

    def __init__(self, q_max: float, k_ads: float = 100.0):
        self.q_max = q_max
        self.k_ads = k_ads
        self.K_L = 1e15  # effectively infinite affinity

    def equilibrium_loading(self, C):
        """At equilibrium, all sites are occupied (irreversible).

        Array-safe: handles both scalar and NumPy array inputs.
        """
        if isinstance(C, np.ndarray):
            return np.where(C > 0, self.q_max, 0.0)
        return self.q_max if C > 0 else 0.0

    def kinetic_rate(self, C, q):
        """One-way adsorption rate: k_ads * C * (q_max - q)."""
        if isinstance(C, np.ndarray):
            return self.k_ads * C * np.maximum(self.q_max - q, 0.0)
        return self.k_ads * C * max(self.q_max - q, 0.0)


def select_isotherm_from_fmc(fmc, process_state: dict | None = None):
    """Select and configure an isotherm based on FunctionalMediaContract.

    v5.9.5 H1: Fixed SMA constructor to use correct parameter names
    (Lambda, z, sigma, K_eq as arrays).

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
        # v5.9.5 H1: SMA uses Lambda (ionic capacity), z/sigma/K_eq (arrays)
        # q_max from FMC maps to Lambda (total ionic capacity)
        try:
            from ..isotherms.sma import SMAIsotherm
            Lambda = q_max  # FMC q_max represents ionic capacity for IEX
            z = np.array([state.get("z", 5.0)])          # characteristic charge
            sigma = np.array([state.get("sigma", 50.0)])  # steric factor
            K_eq = np.array([state.get("K_eq", 1e-3)])    # equilibrium constant
            iso = SMAIsotherm(Lambda=Lambda, z=z, sigma=sigma, K_eq=K_eq)
            logger.info("M3 routing: charge_exchange -> SMA (Lambda=%.1f)", Lambda)
            return EquilibriumAdapter(iso, state)
        except (ImportError, Exception) as e:
            logger.warning("SMA isotherm failed (%s). Falling back to Langmuir.", e)
            return LangmuirIsotherm(q_max=q_max, K_L=1000.0)

    elif hint == "metal_chelation":
        # IMAC competitive isotherm
        try:
            from ..isotherms.imac import IMACCompetitionIsotherm
            iso = IMACCompetitionIsotherm(q_max=q_max, K_protein=1e4, K_imidazole=10.0)
            logger.info("M3 routing: metal_chelation -> IMAC competitive")
            return EquilibriumAdapter(iso, state)
        except (ImportError, Exception) as e:
            logger.warning("IMAC isotherm failed (%s). Falling back to Langmuir.", e)
            return LangmuirIsotherm(q_max=q_max, K_L=1e4)

    elif hint == "near_irreversible":
        # Irreversible adsorption (streptavidin-biotin)
        logger.info("M3 routing: near_irreversible -> IrreversibleAdsorption (q_max=%.2f)", q_max)
        return IrreversibleAdsorptionIsotherm(q_max=q_max, k_ads=100.0)

    elif hint == "salt_promoted":
        # HIC - Langmuir placeholder
        # TODO v6.0-beta: Replace with HICIsotherm(K_0, m_salt) when calibration framework exists
        logger.info("M3 routing: salt_promoted -> Langmuir placeholder (HIC)")
        return LangmuirIsotherm(q_max=q_max, K_L=100.0)

    elif hint in ("fc_affinity", "kappa_light_chain_affinity",
                   "gst_glutathione", "mixed_mode"):
        # Generic affinity - Langmuir
        K_L = state.get("K_affinity", 1e6)
        logger.info("M3 routing: %s -> Langmuir (K_L=%.0e)", hint, K_L)
        return LangmuirIsotherm(q_max=q_max, K_L=K_L)

    elif hint.startswith("lectin_"):
        # Lectin - Langmuir placeholder, requires calibration
        # TODO v6.0-rc: Replace with CompetitiveAffinityIsotherm when available
        logger.warning("M3 routing: %s -> Langmuir placeholder (requires_user_calibration)", hint)
        return LangmuirIsotherm(q_max=q_max, K_L=1e3)

    else:
        # Default: generic Langmuir
        logger.info("M3 routing: '%s' -> default Langmuir", hint)
        return LangmuirIsotherm(q_max=q_max, K_L=1000.0)


# TODO v6.0-beta H6: Gradient-dependent equilibrium
# Current gradient code passes a scalar GradientProgram value but does not
# propagate it into the equilibrium adapter's process_state during LRM time
# integration. This requires extending solve_lrm() to accept a time-varying
# ProcessState callback. Planned for v6.0-beta (doc 33 Section 5.1).
