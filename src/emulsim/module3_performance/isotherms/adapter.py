"""Isotherm adapters for M3 routing from FunctionalMediaContract.

v5.9.0: Enables M3 to consume M2's binding_model_hint and select the
appropriate isotherm with correct process-state parameters.

v5.9.5 hardening: Fixed SMA constructor (H1), array-safe irreversible (H2).
v6.0-beta: ProcessState dataclass, HIC isotherm routing.
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

    def __init__(self, isotherm, process_state=None):
        self._isotherm = isotherm
        # Accept both dict and ProcessState (v6.0-beta backward compat)
        if process_state is None:
            self._state = {}
        elif hasattr(process_state, 'to_dict'):
            self._state = process_state.to_dict()
        elif isinstance(process_state, dict):
            self._state = process_state
        else:
            self._state = {}

    def equilibrium_loading(self, C):
        """Compute equilibrium loading, passing process state to multi-param isotherms."""
        _cls_name = type(self._isotherm).__name__

        if _cls_name == "SMAIsotherm":
            salt = self._state.get("salt_concentration", 100.0)
            C_arr = np.atleast_1d(np.asarray(C, dtype=float))
            return self._isotherm.equilibrium_loading(C_arr, salt)
        elif _cls_name == "IMACCompetitionIsotherm":
            imidazole = self._state.get("imidazole", 0.0)
            # IMAC returns (q_protein, q_imidazole) tuple; LRM only needs q_protein
            q_protein, _q_imidazole = self._isotherm.equilibrium_loading(C, imidazole)
            return q_protein
        elif _cls_name == "HICIsotherm":
            salt = self._state.get("salt_concentration", 0.0)
            return self._isotherm.equilibrium_loading(C, salt)
        elif _cls_name == "CompetitiveAffinityIsotherm":
            competitor = self._state.get("sugar_competitor", 0.0)
            return self._isotherm.equilibrium_loading(C, competitor)
        else:
            return self._isotherm.equilibrium_loading(C)

    def update_process_state(self, field: str, value: float) -> None:
        """Update a single process_state field (v6.0 H6 gradient support).

        Args:
            field: ProcessState field name (e.g., "salt_concentration", "imidazole").
            value: New value for the field.
        """
        self._state[field] = value

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


def select_isotherm_from_fmc(fmc, process_state=None):
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
    # v6.0-beta: Accept ProcessState or dict
    if process_state is None:
        state = {}
    elif hasattr(process_state, 'to_dict'):
        state = process_state.to_dict()
    elif isinstance(process_state, dict):
        state = process_state
    else:
        state = {}

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
        # v6.0-beta: HIC isotherm when K_0 and m_salt are calibrated
        K_0 = state.get("K_0", 0.0)
        m_salt = state.get("m_salt", 0.0)
        if K_0 > 0 and m_salt > 0:
            from ..isotherms.hic import HICIsotherm
            iso = HICIsotherm(q_max=q_max, K_0=K_0, m_salt=m_salt)
            logger.info("M3 routing: salt_promoted -> HIC (K_0=%.3e, m=%.3e)", K_0, m_salt)
            return EquilibriumAdapter(iso, state)
        else:
            logger.info("M3 routing: salt_promoted -> Langmuir (HIC K_0/m_salt not calibrated)")
            return LangmuirIsotherm(q_max=q_max, K_L=100.0)

    elif hint in ("fc_affinity", "kappa_light_chain_affinity",
                   "gst_glutathione", "mixed_mode"):
        # Generic affinity - Langmuir
        K_L = state.get("K_affinity", 1e6)
        logger.info("M3 routing: %s -> Langmuir (K_L=%.0e)", hint, K_L)
        return LangmuirIsotherm(q_max=q_max, K_L=K_L)

    elif hint.startswith("lectin_"):
        # v6.0-rc: CompetitiveAffinityIsotherm for lectin elution
        K_competitor = state.get("K_competitor", 0.0)
        if K_competitor > 0:
            from ..isotherms.competitive_affinity import CompetitiveAffinityIsotherm
            K_protein = state.get("K_affinity", 1e5)
            sugar = state.get("sugar_type", "sugar")
            iso = CompetitiveAffinityIsotherm(
                q_max=q_max, K_protein=K_protein,
                K_competitor=K_competitor, competitor_name=sugar,
            )
            logger.info("M3 routing: %s -> CompetitiveAffinity (K_comp=%.0e)", hint, K_competitor)
            return EquilibriumAdapter(iso, state)
        else:
            logger.warning("M3 routing: %s -> Langmuir (requires_user_calibration, no K_competitor)", hint)
            return LangmuirIsotherm(q_max=q_max, K_L=1e3)

    else:
        # Default: generic Langmuir
        logger.info("M3 routing: '%s' -> default Langmuir", hint)
        return LangmuirIsotherm(q_max=q_max, K_L=1000.0)


# v6.0 H6: Gradient-dependent equilibrium — IMPLEMENTED
# solve_lrm() now accepts optional gradient_program + equilibrium_adapter params.
# When provided, _build_rhs() updates adapter._state[gradient_field] at each
# time step with gradient_program.value_at_time(t), enabling time-varying
# equilibrium for SMA, HIC, IMAC, and lectin isotherms during gradient elution.
