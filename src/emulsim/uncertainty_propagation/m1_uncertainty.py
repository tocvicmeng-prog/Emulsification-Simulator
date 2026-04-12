"""M1 uncertainty contract for Monte Carlo propagation.

v6.0-alpha: Coefficients of variation for M1 export parameters.
Two tiers: "measured" (user-supplied) and "assumed" (screening defaults).
Audit F4: default CVs are clearly labeled as assumptions.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class M1UncertaintyContract:
    """Uncertainty specification for M1ExportContract parameters.

    Each cv_* field is a coefficient of variation [0, 1].
    CV = 0 means no uncertainty (deterministic).

    Tier "measured": user has experimental replicate data.
    Tier "assumed": using literature-typical screening defaults.

    Attributes:
        cv_bead_d50: CV on median bead diameter.
        cv_porosity: CV on particle porosity.
        cv_pore_size: CV on mean pore diameter.
        cv_nh2_bulk: CV on NH2 bulk concentration.
        cv_oh_bulk: CV on OH bulk concentration.
        tier: "measured" or "assumed".
    """
    cv_bead_d50: float = 0.0
    cv_porosity: float = 0.0
    cv_pore_size: float = 0.0
    cv_nh2_bulk: float = 0.0
    cv_oh_bulk: float = 0.0
    tier: str = "assumed"

    @classmethod
    def screening_defaults(cls) -> M1UncertaintyContract:
        """Literature-typical CVs for screening (not measured).

        These are ASSUMED values for exploratory analysis.
        Do not present as measured uncertainty.
        """
        return cls(
            cv_bead_d50=0.10,   # 10% CV typical for emulsification
            cv_porosity=0.05,   # 5% CV
            cv_pore_size=0.10,  # 10% CV
            cv_nh2_bulk=0.08,   # 8% CV (DDA variability)
            cv_oh_bulk=0.05,    # 5% CV
            tier="assumed",
        )

    def is_deterministic(self) -> bool:
        """True if all CVs are zero (no uncertainty)."""
        return all(cv == 0.0 for cv in [
            self.cv_bead_d50, self.cv_porosity, self.cv_pore_size,
            self.cv_nh2_bulk, self.cv_oh_bulk,
        ])
