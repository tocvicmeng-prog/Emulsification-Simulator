"""Competitive affinity isotherm for lectin and similar systems.

v6.0-rc: Generalized competitive binding between target protein
and a small-molecule competitor for immobilized ligand sites.

    q = q_max * K_protein * C_protein / (1 + K_protein*C_protein + K_competitor*C_competitor)

Used for:
- Con A: protein vs mannose/glucose competition
- WGA: protein vs GlcNAc competition
- Generic sugar/metabolite affinity

Generalized from IMACCompetitionIsotherm (imidazole competition).
IMAC and lectin share the same mathematical model but different
physical contexts (metal chelation vs carbohydrate recognition).

requires_user_calibration unless K_competitor is provided.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class CompetitiveAffinityIsotherm:
    """Competitive binding isotherm for affinity chromatography.

    Attributes:
        q_max: Maximum binding capacity [mol/m3 solid].
        K_protein: Protein-ligand affinity [m3/mol].
        K_competitor: Competitor-ligand affinity [m3/mol].
        competitor_name: Identity of competitor (e.g., "mannose", "imidazole").
    """
    q_max: float = 50.0
    K_protein: float = 1e6
    K_competitor: float = 1e3
    competitor_name: str = "competitor"

    @property
    def gradient_sensitive(self) -> bool:
        """Lectin equilibrium depends on competitor (sugar) concentration."""
        return True

    @property
    def gradient_field(self) -> str:
        """ProcessState field that modulates this isotherm."""
        return "sugar_competitor"

    def equilibrium_loading(
        self,
        C_protein: float | np.ndarray,
        C_competitor: float = 0.0,
    ) -> float | np.ndarray:
        """Compute competitive equilibrium loading.

        At C_competitor=0: standard Langmuir.
        At high C_competitor: protein displaced (elution).

        Args:
            C_protein: Protein concentration [mol/m3].
            C_competitor: Competitor concentration [mol/m3].

        Returns:
            Protein loading [mol/m3 solid].
        """
        if isinstance(C_protein, np.ndarray):
            C_p = np.maximum(C_protein, 0.0)
        else:
            C_p = max(C_protein, 0.0)

        C_c = max(C_competitor, 0.0)
        denom = 1.0 + self.K_protein * C_p + self.K_competitor * C_c

        if isinstance(C_p, np.ndarray):
            return self.q_max * self.K_protein * C_p / np.maximum(denom, 1e-30)
        return self.q_max * self.K_protein * C_p / max(denom, 1e-30)

    def validate(self) -> list[str]:
        """Check parameter validity."""
        errors = []
        if self.q_max <= 0:
            errors.append(f"q_max must be positive, got {self.q_max}")
        if self.K_protein <= 0:
            errors.append(f"K_protein must be positive, got {self.K_protein}")
        if self.K_competitor < 0:
            errors.append(f"K_competitor must be non-negative, got {self.K_competitor}")
        return errors
