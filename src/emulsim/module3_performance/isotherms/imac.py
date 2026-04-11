"""IMAC (Immobilized Metal Affinity Chromatography) isotherm with imidazole competition.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

IMAC separates His-tagged proteins via competitive binding at metal-chelate
sites.  Both the target protein and imidazole (the eluent) compete for
the same immobilized metal (Ni2+, Co2+, etc.) sites via a shared
competitive Langmuir mechanism:

    q_protein = q_max * K_protein * C_protein
                / (1 + K_protein * C_protein + K_imidazole * C_imidazole)

    q_imidazole = q_max * K_imidazole * C_imidazole
                  / (1 + K_protein * C_protein + K_imidazole * C_imidazole)

This is physically equivalent to competitive Langmuir between two species.
Imidazole is typically present at much higher concentrations than protein,
so small changes in [imidazole] can dramatically reduce protein retention.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class IMACCompetitionIsotherm:
    """Competitive Langmuir between His-tagged protein and imidazole on IMAC resin.

    Attributes:
        q_max: Total binding site capacity [mol/m^3 solid].
        K_protein: Protein-metal affinity constant [m^3/mol].
        K_imidazole: Imidazole-metal affinity constant [m^3/mol].
    """

    q_max: float = 50.0        # [mol/m^3 solid]
    K_protein: float = 1e4     # [m^3/mol]  (high affinity His-tag)
    K_imidazole: float = 10.0  # [m^3/mol]  (weak competitor)

    def equilibrium_loading(
        self,
        C_protein: float | np.ndarray,
        C_imidazole: float | np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compute equilibrium loading for protein and imidazole.

        Args:
            C_protein: Protein concentration in mobile phase [mol/m^3].
                       Accepts scalar or array.
            C_imidazole: Imidazole concentration in mobile phase [mol/m^3].
                         Accepts scalar or array (broadcast-compatible with C_protein).

        Returns:
            Tuple of (q_protein, q_imidazole), each [mol/m^3 solid],
            same shape as the broadcast of C_protein and C_imidazole.
        """
        C_p = np.asarray(C_protein, dtype=float)
        C_i = np.asarray(C_imidazole, dtype=float)

        # Clip negative concentrations
        C_p = np.maximum(C_p, 0.0)
        C_i = np.maximum(C_i, 0.0)

        denominator = 1.0 + self.K_protein * C_p + self.K_imidazole * C_i

        q_protein = self.q_max * self.K_protein * C_p / denominator
        q_imidazole = self.q_max * self.K_imidazole * C_i / denominator

        return q_protein, q_imidazole

    def protein_loading_only(
        self,
        C_protein: float | np.ndarray,
        C_imidazole: float | np.ndarray,
    ) -> np.ndarray:
        """Return only the protein equilibrium loading.

        Convenience wrapper for orchestrator use.

        Args:
            C_protein: Protein concentration [mol/m^3].
            C_imidazole: Imidazole concentration [mol/m^3].

        Returns:
            Protein equilibrium loading [mol/m^3 solid].
        """
        q_prot, _ = self.equilibrium_loading(C_protein, C_imidazole)
        return q_prot

    def imidazole_ic50(self) -> float:
        """Imidazole concentration that halves protein binding at low load.

        At the low-protein limit, q_protein ~ q_max * K_p * C_p / (1 + K_i * [Im]).
        IC50 is where (1 + K_i * [Im]) = 2, i.e., K_i * [Im] = 1.

        Returns:
            IC50 for imidazole [mol/m^3].
        """
        return 1.0 / self.K_imidazole

    def validate(self) -> list[str]:
        """Check parameter validity.

        Returns:
            List of violation messages (empty = OK).
        """
        errors: list[str] = []
        if self.q_max <= 0:
            errors.append(f"q_max must be positive, got {self.q_max}")
        if self.K_protein <= 0:
            errors.append(f"K_protein must be positive, got {self.K_protein}")
        if self.K_imidazole <= 0:
            errors.append(f"K_imidazole must be positive, got {self.K_imidazole}")
        if self.K_imidazole >= self.K_protein:
            errors.append(
                f"K_imidazole ({self.K_imidazole:.2g}) >= K_protein ({self.K_protein:.2g}); "
                "imidazole would bind more strongly than protein — check units."
            )
        return errors
