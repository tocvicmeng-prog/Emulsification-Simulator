"""Single-component Langmuir isotherm.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase C.

The Langmuir isotherm describes reversible monolayer adsorption:

    q* = q_max * K_L * C / (1 + K_L * C)

where:
    q*    : equilibrium loading [mol/m^3 solid]
    q_max : maximum binding capacity [mol/m^3 solid]
    K_L   : Langmuir equilibrium constant [m^3/mol]
    C     : mobile-phase concentration [mol/m^3]
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class LangmuirIsotherm:
    """Single-component Langmuir isotherm parameters.

    Attributes:
        q_max: Maximum binding capacity [mol/m^3 solid].
        K_L: Langmuir equilibrium constant [m^3/mol].
    """

    q_max: float = 100.0   # [mol/m^3 solid]
    K_L: float = 1e3       # [m^3/mol]

    @property
    def K_d(self) -> float:
        """Dissociation constant K_d = 1/K_L [mol/m^3]."""
        return 1.0 / self.K_L

    def equilibrium_loading(self, C: np.ndarray | float) -> np.ndarray | float:
        """Compute equilibrium loading q*(C).

        q* = q_max * K_L * C / (1 + K_L * C)

        Args:
            C: Mobile-phase concentration [mol/m^3].
               Accepts scalar or array.

        Returns:
            Equilibrium loading [mol/m^3 solid], same shape as C.
        """
        C = np.asarray(C, dtype=float)
        # Protect against negative concentrations
        C_safe = np.maximum(C, 0.0)
        return self.q_max * self.K_L * C_safe / (1.0 + self.K_L * C_safe)

    def jacobian(self, C: np.ndarray | float) -> np.ndarray | float:
        """Analytical derivative dq*/dC.

        dq*/dC = q_max * K_L / (1 + K_L * C)^2

        Args:
            C: Mobile-phase concentration [mol/m^3].

        Returns:
            dq*/dC [m^3_solid/m^3], same shape as C.
        """
        C = np.asarray(C, dtype=float)
        C_safe = np.maximum(C, 0.0)
        denom = (1.0 + self.K_L * C_safe) ** 2
        return self.q_max * self.K_L / denom

    def validate(self) -> list[str]:
        """Check isotherm parameter validity.

        Returns:
            List of violation messages (empty = OK).
        """
        errors: list[str] = []
        if self.q_max <= 0:
            errors.append(f"q_max must be positive, got {self.q_max}")
        if self.K_L <= 0:
            errors.append(f"K_L must be positive, got {self.K_L}")
        return errors
