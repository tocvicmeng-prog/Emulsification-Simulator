"""pH-dependent Protein A affinity isotherm.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

Protein A chromatography is used for IgG purification.  The IgG-Protein A
interaction is strongly pH-dependent:
  - At neutral pH (7.0-8.0): high affinity, tight binding (loading step).
  - Below pH ~3.5:           near-zero affinity, protein elutes (elution step).

The pH dependence is modelled as a sigmoid transition of the effective
affinity constant K_a(pH):

    K_a(pH) = K_a_max / (1 + exp(-steepness * (pH - pH_transition)))

So at high pH, K_a -> K_a_max; at low pH, K_a -> 0.

The equilibrium loading follows single-component Langmuir with K_a(pH):

    q* = q_max * K_a(pH) * C / (1 + K_a(pH) * C)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass
class ProteinAIsotherm:
    """pH-dependent Protein A affinity isotherm for IgG purification.

    Attributes:
        q_max: Maximum binding capacity [mol/m^3 solid].
        K_a_max: Affinity constant at optimal (neutral) pH [m^3/mol].
        pH_transition: pH at which K_a = K_a_max / 2 (inflection point).
        steepness: Sigmoid steepness parameter [pH^-1].
            Higher values = sharper pH transition.
    """

    q_max: float = 60.0          # [mol/m^3 solid]  (typical Protein A resin)
    K_a_max: float = 1e5         # [m^3/mol]         (very high affinity at neutral pH)
    pH_transition: float = 3.5   # pH at half-maximal K_a
    steepness: float = 5.0       # sigmoid steepness [1/pH unit]

    @property
    def gradient_sensitive(self) -> bool:
        """Protein A equilibrium depends on pH (gradient-sensitive)."""
        return True

    @property
    def gradient_field(self) -> str:
        """ProcessState field that modulates this isotherm."""
        return "ph"

    def K_a(self, pH: float) -> float:
        """Effective affinity constant at a given pH.

        K_a(pH) = K_a_max / (1 + exp(-steepness * (pH - pH_transition)))

        Args:
            pH: Solution pH.

        Returns:
            Effective affinity constant [m^3/mol].
        """
        # Sigmoid: approaches K_a_max for pH >> pH_transition,
        #          approaches 0     for pH << pH_transition.
        exponent = -self.steepness * (pH - self.pH_transition)
        # Clip to avoid overflow
        exponent = max(min(exponent, 500.0), -500.0)
        return self.K_a_max / (1.0 + math.exp(exponent))

    def equilibrium_loading(
        self,
        C: float | np.ndarray,
        pH: float,
    ) -> np.ndarray:
        """Compute equilibrium loading at a given pH.

        q* = q_max * K_a(pH) * C / (1 + K_a(pH) * C)

        Args:
            C: Mobile-phase IgG concentration [mol/m^3].
               Accepts scalar or array.
            pH: Solution pH (uniform along the column for isotherm evaluation).

        Returns:
            Equilibrium loading [mol/m^3 solid], same shape as C.
        """
        C = np.asarray(C, dtype=float)
        C_safe = np.maximum(C, 0.0)
        Ka = self.K_a(pH)
        return self.q_max * Ka * C_safe / (1.0 + Ka * C_safe)

    def elution_pH(self, fraction: float = 0.01) -> float:
        """Find the pH at which K_a drops to a given fraction of K_a_max.

        Useful for designing elution conditions.

        Args:
            fraction: Target K_a / K_a_max (default 0.01 = 1% residual affinity).

        Returns:
            pH at which K_a(pH) = fraction * K_a_max.
        """
        if not 0 < fraction < 1:
            raise ValueError(f"fraction must be in (0, 1), got {fraction}")
        # Solve sigmoid: fraction = 1 / (1 + exp(-s*(pH - pH_t)))
        # => 1/fraction - 1 = exp(-s*(pH - pH_t))
        # => pH = pH_t - ln(1/fraction - 1) / s
        return self.pH_transition - math.log(1.0 / fraction - 1.0) / self.steepness

    def validate(self) -> list[str]:
        """Check parameter validity.

        Returns:
            List of violation messages (empty = OK).
        """
        errors: list[str] = []
        if self.q_max <= 0:
            errors.append(f"q_max must be positive, got {self.q_max}")
        if self.K_a_max <= 0:
            errors.append(f"K_a_max must be positive, got {self.K_a_max}")
        if not (0 < self.pH_transition < 14):
            errors.append(
                f"pH_transition must be in (0, 14), got {self.pH_transition}"
            )
        if self.steepness <= 0:
            errors.append(f"steepness must be positive, got {self.steepness}")
        return errors
