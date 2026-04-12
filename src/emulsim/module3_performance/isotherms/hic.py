"""Hydrophobic Interaction Chromatography (HIC) isotherm.

v6.0-beta: Salt-modulated Langmuir model.

    K_eff = K_0 * exp(m_salt * C_salt)
    q = q_max * K_eff * C / (1 + K_eff * C)

K_0 and m_salt MUST be user-calibrated. Default values are placeholders
marked as requires_user_calibration (audit F7).

The exponential salt dependence is the standard empirical form from
Melander & Horvath (1977) and Fausnaugh & Regnier (1986). The parameter
m_salt is the molal surface tension increment, which is protein-specific
and salt-specific. It cannot be predicted from ligand density alone.
"""

from __future__ import annotations

import math
import logging
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class HICIsotherm:
    """Salt-modulated Langmuir isotherm for HIC.

    Attributes:
        q_max: Maximum binding capacity [mol/m3 solid].
        K_0: Base affinity at zero salt [m3/mol]. User-calibrated.
        m_salt: Salt sensitivity coefficient [m3/mol]. User-calibrated.
        salt_type: Salt identity for validation.
    """
    q_max: float = 50.0
    K_0: float = 0.01       # ASSUMPTION: placeholder; requires_user_calibration
    m_salt: float = 0.005   # ASSUMPTION: placeholder; protein and salt dependent
    salt_type: str = "ammonium_sulfate"

    def equilibrium_loading(
        self,
        C: float | np.ndarray,
        salt_concentration: float,
    ) -> float | np.ndarray:
        """Compute HIC equilibrium loading.

        Higher salt -> higher K_eff -> more binding (HIC characteristic).
        At zero salt, binding is minimal (K_eff = K_0, typically small).

        Args:
            C: Protein concentration in mobile phase [mol/m3].
            salt_concentration: Salt concentration [mol/m3].

        Returns:
            Equilibrium loading [mol/m3 solid].
        """
        # Guard against overflow for extreme salt
        exponent = self.m_salt * salt_concentration
        if exponent > 500:
            K_eff = self.K_0 * 1e200  # cap
        elif exponent < -500:
            K_eff = 0.0
        else:
            K_eff = self.K_0 * math.exp(exponent)

        if isinstance(C, np.ndarray):
            C_safe = np.maximum(C, 0.0)
            return self.q_max * K_eff * C_safe / (1.0 + K_eff * C_safe)
        else:
            C_safe = max(C, 0.0)
            denom = 1.0 + K_eff * C_safe
            return self.q_max * K_eff * C_safe / denom if denom > 0 else 0.0

    def validate(self) -> list[str]:
        """Check parameter validity."""
        errors = []
        if self.q_max <= 0:
            errors.append(f"q_max must be positive, got {self.q_max}")
        if self.K_0 <= 0:
            errors.append(f"K_0 must be positive, got {self.K_0}")
        if self.m_salt <= 0:
            errors.append(f"m_salt should be positive for HIC, got {self.m_salt}")
        return errors
