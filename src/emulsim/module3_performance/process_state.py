"""Process state for M3 chromatography simulations.

v6.0-beta: Typed dataclass replacing loose dict for process conditions.
Carries time-varying mobile-phase composition that isotherms need
for salt-dependent (HIC, IEX), competitor-dependent (IMAC, lectin),
and pH-dependent (Protein A) binding calculations.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ProcessState:
    """Time-varying process conditions for M3 equilibrium calculations.

    All concentrations in [mol/m3] for consistency with M2/M3 internals.
    1 M = 1000 mol/m3.

    Attributes:
        salt_concentration: Mobile-phase salt [mol/m3]. For IEX (SMA) and HIC.
        salt_type: Salt identity (e.g., "NaCl", "ammonium_sulfate").
        ph: Mobile-phase pH.
        imidazole: Imidazole concentration [mol/m3]. For IMAC elution.
        sugar_competitor: Sugar competitor [mol/m3]. For lectin elution.
        sugar_type: Sugar identity (e.g., "mannose", "GlcNAc").
        temperature: Temperature [K].
        conductivity: Conductivity [mS/cm] (alternative to salt conc).
        K_0: HIC base affinity (user-calibrated). 0 = not provided.
        m_salt: HIC salt sensitivity [m3/mol] (user-calibrated). 0 = not provided.
        K_affinity: Generic affinity constant [m3/mol] for Langmuir routing.
        Lambda: SMA ionic capacity override [mol/m3 solid].
        z: SMA characteristic charge override.
        sigma: SMA steric factor override.
        K_eq: SMA equilibrium constant override.
    """
    salt_concentration: float = 0.0
    salt_type: str = ""
    ph: float = 7.0
    imidazole: float = 0.0
    sugar_competitor: float = 0.0
    sugar_type: str = ""
    temperature: float = 298.15
    conductivity: float = 0.0
    # HIC calibrated parameters (user-supplied)
    K_0: float = 0.0
    m_salt: float = 0.0
    # Generic affinity
    K_affinity: float = 1e6
    # SMA overrides
    Lambda: float = 0.0
    z: float = 5.0
    sigma: float = 50.0
    K_eq: float = 1e-3

    def to_dict(self) -> dict:
        """Convert to dict for backward compatibility with dict-based adapters."""
        return {
            "salt_concentration": self.salt_concentration,
            "salt_type": self.salt_type,
            "ph": self.ph,
            "imidazole": self.imidazole,
            "sugar_competitor": self.sugar_competitor,
            "sugar_type": self.sugar_type,
            "temperature": self.temperature,
            "conductivity": self.conductivity,
            "K_0": self.K_0,
            "m_salt": self.m_salt,
            "K_affinity": self.K_affinity,
            "Lambda": self.Lambda,
            "z": self.z,
            "sigma": self.sigma,
            "K_eq": self.K_eq,
        }
