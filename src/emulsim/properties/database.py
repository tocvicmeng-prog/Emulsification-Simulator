"""Centralized material property database with temperature/concentration interpolation."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from ..datatypes import MaterialProperties


class PropertyDatabase:
    """Material property storage and retrieval with interpolation.

    Provides temperature- and concentration-dependent material properties
    for the emulsification simulation pipeline.
    """

    def __init__(self, props: Optional[MaterialProperties] = None):
        self.props = props or MaterialProperties()

    @classmethod
    def from_toml(cls, path: Path) -> PropertyDatabase:
        """Load from a TOML property file."""
        from ..config import load_properties
        return cls(load_properties(path))

    def oil_density(self, T: float) -> float:
        """Paraffin oil density [kg/m³] at temperature T [K].

        Linear approximation: rho(T) = rho_ref + alpha*(T - T_ref)
        alpha ~ -0.65 kg/(m³·K) for liquid paraffin.
        """
        T_ref = 293.15  # 20 degC reference
        alpha = -0.65   # kg/(m³·K)
        return self.props.rho_oil + alpha * (T - T_ref)

    def oil_viscosity(self, T: float) -> float:
        """Paraffin oil dynamic viscosity [Pa·s] at temperature T [K].

        Arrhenius model: mu(T) = A*exp(E_a/(R*T))
        Calibrated to mu(90 degC) = 5 mPa·s.
        """
        # Calibrate: at T_ref=363.15K, mu=0.005 Pa·s
        # At 20 degC (293.15K), mu ~ 0.03 Pa·s (typical light paraffin)
        R = 8.314
        T = max(T, 100.0)  # prevent division issues
        T_ref = 363.15
        mu_ref = self.props.mu_oil
        # From two-point fit: E_a/R ~ 2500 K
        E_a_over_R = 2500.0
        return mu_ref * np.exp(E_a_over_R * (1.0 / T - 1.0 / T_ref))

    def dispersed_viscosity(self, T: float, c_agarose: float, c_chitosan: float,
                            shear_rate: float = 0.0) -> float:
        """Dispersed phase (agarose + chitosan solution) viscosity [Pa·s].

        Uses Mark-Houwink for agarose intrinsic viscosity, Huggins equation
        for concentration dependence, and logarithmic mixing rule for blend.
        When a non-zero ``shear_rate`` is supplied the zero-shear viscosity
        is corrected for shear-thinning via the Cross model.

        Parameters
        ----------
        T : float
            Temperature [K]. Must be above T_gel for sol state.
        c_agarose : float
            Agarose concentration [kg/m³].
        c_chitosan : float
            Chitosan concentration [kg/m³].
        shear_rate : float, optional
            Shear rate [1/s].  If > 0, apply Cross model correction.
        """
        from .viscosity import solution_viscosity
        mu_0 = solution_viscosity(T, c_agarose, c_chitosan,
                                   eta_intr_chit=self.props.eta_intr_chit)
        if shear_rate > 0:
            from .viscosity import cross_model_correction
            return cross_model_correction(mu_0, shear_rate)
        return mu_0

    def interfacial_tension(self, T: float, c_span80: float) -> float:
        """Interfacial tension [N/m] for Span-80 stabilized W/O interface.

        Szyszkowski-Langmuir isotherm.

        Parameters
        ----------
        T : float
            Temperature [K].
        c_span80 : float
            Span-80 concentration in oil phase [kg/m³].
        """
        from .interfacial import interfacial_tension_span80
        return interfacial_tension_span80(T, c_span80)

    def chi_parameter(self, T: float) -> float:
        """Flory-Huggins chi parameter for agarose-water at temperature T [K].

        chi(T) = A/T + B
        """
        from .thermodynamic import chi_flory_huggins
        A, B = self.props.chi_T_coeffs
        return chi_flory_huggins(T, A, B)

    def update_for_conditions(self, T_oil: float, c_agarose: float,
                               c_chitosan: float, c_span80: float) -> MaterialProperties:
        """Return a MaterialProperties with values interpolated to actual conditions."""
        import copy
        props = copy.deepcopy(self.props)
        props.rho_oil = self.oil_density(T_oil)
        props.mu_oil = self.oil_viscosity(T_oil)
        props.mu_d = self.dispersed_viscosity(T_oil, c_agarose, c_chitosan)
        props.sigma = self.interfacial_tension(T_oil, c_span80)
        props.chi_0 = self.chi_parameter(T_oil)
        return props
