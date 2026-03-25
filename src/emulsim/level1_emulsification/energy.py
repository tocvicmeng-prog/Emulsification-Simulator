"""Energy dissipation rate models for rotor-stator mixers."""

from __future__ import annotations

import numpy as np

from ..datatypes import MixerGeometry


def power_draw(mixer: MixerGeometry, rpm: float, rho_emulsion: float) -> float:
    """Power draw of rotor-stator mixer [W].

    P = N_p · ρ · N³ · D⁵

    Parameters
    ----------
    mixer : MixerGeometry
        Mixer geometry specification.
    rpm : float
        Rotational speed [rev/min].
    rho_emulsion : float
        Emulsion density [kg/m³].

    Returns
    -------
    float
        Power draw [W].
    """
    N = rpm / 60.0  # rev/s
    return mixer.power_number * rho_emulsion * N**3 * mixer.rotor_diameter**5


def average_dissipation(mixer: MixerGeometry, rpm: float, rho_emulsion: float) -> float:
    """Average turbulent energy dissipation rate [m²/s³].

    ε_avg = P / (ρ · V_tank)
    """
    P = power_draw(mixer, rpm, rho_emulsion)
    return P / (rho_emulsion * mixer.tank_volume)


def max_dissipation(mixer: MixerGeometry, rpm: float, rho_emulsion: float) -> float:
    """Maximum energy dissipation rate in rotor-stator gap [m²/s³].

    ε_max = k_ε · ε_avg
    where k_ε is the dissipation ratio (typically 10-100 for rotor-stator).
    """
    return mixer.dissipation_ratio * average_dissipation(mixer, rpm, rho_emulsion)


def gap_shear_rate(mixer: MixerGeometry, rpm: float) -> float:
    """Characteristic shear rate in the rotor-stator gap [1/s].

    γ̇ = π·N·D / δ
    where δ is the gap width.
    """
    N = rpm / 60.0
    return np.pi * N * mixer.rotor_diameter / mixer.gap_width


def kolmogorov_length_scale(epsilon: float, nu_c: float) -> float:
    """Kolmogorov microscale [m].

    η_K = (ν³/ε)^(1/4)

    Parameters
    ----------
    epsilon : float
        Energy dissipation rate [m²/s³].
    nu_c : float
        Kinematic viscosity of continuous phase [m²/s].
    """
    if epsilon <= 0:
        return np.inf
    return (nu_c**3 / epsilon)**0.25


def emulsion_density(rho_c: float, rho_d: float, phi_d: float) -> float:
    """Effective emulsion density [kg/m³].

    ρ_eff = ρ_c·(1 - φ_d) + ρ_d·φ_d
    """
    return rho_c * (1.0 - phi_d) + rho_d * phi_d
