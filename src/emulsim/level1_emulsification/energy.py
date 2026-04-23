"""Energy dissipation rate models for emulsification mixers.

Supports:
  - Legacy rotor-stator homogeniser (MixerGeometry)
  - Pitched-blade impeller / Stirrer A (StirrerGeometry)
  - Small rotor-stator / Stirrer B (StirrerGeometry)

References
----------
- Metzner & Otto (1957), AIChE J 3(1):3-10
- Utomo et al. (2009), Chem. Eng. Sci. 64(21):4426-4439
- Padron (2005), PhD Thesis, University of Maryland
"""

from __future__ import annotations

from typing import Union

import numpy as np

from ..datatypes import MixerGeometry, StirrerGeometry

# Type alias for any geometry object accepted by the generic functions.
Geometry = Union[MixerGeometry, StirrerGeometry]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_diameter(geometry: Geometry) -> float:
    """Return the active impeller/rotor diameter [m] from any geometry."""
    # StirrerGeometry uses `impeller_diameter`; MixerGeometry uses `rotor_diameter`.
    d = getattr(geometry, "impeller_diameter",
                getattr(geometry, "rotor_diameter", None))
    if d is None:
        raise AttributeError(
            f"Geometry {type(geometry).__name__} has neither impeller_diameter "
            f"nor rotor_diameter"
        )
    return float(d)


def _get_tank_volume(geometry: Geometry) -> float:
    """Return the tank/working volume [m3] from any geometry."""
    v = getattr(geometry, "working_volume",
                getattr(geometry, "tank_volume", None))
    if v is None:
        raise AttributeError(
            f"Geometry {type(geometry).__name__} has neither working_volume "
            f"nor tank_volume"
        )
    return float(v)


# ---------------------------------------------------------------------------
# Power draw & dissipation (universal)
# ---------------------------------------------------------------------------

def power_draw(geometry: Geometry, rpm: float, rho_emulsion: float) -> float:
    """Power draw of an impeller or rotor-stator [W].

    P = Np * rho * N^3 * D^5

    Parameters
    ----------
    geometry : MixerGeometry | StirrerGeometry
        Equipment geometry.  Must expose ``power_number`` and either
        ``rotor_diameter`` (legacy) or ``impeller_diameter`` (new).
    rpm : float
        Rotational speed [rev/min].
    rho_emulsion : float
        Emulsion density [kg/m3].

    Returns
    -------
    float
        Power draw [W].
    """
    N = rpm / 60.0  # [rev/s]
    D = _get_diameter(geometry)
    Np = geometry.power_number
    return Np * rho_emulsion * N**3 * D**5


def average_dissipation(
    geometry: Geometry,
    rpm: float,
    rho_emulsion: float,
    tank_volume: float | None = None,
) -> float:
    """Average turbulent energy dissipation rate [m2/s3].

    eps_avg = P / (rho * V_tank)

    Parameters
    ----------
    geometry : MixerGeometry | StirrerGeometry
        Equipment geometry.
    rpm : float
        Rotational speed [rev/min].
    rho_emulsion : float
        Emulsion density [kg/m3].
    tank_volume : float, optional
        Override tank volume [m3]. Required when *geometry* is a
        ``StirrerGeometry`` (which does not carry volume information).
        If ``None``, volume is read from the geometry object.

    Returns
    -------
    float
        Volume-averaged dissipation rate [m2/s3].
    """
    P = power_draw(geometry, rpm, rho_emulsion)
    V = tank_volume if tank_volume is not None else _get_tank_volume(geometry)
    if V is None:
        raise ValueError(
            "Cannot determine tank volume from geometry. "
            "Pass tank_volume explicitly when using StirrerGeometry."
        )
    return P / (rho_emulsion * V)


def max_dissipation(
    geometry: Geometry,
    rpm: float,
    rho_emulsion: float,
    tank_volume: float | None = None,
) -> float:
    """Maximum energy dissipation rate in the high-shear zone [m2/s3].

    eps_max = k_eps * eps_avg

    where k_eps is ``geometry.dissipation_ratio`` (typically 5 for open
    impellers, 25-50 for rotor-stators).

    Parameters
    ----------
    geometry : MixerGeometry | StirrerGeometry
        Equipment geometry.
    rpm : float
        Rotational speed [rev/min].
    rho_emulsion : float
        Emulsion density [kg/m3].
    tank_volume : float, optional
        Override tank volume [m3] (see ``average_dissipation``).

    Returns
    -------
    float
        Maximum local dissipation rate [m2/s3].
    """
    return geometry.dissipation_ratio * average_dissipation(
        geometry, rpm, rho_emulsion, tank_volume=tank_volume
    )


# ---------------------------------------------------------------------------
# Reynolds number & corrected power number
# ---------------------------------------------------------------------------

def impeller_reynolds_number(
    rpm: float,
    impeller_diameter: float,
    rho: float,
    mu: float,
) -> float:
    """Impeller Reynolds number [-].

    Re_imp = rho * N * D^2 / mu

    Parameters
    ----------
    rpm : float
        Rotational speed [rev/min].
    impeller_diameter : float
        Impeller diameter [m].
    rho : float
        Fluid density [kg/m3].
    mu : float
        Dynamic viscosity [Pa*s].

    Returns
    -------
    float
        Impeller Reynolds number [-].
    """
    N = rpm / 60.0  # [rev/s]
    return rho * N * impeller_diameter**2 / mu


def power_number_corrected(
    stirrer: StirrerGeometry,
    Re: float,
) -> float:
    """Reynolds-dependent power number [-].

    In the turbulent regime (Re > 10 000) Np equals the geometric constant
    ``stirrer.power_number``.  In the transitional regime (10 < Re < 10 000)
    Np increases following:

        Np(Re) = Np_turb * (1 + K_p / Re^0.5)

    Parameters
    ----------
    stirrer : StirrerGeometry
        Stirrer geometry (provides ``power_number``).
    Re : float
        Impeller Reynolds number [-].

    Returns
    -------
    float
        Corrected power number [-].

    Notes
    -----
    - ASSUMPTION: K_p ~ 25 for pitched-blade turbines (literature range 20-30).
    - ASSUMPTION: K_p ~ 15 for rotor-stator types (higher baseline Np, weaker
      correction).
    """
    Np_turb = stirrer.power_number

    if Re > 10_000:
        return Np_turb

    # ASSUMPTION: K_p depends on stirrer family
    if getattr(stirrer.stirrer_type, 'value', stirrer.stirrer_type) == "pitched_blade":
        K_p = 25.0
    else:
        K_p = 15.0

    return Np_turb * (1.0 + K_p / Re**0.5)


# ---------------------------------------------------------------------------
# Swept volume
# ---------------------------------------------------------------------------

def swept_volume(stirrer: StirrerGeometry) -> float:
    """Volume of the impeller-swept region [m3].

    V_swept = pi/4 * D^2 * blade_height

    Parameters
    ----------
    stirrer : StirrerGeometry
        Stirrer geometry.

    Returns
    -------
    float
        Swept volume [m3].
    """
    D = stirrer.impeller_diameter
    return np.pi / 4.0 * D**2 * stirrer.blade_height


# ---------------------------------------------------------------------------
# Shear rate models
# ---------------------------------------------------------------------------

def metzner_otto_shear_rate(
    rpm: float,
    stirrer: StirrerGeometry,
) -> float:
    """Average shear rate in the vessel via Metzner-Otto correlation [1/s].

    gamma_dot = k_s * N

    Parameters
    ----------
    rpm : float
        Rotational speed [rev/min].
    stirrer : StirrerGeometry
        Stirrer geometry (determines k_s).

    Returns
    -------
    float
        Average shear rate [1/s].

    Notes
    -----
    - ASSUMPTION: k_s ~ 11.5 for pitched-blade turbines (Metzner & Otto 1957).
    - ASSUMPTION: k_s ~ 30 for rotor-stator types (Utomo et al. 2009).
    """
    N = rpm / 60.0  # [rev/s]

    # ASSUMPTION: k_s by stirrer family
    if getattr(stirrer.stirrer_type, 'value', stirrer.stirrer_type) == "pitched_blade":
        k_s = 11.5
    else:
        k_s = 30.0

    return k_s * N


def gap_shear_rate_extended(
    rpm: float,
    stirrer: StirrerGeometry,
) -> float:
    """Shear rate in the rotor-stator gap or average vessel shear rate [1/s].

    For stirrers with a stator (has_stator=True):
        gamma_dot = pi * N * D / delta

    For open impellers (has_stator=False):
        Falls back to ``metzner_otto_shear_rate``.

    Parameters
    ----------
    rpm : float
        Rotational speed [rev/min].
    stirrer : StirrerGeometry
        Stirrer geometry.

    Returns
    -------
    float
        Characteristic shear rate [1/s].
    """
    if not stirrer.has_stator:
        return metzner_otto_shear_rate(rpm, stirrer)

    N = rpm / 60.0
    return np.pi * N * stirrer.impeller_diameter / stirrer.gap_width


# ---------------------------------------------------------------------------
# Legacy gap shear rate (backward compatibility)
# ---------------------------------------------------------------------------

def gap_shear_rate(mixer: MixerGeometry, rpm: float) -> float:
    """Characteristic shear rate in the rotor-stator gap [1/s].

    gamma_dot = pi * N * D / delta
    where delta is the gap width.

    Parameters
    ----------
    mixer : MixerGeometry
        Legacy mixer geometry.
    rpm : float
        Rotational speed [rev/min].

    Returns
    -------
    float
        Gap shear rate [1/s].
    """
    N = rpm / 60.0
    return np.pi * N * mixer.rotor_diameter / mixer.gap_width


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def kolmogorov_length_scale(epsilon: float, nu_c: float) -> float:
    """Kolmogorov microscale [m].

    eta_K = (nu^3 / eps)^(1/4)

    Parameters
    ----------
    epsilon : float
        Energy dissipation rate [m2/s3].
    nu_c : float
        Kinematic viscosity of continuous phase [m2/s].

    Returns
    -------
    float
        Kolmogorov length scale [m].
    """
    if epsilon <= 0:
        return np.inf
    return (nu_c**3 / epsilon)**0.25


def emulsion_density(rho_c: float, rho_d: float, phi_d: float) -> float:
    """Effective emulsion density [kg/m3].

    rho_eff = rho_c * (1 - phi_d) + rho_d * phi_d

    Parameters
    ----------
    rho_c : float
        Continuous phase density [kg/m3].
    rho_d : float
        Dispersed phase density [kg/m3].
    phi_d : float
        Dispersed phase volume fraction [-].

    Returns
    -------
    float
        Effective density [kg/m3].
    """
    return rho_c * (1.0 - phi_d) + rho_d * phi_d
