"""Node F1-a Phase 2b (v8.0): L4 alginate mechanical-properties solver.

Modulus from egg-box junction density. The empirical scaling (Kong et
al. 2004 *Macromolecules* 37:6838) for Ca²⁺-alginate gels is::

    G_DN = K_alg · (c_alginate · f_G)^n_alg · (X_mean / X_max)

where:
    c_alginate       — [kg/m³] alginate concentration
    f_G              — [-]     guluronate fraction (G-block)
    X_mean           — [mol/m³] volume-averaged egg-box junction density
                                 (from the L2 ionic_ca solver)
    X_max            — [mol/m³] stoichiometric ceiling (½·c_alginate·f_G)
    K_alg, n_alg     — empirical prefactor / exponent, from the
                       ``MaterialProperties`` bundle (Node F1-a defaults
                       K_alg=30 kPa, n_alg=2.0).

Alginate has no separate L3 (ionic gelation IS the crosslinking), so
``solve_mechanical_alginate`` takes only the gelation result and the
polymer bundle — no ``CrosslinkingResult`` dependency.

Emitted ``MechanicalResult`` follows the same schema as the
agarose/chitosan path so downstream consumers (dossier, optimizer,
UI) are platform-agnostic.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from ..datatypes import (
    GelationResult,
    MaterialProperties,
    MechanicalResult,
    ModelEvidenceTier,
    ModelManifest,
    SimulationParameters,
)
from .solver import effective_youngs_modulus, hertz_contact, ogston_kav

logger = logging.getLogger(__name__)


# Alginate dry density used for the polymer volume fraction → Ogston Kav
# calculation. Matches the repeat-unit density commonly cited
# (≈1600 kg/m³ for Na-alginate).
_RHO_ALGINATE = 1600.0  # [kg/m³]


def alginate_modulus(
    c_alginate: float,
    f_G: float,
    X_mean: float,
    K_alg: float,
    n_alg: float,
) -> float:
    """Kong et al. 2004 empirical modulus law for Ca²⁺-alginate gel.

    Returns the shear modulus [Pa]. Zero alginate, zero guluronate, or
    zero crosslink density all give 0.0 deterministically.
    """
    if c_alginate <= 0.0 or f_G <= 0.0 or X_mean <= 0.0:
        return 0.0
    # Repeat-unit molar mass used for the stoichiometric ceiling. A
    # single Ca²⁺ binds 2 guluronate residues → one junction.
    M_repeat = 0.198  # [kg/mol] (matches the ionic_ca solver)
    G_total = (c_alginate / M_repeat) * f_G                  # [mol/m³]
    X_max = max(0.5 * G_total, 1e-300)                       # [mol/m³]
    conv = max(0.0, min(1.0, X_mean / X_max))
    return float(K_alg * (c_alginate * f_G) ** n_alg * conv)


def solve_mechanical_alginate(
    params: SimulationParameters,
    props: MaterialProperties,
    gelation: GelationResult,
    R_droplet: Optional[float] = None,
) -> MechanicalResult:
    """Compute alginate mechanical properties from the L2 ionic-Ca gel.

    Extracts the mean egg-box junction density ``X_mean`` from the
    gelation result's manifest diagnostics (populated by
    ``solve_ionic_ca_gelation``) and returns a SEMI_QUANTITATIVE
    ``MechanicalResult`` with ``G_agarose = 0`` and
    ``G_chitosan = 0`` (alginate has no agarose/chitosan network).
    """
    # Pull X_mean from the ionic solver's diagnostic bundle.
    diag = {}
    if gelation.model_manifest is not None:
        diag = gelation.model_manifest.diagnostics or {}
    X_mean = float(diag.get("X_mean_final", 0.0))

    c_alg = float(
        params.formulation.c_alginate
        if params.formulation.c_alginate > 0.0
        else params.formulation.c_agarose
    )
    f_G = float(props.f_guluronate)
    K_alg = float(props.K_alg_modulus)
    n_alg = float(props.n_alg_modulus)

    G_DN = alginate_modulus(c_alg, f_G, X_mean, K_alg, n_alg)
    E_star = effective_youngs_modulus(G_DN)

    # Bead radius for Hertz contact + Ogston Kav arrays
    if R_droplet is not None and R_droplet > 0.0:
        R = R_droplet
    elif gelation.L_domain > 0.0:
        R = gelation.L_domain / 2.0
    else:
        R = gelation.r_grid[-1] if gelation.r_grid.size else 50e-6

    delta_arr, F_arr = hertz_contact(E_star, R)

    rh_arr = np.logspace(np.log10(1e-9), np.log10(50e-9), 50)
    # Alginate polymer volume fraction: c / rho_alginate
    phi_fiber = max(c_alg / _RHO_ALGINATE, 0.0)
    Kav_arr = ogston_kav(rh_arr, props.r_fiber, phi_fiber)

    manifest = ModelManifest(
        model_name="L4.Mechanical.AlginateKong2004",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=[
            f"K_alg={K_alg:.2g} Pa",
            f"n_alg={n_alg}",
            f"f_G={f_G}",
            "ionic egg-box junctions (no covalent network)",
            "empirical modulus ∝ (c_alg·f_G)^n with conversion factor X_mean/X_max",
        ],
        diagnostics={
            "X_mean": X_mean,
            "c_alginate_kgm3": c_alg,
            "G_DN_Pa": G_DN,
        },
    )

    return MechanicalResult(
        G_agarose=0.0,
        G_chitosan=0.0,
        G_DN=float(G_DN),
        E_star=float(E_star),
        delta_array=delta_arr,
        F_array=F_arr,
        rh_array=rh_arr,
        Kav_array=Kav_arr,
        pore_size_mean=float(gelation.pore_size_mean),
        xi_mesh=float(gelation.pore_size_mean),  # no separate mesh in alginate
        model_used="alginate_kong2004",
        G_DN_lower=0.0,
        G_DN_upper=0.0,
        model_manifest=manifest,
        network_type="ionic_reinforced",
    )
