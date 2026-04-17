"""Node F1-c Phase 1 (v8.2-alpha): L4 PLGA mechanical-properties solver.

Modulus from polymer volume fraction via a Gibson-Ashby-style
closed-cell foam scaling::

    G = G_glassy · phi_PLGA_mean ^ n_plga

with ``G_glassy ≈ 7 × 10⁸ Pa`` for PLGA 50:50 (from E ≈ 2 GPa,
ν ≈ 0.42; Park et al. 1998) and ``n_plga = 2.0`` for the
closed-cell-foam scaling. A dense microsphere (``phi_PLGA_mean → 1``)
recovers ``G = G_glassy``; a 50 %-porous microsphere yields
``G = 0.25 · G_glassy``.

PLGA is a glassy amorphous polymer below ``T_g`` (40-60 °C depending
on L:G ratio and M_n). Within the normal process / storage
temperature window the modulus is essentially independent of crosslink
density (there are none) — only the polymer volume fraction and
grade-dependent prefactor set the mechanics. Above ``T_g`` the
modulus drops 3-4 decades as the polymer becomes rubbery; this Phase
1 solver does NOT model that transition (flagged in manifest
assumptions).

Emitted ``MechanicalResult`` follows the same schema as every other
L4 solver so downstream consumers (dossier, optimiser, UI) stay
platform-agnostic.
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


# PLGA dry density for downstream Ogston / diffusion calculations.
# ~1200-1300 kg/m³; PLGA 50:50 = 1260 typical.
_RHO_PLGA = 1260.0  # [kg/m³]


def plga_modulus(
    phi_plga_mean: float, G_glassy: float, n_plga: float,
) -> float:
    """Gibson-Ashby closed-cell-foam modulus for PLGA microspheres.

    Returns the shear modulus [Pa]. Zero polymer volume fraction and
    zero prefactor both give 0.0 deterministically.
    """
    if phi_plga_mean <= 0.0 or G_glassy <= 0.0:
        return 0.0
    return float(G_glassy * phi_plga_mean ** n_plga)


def solve_mechanical_plga(
    params: SimulationParameters,
    props: MaterialProperties,
    gelation: GelationResult,
    R_droplet: Optional[float] = None,
) -> MechanicalResult:
    """Compute PLGA mechanical properties from the L2 evaporation result.

    Extracts the volume-averaged polymer fraction ``phi_plga_mean``
    from the L2 manifest diagnostics (populated by
    ``solve_solvent_evaporation``) and returns a SEMI_QUANTITATIVE
    ``MechanicalResult`` with ``G_agarose = 0`` and
    ``G_chitosan = 0``; the PLGA modulus is reported under ``G_DN``
    to reuse the downstream schema.
    """
    diag = {}
    if gelation.model_manifest is not None:
        diag = gelation.model_manifest.diagnostics or {}
    phi_mean = float(diag.get("phi_plga_mean_final", 0.0))

    G_glassy = float(props.G_glassy_plga)
    n_plga = float(props.n_plga_modulus)

    G_DN = plga_modulus(phi_mean, G_glassy, n_plga)
    E_star = effective_youngs_modulus(G_DN)

    if R_droplet is not None and R_droplet > 0.0:
        R = R_droplet
    elif gelation.L_domain > 0.0:
        R = gelation.L_domain / 2.0
    else:
        R = gelation.r_grid[-1] if gelation.r_grid.size else 50e-6

    delta_arr, F_arr = hertz_contact(E_star, R)

    rh_arr = np.logspace(np.log10(1e-9), np.log10(50e-9), 50)
    phi_fiber = max(phi_mean, 0.0)
    Kav_arr = ogston_kav(rh_arr, props.r_fiber, phi_fiber)

    manifest = ModelManifest(
        model_name="L4.Mechanical.PLGAGibsonAshby",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=[
            f"G_glassy={G_glassy:.2g} Pa",
            f"n_plga={n_plga}",
            "closed-cell-foam Gibson-Ashby scaling",
            "glassy amorphous polymer below T_g (no rubbery transition)",
            "no explicit T dependence in v1 (T << T_g assumed)",
        ],
        diagnostics={
            "phi_plga_mean": phi_mean,
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
        xi_mesh=float(gelation.pore_size_mean),
        model_used="plga_gibson_ashby",
        G_DN_lower=0.0,
        G_DN_upper=0.0,
        model_manifest=manifest,
        network_type="glassy_polymer",
    )
