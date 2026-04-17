"""Node F1-b Phase 1 (v8.1-alpha): L4 cellulose mechanical-properties solver.

Modulus from cellulose volume fraction. The empirical scaling (Zhang
et al. 2020 *Cellulose* 27:1071; Rubinstein & Colby 2003
*Polymer Physics*) for regenerated-cellulose hydrogels is::

    G = K_cell · phi_cellulose_mean^alpha_cell

where:
    phi_cellulose_mean — [-]   volume-averaged polymer fraction (from L2 NIPS)
    K_cell, alpha_cell — empirical prefactor / exponent, from the
                         ``MaterialProperties`` bundle. NaOH/urea default:
                         K_cell = 5 × 10⁵ Pa, alpha_cell = 2.25 (entangled
                         regime). Semi-dilute regime uses alpha_cell = 2.0.

Cellulose gels have no separate L3 (NIPS IS the crosslinking for the
default path), so ``solve_mechanical_cellulose`` takes only the L2
gelation result and the polymer bundle — no ``CrosslinkingResult``
dependency.

Emitted ``MechanicalResult`` follows the same schema as the
agarose/chitosan and alginate paths so downstream consumers (dossier,
optimizer, UI) are platform-agnostic.
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


# Cellulose dry density for the polymer volume fraction → Ogston Kav
# calculation. ~1500 kg/m³ for regenerated cellulose.
_RHO_CELLULOSE = 1500.0  # [kg/m³]


def cellulose_modulus(
    phi_cell_mean: float, K_cell: float, alpha_cell: float,
) -> float:
    """Empirical modulus law for regenerated cellulose hydrogels.

    Returns the shear modulus [Pa]. Zero volume fraction gives 0.0
    deterministically. No implicit temperature dependence in v1 — the
    prefactor ``K_cell`` already bakes in typical processing
    conditions (25 °C for NaOH/urea).
    """
    if phi_cell_mean <= 0.0 or K_cell <= 0.0:
        return 0.0
    return float(K_cell * phi_cell_mean ** alpha_cell)


def solve_mechanical_cellulose(
    params: SimulationParameters,
    props: MaterialProperties,
    gelation: GelationResult,
    R_droplet: Optional[float] = None,
) -> MechanicalResult:
    """Compute cellulose mechanical properties from the L2 NIPS gel.

    Extracts the volume-averaged cellulose fraction ``phi_mean`` from
    the NIPS gelation result's manifest diagnostics (populated by
    ``solve_nips_cellulose``) and returns a SEMI_QUANTITATIVE
    ``MechanicalResult`` with ``G_agarose = 0`` and
    ``G_chitosan = 0``; the cellulose modulus is reported under
    ``G_DN`` to reuse the downstream schema.
    """
    diag = {}
    if gelation.model_manifest is not None:
        diag = gelation.model_manifest.diagnostics or {}
    phi_mean = float(diag.get("phi_mean_final", 0.0))

    K_cell = float(props.K_cell_modulus)
    alpha_cell = float(props.alpha_cell_modulus)

    G_DN = cellulose_modulus(phi_mean, K_cell, alpha_cell)
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
    phi_fiber = max(phi_mean, 0.0)  # NIPS phi is already a volume fraction
    Kav_arr = ogston_kav(rh_arr, props.r_fiber, phi_fiber)

    manifest = ModelManifest(
        model_name="L4.Mechanical.CelluloseZhang2020",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=[
            f"K_cell={K_cell:.2g} Pa",
            f"alpha_cell={alpha_cell}",
            "entangled-regime power law",
            "NIPS-derived phi_cellulose_mean (no covalent crosslinking)",
            "no anisotropy correction",
        ],
        diagnostics={
            "phi_cell_mean": phi_mean,
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
        model_used="cellulose_zhang2020",
        G_DN_lower=0.0,
        G_DN_upper=0.0,
        model_manifest=manifest,
        network_type="physical_entangled",
    )
