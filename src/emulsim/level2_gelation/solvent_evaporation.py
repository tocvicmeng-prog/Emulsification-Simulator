"""Node F1-c Phase 1 (v8.2-alpha): L2 solvent-evaporation solver for
PLGA microspheres (O/W single-emulsion route).

Simulates the DCM (or ethyl acetate) depletion inside a dispersed
PLGA/DCM droplet suspended in an aqueous continuous phase. DCM
diffuses from the droplet interior to the droplet surface, partitions
into the water phase (approximated as a perfect sink), and eventually
evaporates to the gas headspace — from the droplet's perspective, the
outer boundary is simply a low-DCM Dirichlet condition.

Governing equation
------------------
1D spherical Fickian diffusion on ``r ∈ [0, R_droplet]``::

    ∂phi_DCM/∂t = D · (1/r²) · ∂/∂r(r² · ∂phi_DCM/∂r)

    phi_PLGA(r, t) = 1 − phi_DCM(r, t)

Boundary conditions
-------------------
- ``∂phi_DCM/∂r(0, t) = 0``                    (spherical symmetry)
- ``phi_DCM(R, t) = phi_DCM_eq``               (aqueous-phase sink;
                                                ~0.005 for DCM/water)

Initial condition
-----------------
``phi_DCM(r, 0) = 1 − phi_PLGA_0`` (uniform).

Approximations (Phase 1)
------------------------
- **Fixed droplet radius**: the droplet should shrink as DCM leaves
  (final volume = ``phi_PLGA_0 × V_0``), but that requires a
  Lagrangian / ALE moving-boundary scheme. v1 keeps ``R`` constant
  and reports ``R_final = R_0 · phi_PLGA_0^(1/3)`` as a
  mass-conservation post-processing diagnostic.
- **Constant D**: real DCM diffusion in PLGA/DCM follows
  ``D(phi) = D_0 · exp(−β · phi_PLGA)``; v1 uses a single ``D_DCM``
  representative of the early-to-mid regime (``phi_PLGA < 0.6``).
  Late-time (``phi_PLGA > 0.8``) dynamics will be slower than
  predicted — vitrification times are order-of-magnitude right,
  not quantitative.

See :file:`docs/f1c_plga_protocol.md` for the full scientific basis.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from ..datatypes import (
    GelationResult,
    MaterialProperties,
    ModelEvidenceTier,
    ModelManifest,
    SimulationParameters,
)

logger = logging.getLogger(__name__)


# Solidification threshold: volume fraction above which the local
# polymer is treated as glassy / vitrified. ~0.85 corresponds to the
# PLGA/DCM composition where the mixture T_g crosses ambient T.
_VITRIFICATION_PHI = 0.85


def _spherical_laplacian_dirichlet(
    f: np.ndarray, bath: float,
    r_face: np.ndarray, V_shells: np.ndarray, dr: float,
) -> np.ndarray:
    """Spherical Laplacian (1/r²) d/dr(r² df/dr) with Dirichlet outer
    boundary (ghost value = ``bath``) and symmetry at ``r = 0``.

    Shares the same numerical convention as the NIPS-cellulose solver
    (node-centred uniform grid, face-flux divergence). Kept local
    here rather than importing from ``nips_cellulose`` to avoid a
    cross-module dependency between two independent platforms.
    """
    n = f.size
    grad = np.empty(n)
    grad[:-1] = (f[1:] - f[:-1]) / dr
    # Outer face at r = R: ghost value = bath, half-cell distance
    grad[-1] = (bath - f[-1]) / (0.5 * dr)

    lap = np.empty(n)
    # Cell 0: inner face at r = 0 has zero flux (symmetry)
    lap[0] = (r_face[0] ** 2) * grad[0] / V_shells[0]
    for i in range(1, n - 1):
        lap[i] = (
            (r_face[i] ** 2) * grad[i] - (r_face[i - 1] ** 2) * grad[i - 1]
        ) / V_shells[i]
    lap[-1] = (
        (r_face[-1] ** 2) * grad[-1] - (r_face[-2] ** 2) * grad[-2]
    ) / V_shells[-1]
    return lap


def solve_solvent_evaporation(
    params: SimulationParameters,
    props: MaterialProperties,
    *,
    R_droplet: float,
    n_r: int = 40,
    time: Optional[float] = None,
    rtol: float = 1e-6,
    atol: float = 1e-10,
) -> GelationResult:
    """Solve 1D spherical DCM depletion for a PLGA/DCM droplet.

    Parameters
    ----------
    params : SimulationParameters
        Uses ``formulation.phi_PLGA_0`` as the initial polymer volume
        fraction and ``formulation.t_crosslink`` as the total
        evaporation time.
    props : MaterialProperties
        Reads ``D_DCM_plga`` and ``phi_DCM_eq``. Recommended: call
        ``properties.plga_defaults.apply_preset(props, "50_50")``
        to bulk-load a consistent grade.
    R_droplet : float
        Droplet radius [m].
    n_r : int
        Radial grid cells.
    time : float, optional
        Total evaporation time [s]. Defaults to
        ``params.formulation.t_crosslink``.
    rtol, atol : float
        ``solve_ivp`` tolerances.

    Returns
    -------
    GelationResult
        Schema-compatible with the other L2 solvers. Manifest tagged
        ``L2.Gelation.SolventEvaporationPLGA`` at
        SEMI_QUANTITATIVE tier with diagnostics
        ``phi_plga_mean_final``, ``t_vitrification``,
        ``skin_thickness_proxy``, ``core_porosity_proxy``,
        ``R_shrunk_m``.
    """
    if R_droplet <= 0:
        raise ValueError(f"R_droplet must be > 0, got {R_droplet}")
    if n_r < 8:
        raise ValueError(f"n_r must be >= 8 for stability, got {n_r}")

    t_end = float(time) if time is not None else params.formulation.t_crosslink
    if t_end < 0:
        raise ValueError(f"time must be >= 0, got {t_end}")

    phi_0 = float(params.formulation.phi_PLGA_0)
    if phi_0 < 0 or phi_0 > 1:
        raise ValueError(f"phi_PLGA_0 must be in [0, 1], got {phi_0}")
    if phi_0 <= 0:
        return _zero_plga_result(R_droplet, n_r, t_end)

    D_DCM = float(props.D_DCM_plga)
    phi_eq = float(props.phi_DCM_eq)
    if D_DCM <= 0:
        raise ValueError(f"D_DCM_plga must be > 0, got {D_DCM}")
    if phi_eq < 0 or phi_eq > 1:
        raise ValueError(f"phi_DCM_eq must be in [0, 1], got {phi_eq}")

    # Grid
    dr = R_droplet / n_r
    r = np.linspace(0.5 * dr, R_droplet - 0.5 * dr, n_r)
    r_face = np.linspace(dr, R_droplet, n_r)
    r_in_arr = np.concatenate([[0.0], r_face[:-1]])
    V_shells = (r_face ** 3 - r_in_arr ** 3) / 3.0
    V_total = float(V_shells.sum())

    # Initial condition: uniform DCM = 1 - phi_0 everywhere
    phi_DCM_0 = 1.0 - phi_0
    y0 = np.full(n_r, phi_DCM_0)

    def rhs(_t, y):
        return D_DCM * _spherical_laplacian_dirichlet(
            y, phi_eq, r_face, V_shells, dr,
        )

    # Dense output to probe t_vitrification
    sol = solve_ivp(
        rhs, (0.0, t_end), y0,
        method="BDF", rtol=rtol, atol=atol,
        dense_output=True,
    )

    if not sol.success:
        logger.warning("Solvent evaporation solver failed: %s", sol.message)

    phi_DCM_final = np.clip(sol.y[:, -1], 0.0, 1.0)
    phi_PLGA_final = 1.0 - phi_DCM_final

    phi_plga_mean = float((phi_PLGA_final * V_shells).sum() / V_total)
    phi_dcm_mean = float((phi_DCM_final * V_shells).sum() / V_total)

    # Mass-conservation post-processing: the real physical sphere
    # shrinks to volume V_0 · phi_0 / phi_plga_mean. R_shrunk below
    # assumes a uniform dense microsphere at phi_plga_mean.
    R_shrunk = R_droplet * (phi_0 / max(phi_plga_mean, 1e-12)) ** (1.0 / 3.0)
    R_shrunk = float(min(R_shrunk, R_droplet))

    # Skin-thickness proxy: radial distance from the outer surface to
    # the first interior cell whose phi_PLGA has crossed the
    # vitrification threshold. Returns 0 if no cell has vitrified.
    vitrified_mask = phi_PLGA_final >= _VITRIFICATION_PHI
    if vitrified_mask.any():
        innermost_vitrified_r = float(r[vitrified_mask].min())
        skin_thickness = float(R_droplet - innermost_vitrified_r)
    else:
        skin_thickness = 0.0

    # Core porosity proxy: non-polymer fraction at r = r[0] (centre
    # cell). High core-porosity = trapped DCM after skin formed.
    core_porosity_proxy = float(1.0 - phi_PLGA_final[0])

    # t_vitrification: first time where volume-averaged phi_PLGA
    # crosses the threshold. Use dense solver output.
    t_vitrification = _find_mean_vitrification_time(
        sol, V_shells, V_total, t_end,
    )

    # Porosity (final) = residual DCM fraction (for fixed-R; for real
    # shrunk sphere it's much lower, but the fixed-R porosity is still
    # a useful process-design KPI).
    porosity = float(max(0.0, min(1.0, phi_dcm_mean)))

    # Characteristic pore size proxy: cell spacing × residual DCM
    pore_size_mean = max(dr * max(phi_dcm_mean, 0.01), 1.0e-9)
    pore_size_std = 0.3 * pore_size_mean
    chord_samples = np.array(
        [pore_size_mean * 0.6, pore_size_mean, pore_size_mean * 1.5]
    )

    manifest = ModelManifest(
        model_name="L2.Gelation.SolventEvaporationPLGA",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=[
            "1D spherical Fickian diffusion of DCM",
            "fixed droplet radius (no moving boundary)",
            "constant D_DCM (concentration-independent)",
            "perfect aqueous sink at r=R (Dirichlet at phi_DCM_eq)",
            f"D_DCM={D_DCM:.2e} m²/s",
            f"phi_DCM_eq={phi_eq}",
            f"phi_PLGA_0={phi_0}",
            f"vitrification threshold = {_VITRIFICATION_PHI}",
        ],
        diagnostics={
            "phi_plga_mean_final": phi_plga_mean,
            "phi_dcm_mean_final": phi_dcm_mean,
            "t_vitrification": t_vitrification,
            "skin_thickness_proxy": skin_thickness,
            "core_porosity_proxy": core_porosity_proxy,
            "R_shrunk_m": R_shrunk,
            "solver_success": bool(sol.success),
            "n_steps": int(sol.t.size),
        },
    )

    return GelationResult(
        r_grid=r,
        phi_field=phi_PLGA_final,
        pore_size_mean=pore_size_mean,
        pore_size_std=pore_size_std,
        pore_size_distribution=chord_samples,
        porosity=porosity,
        alpha_final=phi_plga_mean,      # "conversion" = densification
        char_wavelength=pore_size_mean,
        T_history=None,
        phi_snapshots=None,
        L_domain=R_droplet,
        grid_spacing=dr,
        bicontinuous_score=0.0,
        anisotropy=0.0,
        connectivity=1.0,
        chord_skewness=0.0,
        model_tier="mechanistic",
        model_manifest=manifest,
    )


def _find_mean_vitrification_time(
    sol, V_shells: np.ndarray, V_total: float, t_end: float,
) -> float:
    """Bisect-ish scan of dense solver output for the first time where
    the volume-averaged phi_PLGA crosses the vitrification threshold.

    Returns ``t_end`` if the mean never crosses.
    """
    # Scan 200 points in [0, t_end] using dense output. 200 is plenty
    # even for short t_end (probe granularity = t_end/200); refine
    # further with a Brent root-find if callers need sub-probe accuracy.
    probe = np.linspace(0.0, t_end, 201)
    prev_mean = 1.0 - sol.y[0, 0]  # initial phi_PLGA mean (uniform)
    for t_k in probe[1:]:
        y_k = sol.sol(t_k)
        phi_plga_k = 1.0 - np.clip(y_k, 0.0, 1.0)
        mean_k = float((phi_plga_k * V_shells).sum() / V_total)
        if mean_k >= _VITRIFICATION_PHI and prev_mean < _VITRIFICATION_PHI:
            return float(t_k)
        prev_mean = mean_k
    return float(t_end)


def _zero_plga_result(
    R_droplet: float, n_r: int, t_end: float,
) -> GelationResult:
    """Empty-result fallback for zero-PLGA input."""
    r = np.linspace(0.0, R_droplet, n_r)
    return GelationResult(
        r_grid=r,
        phi_field=np.zeros(n_r),
        pore_size_mean=R_droplet,
        pore_size_std=0.0,
        pore_size_distribution=np.array([R_droplet]),
        porosity=1.0,
        alpha_final=0.0,
        char_wavelength=R_droplet,
        T_history=None,
        phi_snapshots=None,
        L_domain=R_droplet,
        grid_spacing=R_droplet / n_r,
        bicontinuous_score=0.0,
        anisotropy=0.0,
        connectivity=0.0,
        chord_skewness=0.0,
        model_tier="mechanistic",
        model_manifest=ModelManifest(
            model_name="L2.Gelation.SolventEvaporationPLGA",
            evidence_tier=ModelEvidenceTier.UNSUPPORTED,
            assumptions=["zero PLGA input -> no microsphere"],
            diagnostics={"phi_PLGA_0": 0.0},
        ),
    )
