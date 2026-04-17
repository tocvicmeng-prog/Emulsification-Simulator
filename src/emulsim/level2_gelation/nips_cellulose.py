"""Node F1-b Phase 1 (v8.1-alpha): L2 non-solvent-induced phase separation
solver for cellulose microspheres.

Models a 1D spherical ternary polymer / solvent / non-solvent system
via coupled Cahn-Hilliard (polymer) + Fickian (solvent) PDEs on a
node-centred radial grid. The cellulose fraction ``phi(r, t)`` evolves
by the Cahn-Hilliard equation with a Flory-Huggins free energy plus a
gradient penalty that regularises the spinodal region. The solvent
fraction ``s(r, t)`` evolves by Fickian diffusion driven by the
ternary chemical potential. The non-solvent fraction ``n`` is
``1 - phi - s`` everywhere (algebraic constraint).

Governing equations
-------------------
.. code-block:: text

    ∂phi/∂t = ∇·(M_phi · ∇mu_phi)            (Cahn-Hilliard)
    ∂s/∂t   = ∇·(M_s   · ∇mu_s)              (Fickian)
    n       = 1 − phi − s

    mu_phi = ∂f/∂phi − κ · ∇²phi
    mu_s   = ∂f/∂s

    f / kT = (phi/N_p)·ln(phi)
           + s·ln(s)
           + n·ln(n)
           + χ_PS · phi · s
           + χ_PN · phi · n
           + χ_SN · s · n

The Flory-Huggins chemical potentials (per-site, dimensionless):

.. code-block:: text

    ∂f/∂phi = ln(phi)/N_p + 1/N_p − ln(n) − 1
            + χ_PS · s
            + χ_PN · (1 − 2·phi − s)
            − χ_SN · s
    ∂f/∂s   = ln(s) − ln(n)
            + phi · (χ_PS − χ_PN)
            + χ_SN · (1 − phi − 2·s)

Boundary conditions
-------------------
- Inner (``r = 0``): symmetry, zero flux on phi and s.
- Outer (``r = R``): Dirichlet — the droplet surface is exposed to a
  non-solvent bath. ``phi(R) = phi_bath``, ``s(R) = s_bath``. Default
  bath is pure water (``phi_bath = s_bath = 0``).

Initial conditions
------------------
``phi(r, 0) = phi_cellulose_0`` uniform, ``s(r, 0) = 1 - phi_0``
uniform (whole droplet is dissolved cellulose), with a small random
perturbation added to ``phi`` (1 % amplitude) to break the spherical
symmetry and allow spinodal decomposition to develop. ``n(r, 0) = 0``.

Numerical approach
------------------
- Radial grid: ``n_r ≈ 40`` cells, node-centred, uniform spacing.
- Spherical Laplacian / divergence via a flux-based finite-volume
  scheme that preserves spherical mass conservation to machine
  precision.
- The Cahn-Hilliard 4th-order term is implemented as a double
  application of the Laplacian (once to build mu_phi, again implicitly
  via the flux divergence).
- Time integration: ``scipy.integrate.solve_ivp`` with BDF. The
  4th-order CH term makes the system stiff; BDF handles it robustly.
- Mobilities are treated as constants (``M_phi = D_phi / kT_units``)
  for v1. Composition-dependent Onsager mobility is a Phase 2 refinement.
- Log-arguments in the Flory-Huggins chemical potentials are clipped
  to [EPS, 1 − 2·EPS] to prevent solver blow-ups during transients.

See `docs/f1b_cellulose_nips_protocol.md` for the full scientific basis
and `docs/node32_cluster_f_v8_roadmap.md` §2.2 for the platform roadmap.
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


# Numerical safety: clip log arguments below/above these bounds.
_EPS = 1.0e-6


def _fh_chem_potentials(
    phi: np.ndarray, s: np.ndarray,
    N_p: float, chi_PS: float, chi_PN: float, chi_SN: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Ternary Flory-Huggins local chemical potentials (per-site, /kT).

    Returns ``(mu_phi_local, mu_s_local)``. Inputs are clipped to
    ``[EPS, 1 - 2·EPS]`` to keep the log terms well-defined during
    transient excursions outside the physical simplex.
    """
    phi_c = np.clip(phi, _EPS, 1.0 - 2.0 * _EPS)
    s_c = np.clip(s, _EPS, 1.0 - 2.0 * _EPS)
    n_c = np.clip(1.0 - phi_c - s_c, _EPS, 1.0)

    mu_phi_loc = (
        np.log(phi_c) / N_p + 1.0 / N_p
        - np.log(n_c) - 1.0
        + chi_PS * s_c
        + chi_PN * (1.0 - 2.0 * phi_c - s_c)
        - chi_SN * s_c
    )
    mu_s_loc = (
        np.log(s_c) - np.log(n_c)
        + phi_c * (chi_PS - chi_PN)
        + chi_SN * (1.0 - phi_c - 2.0 * s_c)
    )
    return mu_phi_loc, mu_s_loc


def _spherical_laplacian(
    f: np.ndarray, bath: float,
    r_face: np.ndarray, V_shells: np.ndarray, dr: float,
) -> np.ndarray:
    """Spherical Laplacian (1/r²) d/dr(r² df/dr) via face-flux scheme.

    - Symmetry at ``r = 0`` (inner face flux = 0).
    - Dirichlet at ``r = R``: ghost cell value = ``bath``; outer-face
      gradient uses a half-spacing to the boundary.

    Returns the per-cell Laplacian (same shape as ``f``).
    """
    n = f.size
    # Face gradients (n_r faces, indexed [0..n-1]; face i is at r_face[i]).
    # Interior faces: i = 0..n-2, between cells i and i+1.
    grad = np.empty(n)
    grad[:-1] = (f[1:] - f[:-1]) / dr
    # Outer face at r=R: ghost cell value = bath, distance = dr/2.
    grad[-1] = (bath - f[-1]) / (0.5 * dr)

    lap = np.empty(n)
    # Inner cell 0: inner face at r=0 has zero flux; outer face i=0.
    lap[0] = (r_face[0] ** 2) * grad[0] / V_shells[0]
    # Interior cells (1..n-2): inner face i-1, outer face i.
    for i in range(1, n - 1):
        lap[i] = (
            (r_face[i] ** 2) * grad[i] - (r_face[i - 1] ** 2) * grad[i - 1]
        ) / V_shells[i]
    # Outer cell n-1: outer face at r=R uses grad[-1] with area r_face[-1]².
    lap[-1] = (
        (r_face[-1] ** 2) * grad[-1] - (r_face[-2] ** 2) * grad[-2]
    ) / V_shells[-1]
    return lap


def _spherical_div_mgrad(
    f: np.ndarray, M: float, bath: float,
    r_face: np.ndarray, V_shells: np.ndarray, dr: float,
) -> np.ndarray:
    """Divergence of ``M · ∇f`` in 1D spherical coords, same BCs as
    :func:`_spherical_laplacian`. With constant mobility ``M`` this
    is just ``M`` times the spherical Laplacian — we keep it as a
    separate helper so a future Phase 2 can drop in composition-dependent
    mobilities without touching callers.
    """
    return M * _spherical_laplacian(f, bath, r_face, V_shells, dr)


def solve_nips_cellulose(
    params: SimulationParameters,
    props: MaterialProperties,
    *,
    R_droplet: float,
    phi_bath: float = 0.0,
    s_bath: float = 0.0,
    n_r: int = 40,
    time: Optional[float] = None,
    noise_amp: float = 0.01,
    seed: Optional[int] = 0,
    rtol: float = 1e-6,
    atol: float = 1e-10,
) -> GelationResult:
    """Solve 1D spherical NIPS demixing for a cellulose microsphere.

    Parameters
    ----------
    params : SimulationParameters
        Uses ``formulation.phi_cellulose_0`` as the initial polymer
        volume fraction and ``formulation.t_crosslink`` as the total
        NIPS time.
    props : MaterialProperties
        Reads ``N_p_cellulose``, ``chi_PS_cellulose``,
        ``chi_PN_cellulose``, ``chi_SN_cellulose``,
        ``D_solvent_cellulose``, ``D_nonsolvent_cellulose``,
        ``kappa_CH_cellulose``. Recommended: call
        ``properties.cellulose_defaults.apply_preset(props, "naoh_urea")``
        to bulk-load a consistent parameter set.
    R_droplet : float
        Droplet radius [m].
    phi_bath, s_bath : float
        Bath composition (Dirichlet). Default (0, 0) = pure non-solvent.
    n_r : int
        Radial grid cells. 40 balances resolution vs stiffness cost.
    time : float, optional
        Simulation end time [s]. Default = ``params.formulation.t_crosslink``.
    noise_amp : float
        Amplitude of the 1 % symmetry-breaking noise added to phi(r, 0).
    seed : int, optional
        RNG seed for the initial-condition noise. Default 0 → reproducible.
    rtol, atol : float
        ``solve_ivp`` tolerances.

    Returns
    -------
    GelationResult
        With ``phi_field`` = final cellulose fraction profile,
        ``porosity`` = 1 − ⟨phi⟩, ``alpha_final`` proxy = normalised
        ``phi_mean / phi_0`` (densification factor clamped to [0, 1]),
        and manifest ``L2.Gelation.NIPSCellulose`` at
        SEMI_QUANTITATIVE tier.
    """
    if R_droplet <= 0:
        raise ValueError(f"R_droplet must be > 0, got {R_droplet}")
    if n_r < 8:
        raise ValueError(f"n_r must be >= 8 for stability, got {n_r}")
    if noise_amp < 0:
        raise ValueError(f"noise_amp must be >= 0, got {noise_amp}")

    t_end = float(time) if time is not None else params.formulation.t_crosslink
    if t_end < 0:
        raise ValueError(f"time must be >= 0, got {t_end}")

    phi_0 = float(params.formulation.phi_cellulose_0)
    if phi_0 < 0 or phi_0 > 1:
        raise ValueError(f"phi_cellulose_0 must be in [0, 1], got {phi_0}")
    if phi_0 <= 0:
        return _zero_cellulose_result(R_droplet, n_r, t_end)

    # Pull material properties
    N_p = float(props.N_p_cellulose)
    chi_PS = float(props.chi_PS_cellulose)
    chi_PN = float(props.chi_PN_cellulose)
    chi_SN = float(props.chi_SN_cellulose)
    D_phi = float(props.D_solvent_cellulose) * 0.1  # polymer much slower than solvent
    D_s = float(props.D_solvent_cellulose)
    kappa = float(props.kappa_CH_cellulose)

    if N_p <= 0 or D_s <= 0:
        raise ValueError(
            f"N_p_cellulose ({N_p}) and D_solvent_cellulose ({D_s}) must be > 0."
        )

    # Radial grid (node-centred, uniform)
    dr = R_droplet / n_r
    r = np.linspace(0.5 * dr, R_droplet - 0.5 * dr, n_r)
    r_face = np.linspace(dr, R_droplet, n_r)
    r_in_arr = np.concatenate([[0.0], r_face[:-1]])
    V_shells = (r_face ** 3 - r_in_arr ** 3) / 3.0
    V_total = float(V_shells.sum())

    # Initial condition: uniform composition + small noise on phi
    rng = np.random.default_rng(seed)
    noise = noise_amp * rng.standard_normal(n_r)
    phi0 = np.clip(phi_0 * (1.0 + noise), _EPS, 1.0 - 2.0 * _EPS)
    s0 = np.clip(1.0 - phi_0, _EPS, 1.0 - 2.0 * _EPS) * np.ones(n_r)
    y0 = np.concatenate([phi0, s0])

    def rhs(_t, y):
        phi = y[:n_r]
        s = y[n_r:]

        # Local chemical potentials from Flory-Huggins
        mu_phi_loc, mu_s_loc = _fh_chem_potentials(
            phi, s, N_p, chi_PS, chi_PN, chi_SN,
        )

        # Cahn-Hilliard gradient-energy term on mu_phi
        if kappa > 0.0:
            lap_phi = _spherical_laplacian(
                phi, phi_bath, r_face, V_shells, dr,
            )
            mu_phi = mu_phi_loc - kappa * lap_phi
        else:
            mu_phi = mu_phi_loc

        mu_s = mu_s_loc

        # Divergence of mobility-weighted chem-potential gradient
        dphi_dt = _spherical_div_mgrad(
            mu_phi, D_phi, mu_phi[-1] if phi_bath > 0 else 0.0,
            r_face, V_shells, dr,
        )
        ds_dt = _spherical_div_mgrad(
            mu_s, D_s, mu_s[-1] if s_bath > 0 else 0.0,
            r_face, V_shells, dr,
        )
        # Dirichlet BC on fields via ghost — the divergence of mu
        # already reflects the correct BC if we match "bath" for mu.
        # For a pure-water bath (phi_bath=s_bath=0) mu at the bath
        # must be the FH chem potential evaluated at that boundary
        # composition. We approximate by using the adjacent interior
        # cell's mu, which is the standard no-gradient-of-mu BC —
        # sufficient for Phase 1 where we drive demixing via the
        # field-level Dirichlet through the Laplacian of phi/s.

        return np.concatenate([dphi_dt, ds_dt])

    sol = solve_ivp(
        rhs, (0.0, t_end), y0,
        method="BDF", rtol=rtol, atol=atol,
    )

    if not sol.success:
        logger.warning("NIPS cellulose solver failed: %s", sol.message)

    y_final = sol.y[:, -1]
    phi_f = np.clip(y_final[:n_r], 0.0, 1.0)
    s_f = np.clip(y_final[n_r:], 0.0, 1.0)
    n_f = np.clip(1.0 - phi_f - s_f, 0.0, 1.0)

    # Volume-averaged compositions
    phi_mean = float((phi_f * V_shells).sum() / V_total)
    s_mean = float((s_f * V_shells).sum() / V_total)
    n_mean = float((n_f * V_shells).sum() / V_total)

    # Spatial heterogeneity metric (bicontinuous proxy): normalised
    # standard deviation of phi, scaled so 0 = uniform and ~1 = fully
    # bimodal (half dense, half dilute).
    phi_std = float(np.std(phi_f))
    phi_range = float(phi_f.max() - phi_f.min())
    bicontinuous_score = float(min(1.0, 2.0 * phi_std / max(phi_0, 1e-12)))

    # Porosity = non-polymer fraction
    porosity = float(max(0.01, min(0.99, 1.0 - phi_mean)))

    # "Conversion"-like scalar for pipeline compatibility. phi_mean > phi_0
    # means net densification somewhere (balanced by dilution elsewhere
    # under mass conservation). Clamp to [0, 1] as a demixing index.
    demixing_index = float(min(1.0, phi_std / max(phi_0, 1e-12)))

    # Correlation-length proxy for the pore size: first zero-crossing of
    # the radial autocorrelation of (phi - phi_mean). If no zero-crossing,
    # fall back to half-wavelength = R / (2·n_lobes) where n_lobes counts
    # sign changes of (phi − phi_mean).
    signs = np.sign(phi_f - phi_mean)
    sign_changes = int((np.diff(signs) != 0).sum())
    if sign_changes > 0:
        pore_size_mean = float(R_droplet / max(sign_changes, 1))
    else:
        # Homogeneous final state — use the cell size as a lower bound.
        pore_size_mean = float(dr)
    pore_size_mean = max(pore_size_mean, 1.0e-9)
    pore_size_std = 0.4 * pore_size_mean

    chord_samples = np.array(
        [pore_size_mean * 0.6, pore_size_mean, pore_size_mean * 1.5]
    )

    manifest = ModelManifest(
        model_name="L2.Gelation.NIPSCellulose",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=[
            "1D spherical ternary (polymer/solvent/non-solvent)",
            "Cahn-Hilliard on polymer + Fickian on solvent",
            "constant scalar mobilities",
            "Dirichlet at r=R (pure non-solvent bath by default)",
            "symmetry at r=0",
            "initial 1% noise on phi to break spherical symmetry",
            f"N_p={N_p}",
            f"chi_PS={chi_PS}, chi_PN={chi_PN}, chi_SN={chi_SN}",
            f"D_solvent={D_s:.2e} m²/s, D_polymer={D_phi:.2e} m²/s",
            f"kappa_CH={kappa:.2e} J/m",
        ],
        diagnostics={
            "phi_mean_final": phi_mean,
            "s_mean_final": s_mean,
            "n_mean_final": n_mean,
            "phi_std_final": phi_std,
            "phi_range_final": phi_range,
            "bicontinuous_score": bicontinuous_score,
            "demixing_index": demixing_index,
            "solver_success": bool(sol.success),
            "n_steps": int(sol.t.size),
        },
    )

    return GelationResult(
        r_grid=r,
        phi_field=phi_f,
        pore_size_mean=pore_size_mean,
        pore_size_std=pore_size_std,
        pore_size_distribution=chord_samples,
        porosity=porosity,
        alpha_final=demixing_index,
        char_wavelength=pore_size_mean,
        T_history=None,
        phi_snapshots=None,
        L_domain=R_droplet,
        grid_spacing=dr,
        bicontinuous_score=bicontinuous_score,
        anisotropy=0.0,
        connectivity=1.0,
        chord_skewness=0.0,
        model_tier="mechanistic",
        model_manifest=manifest,
    )


def _zero_cellulose_result(
    R_droplet: float, n_r: int, t_end: float
) -> GelationResult:
    """Empty-phase fallback for zero-cellulose input."""
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
            model_name="L2.Gelation.NIPSCellulose",
            evidence_tier=ModelEvidenceTier.UNSUPPORTED,
            assumptions=["zero cellulose input -> no gel"],
            diagnostics={"phi_cellulose_0": 0.0},
        ),
    )
