"""Node F1-a Phase 2a (v8.0): L2 diffusion-limited Ca²⁺ ionic gelation
solver for alginate microspheres.

Shrinking-core 1D spherical PDE for Ca²⁺ transport + "egg-box"
crosslinking on guluronate residues. State fields (all on a radial
grid 0 < r <= R_droplet):

    C(r, t)  — [mol/m³]  free Ca²⁺ concentration
    G(r, t)  — [mol/m³]  unbound guluronate concentration
    X(r, t)  — [mol/m³]  crosslink (egg-box junction) density

Governing equations:

    dC/dt = D·(1/r²)·d/dr(r²·dC/dr) − 2·k·C·G
    dG/dt = −k·C·G
    dX/dt = +½·k·C·G

(Two guluronate residues consumed per Ca²⁺ → one egg-box junction.)

Boundary conditions:
    C(R, t) = C_Ca_bath   (external CaCl₂ bath, Dirichlet)
    dC/dr(0, t) = 0       (spherical symmetry)
    G(r, 0) = f_G · c_alginate / M_alg_repeat
    X(r, 0) = 0
    C(r < R, 0) = 0

Spatial discretisation: finite volume on a uniform radial grid of
n_r ≈ 50 nodes; flux-limited to preserve mass.
Temporal: scipy BDF (stiff; D · t / R² is fast next to reaction).

Output populates a :class:`GelationResult` with the same schema as
other L2 solvers so downstream L3 / L4 code is platform-agnostic.

See `docs/f1a_alginate_protocol.md` §4.1 for the full derivation and
`docs/node32_cluster_f_v8_roadmap.md` §2.1 for the platform roadmap.
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


# Molar mass of an alginate repeat unit (mannuronate or guluronate,
# ~194 g/mol for the anhydrous Na-salt form). Used only when the
# caller passes c_alginate in kg/m³ and we need mol/m³ for the
# guluronate book-keeping.
M_ALG_REPEAT = 0.198  # [kg/mol]


def _guluronate_concentration_mol_m3(
    c_alginate_kgm3: float, f_G: float
) -> float:
    """Convert alginate mass concentration + guluronate fraction to
    guluronate molar concentration in mol/m³."""
    return (c_alginate_kgm3 / M_ALG_REPEAT) * f_G


def solve_ionic_ca_gelation(
    params: SimulationParameters,
    props: MaterialProperties,
    *,
    R_droplet: float,
    C_Ca_bath: float = 100.0,     # [mol/m³] 100 mM CaCl₂
    n_r: int = 50,
    time: Optional[float] = None,
    rtol: float = 1e-6,
    atol: float = 1e-12,
) -> GelationResult:
    """Solve shrinking-core Ca²⁺ ionic gelation for an alginate droplet.

    Parameters
    ----------
    params : SimulationParameters
        Uses ``formulation.c_agarose`` as the alginate concentration
        proxy (the legacy field is reused here until a dedicated
        ``c_alginate`` field lands; F1-a Phase 2b cleanup). Also uses
        ``formulation.t_crosslink`` as the gelation time budget.
    props : MaterialProperties
        Must have ``polymer_family == PolymerFamily.ALGINATE``. Reads
        ``f_guluronate``, ``D_Ca``, ``k_bind_Ca`` for the kinetic
        parameters.
    R_droplet : float
        Droplet radius [m].
    C_Ca_bath : float
        External Ca²⁺ bath concentration [mol/m³]. Default 100 mM.
    n_r : int
        Radial grid node count. 50 is a reasonable default for
        10–500 µm droplets.
    time : float, optional
        Gelation time [s]. Defaults to ``params.formulation.t_crosslink``.
    rtol, atol : float
        scipy.integrate.solve_ivp tolerances.

    Returns
    -------
    GelationResult
        ``pore_size_mean`` proxy = mean residual G fraction × R_droplet /
        n_r (crude estimate of uncrosslinked pocket scale);
        ``porosity`` computed from the polymer volume fraction minus
        the crosslink volume; ``model_manifest`` tagged
        SEMI_QUANTITATIVE.
    """
    if R_droplet <= 0:
        raise ValueError(f"R_droplet must be > 0, got {R_droplet}")
    if C_Ca_bath < 0:
        raise ValueError(f"C_Ca_bath must be >= 0, got {C_Ca_bath}")
    if n_r < 8:
        raise ValueError(f"n_r must be >= 8 for stability, got {n_r}")

    t_end = float(time) if time is not None else params.formulation.t_crosslink
    if t_end < 0:
        raise ValueError(f"time must be >= 0, got {t_end}")

    # Prefer the dedicated c_alginate field (Phase 2b); fall back to
    # c_agarose slot only when c_alginate is zero (Phase 2a compat).
    c_alg = float(
        params.formulation.c_alginate
        if params.formulation.c_alginate > 0.0
        else params.formulation.c_agarose
    )  # [kg/m³]
    f_G = float(props.f_guluronate)
    D = float(props.D_Ca)
    k = float(props.k_bind_Ca)

    G0 = _guluronate_concentration_mol_m3(c_alg, f_G)  # [mol/m³]
    if G0 <= 0:
        # No alginate / no guluronate -> no gel. Return a zero-gel
        # result rather than raising; downstream will detect via
        # porosity=1.0 and G_DN-related downstream computations will
        # land at SEMI_QUANTITATIVE → NO.
        return _zero_alginate_result(R_droplet, n_r, t_end)

    # Radial grid (node-centred, uniform)
    dr = R_droplet / n_r
    r = np.linspace(0.5 * dr, R_droplet - 0.5 * dr, n_r)

    # Initial state: C = 0 inside, G = G0 everywhere, X = 0
    y0 = np.concatenate([
        np.zeros(n_r),     # C
        np.full(n_r, G0),  # G
        np.zeros(n_r),     # X
    ])

    def rhs(_t, y):
        C = y[0:n_r]
        G = y[n_r:2 * n_r]
        # X term doesn't feed back into C or G dynamics; track via rate
        # dC/dt = D * laplacian_spherical(C) - 2*k*C*G
        # dG/dt = -k*C*G
        # dX/dt =  0.5 * k*C*G

        # Spherical Laplacian on node-centred grid.
        # (1/r²) d/dr (r² dC/dr) = via flux at faces
        # Face positions: r_face[i] = r[i] + dr/2 for i in 0..n_r-1
        # Face area ∝ r_face²
        r_face = np.linspace(dr, R_droplet, n_r)  # faces at dr, 2dr, ..., R
        # Gradient at each interior face
        # Outer BC: C(R) = C_Ca_bath. Face n_r-1 is the outer boundary.
        C_padded = np.concatenate([C, [C_Ca_bath]])
        grad = (C_padded[1:] - C_padded[:-1]) / dr      # grad at r_face[0..n_r-1]
        # grad at r=0 (centre) is 0 by symmetry — need a face at r=0 too
        # Incorporate: divergence = (A_out*flux_out - A_in*flux_in) / V
        # Inner face at each cell i is r_face[i-1] for i>=1, and r=0 for i=0
        # Let's compute per-cell divergence explicitly:
        rxn = -2.0 * k * C * G
        dC = np.zeros_like(C)
        # Cell volumes: V_i = (4/3)*pi*(r_out^3 - r_in^3), but ratios suffice
        for i in range(n_r):
            r_out = r_face[i]
            r_in = r_face[i - 1] if i > 0 else 0.0
            # Flux at outer face i: D * A_out * grad[i]; A ∝ r_out²
            flux_out = D * (r_out ** 2) * grad[i]
            # Flux at inner face
            if i == 0:
                flux_in = 0.0  # symmetry
            else:
                flux_in = D * (r_in ** 2) * grad[i - 1]
            V_shell = (r_out ** 3 - r_in ** 3) / 3.0  # omit 4*pi (cancels)
            dC[i] = (flux_out - flux_in) / V_shell + rxn[i]

        dG = -k * C * G
        dX = 0.5 * k * C * G
        return np.concatenate([dC, dG, dX])

    sol = solve_ivp(
        rhs, (0.0, t_end), y0,
        method="BDF", rtol=rtol, atol=atol,
    )

    if not sol.success:
        logger.warning("Ionic Ca gelation solver failed: %s", sol.message)

    y_final = sol.y[:, -1]
    C_f = y_final[0:n_r]
    G_f = y_final[n_r:2 * n_r]
    X_f = y_final[2 * n_r:3 * n_r]

    # Volume-weighted averages
    r_face = np.linspace(dr, R_droplet, n_r)
    r_in = np.concatenate([[0.0], r_face[:-1]])
    V_shells = (r_face ** 3 - r_in ** 3) / 3.0
    V_total = V_shells.sum()

    G_mean_final = float((G_f * V_shells).sum() / V_total)
    X_mean_final = float((X_f * V_shells).sum() / V_total)
    conversion_mean = float(1.0 - G_mean_final / max(G0, 1e-300))

    # Residual-G fraction profile is the "un-gelled pocket" field; its
    # characteristic length is related to the core-to-shell diffusion
    # lag. We use a simple proxy: mean pocket length = dr * residual
    # fraction_mean, scaled to µm. Crude but downstream consumers only
    # need a positive scalar; Phase 2b refines via Ogston-style model.
    fraction_ungelled = G_mean_final / max(G0, 1e-300)
    pore_size_mean = max(dr * max(fraction_ungelled, 0.01), 1.0e-9)
    pore_size_std = 0.3 * pore_size_mean

    # Porosity: polymer occupies ~(alginate mass / density_polymer).
    # Assume alginate dry density 1600 kg/m³; polymer volume fraction
    # phi_p = c_alg / 1600. Porosity = 1 - phi_p * hydration_factor.
    phi_p = c_alg / 1600.0
    porosity = max(0.1, min(0.98, 1.0 - phi_p * 3.0))

    # Characteristic "pore chord length" distribution — in absence of
    # a proper chord analysis, emit a 3-bin histogram centred on the
    # mean so downstream fields are populated.
    chord_samples = np.array([
        pore_size_mean * 0.7, pore_size_mean, pore_size_mean * 1.5,
    ])

    manifest = ModelManifest(
        model_name="L2.Gelation.IonicCaShrinkingCore",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=[
            "1D spherical symmetry",
            "Fickian Ca²⁺ diffusion",
            "second-order Ca²⁺ + 2·G binding (egg-box)",
            f"D_Ca={D:.2e} m²/s",
            f"k_bind={k:.2e} M⁻²·s⁻¹",
            f"f_guluronate={f_G}",
            f"C_Ca_bath={C_Ca_bath:.1f} mol/m³",
            "no competing counter-ion equilibria",
        ],
        diagnostics={
            "conversion_mean": conversion_mean,
            "X_mean_final": X_mean_final,
            "G_mean_final": G_mean_final,
            "solver_success": bool(sol.success),
            "n_steps": int(sol.t.size),
        },
    )

    return GelationResult(
        r_grid=r,
        phi_field=X_f / max(G0, 1e-300),  # crosslink fraction profile
        pore_size_mean=pore_size_mean,
        pore_size_std=pore_size_std,
        pore_size_distribution=chord_samples,
        porosity=porosity,
        alpha_final=conversion_mean,
        char_wavelength=pore_size_mean,
        T_history=None,
        phi_snapshots=None,
        L_domain=R_droplet,
        grid_spacing=dr,
        bicontinuous_score=0.0,  # alginate ionic gel is not bicontinuous
        anisotropy=0.0,
        connectivity=1.0,
        chord_skewness=0.0,
        model_tier="mechanistic",
        model_manifest=manifest,
    )


def solve_internal_gelation(
    params: SimulationParameters,
    props: MaterialProperties,
    *,
    R_droplet: float,
    C_CaCO3_0: float = 20.0,     # [mol/m³] initial CaCO₃ loading
    L_GDL_0: Optional[float] = None,  # default = 2 × C_CaCO3_0 (stoichiometric)
    k_hyd: float = 1.5e-4,        # [1/s] GDL hydrolysis (Draget 1997 25 °C)
    k_diss: float = 1.0e-2,       # [m³/(mol·s)] CaCO₃ + H⁺ pseudo-first-order
    n_r: int = 50,
    time: Optional[float] = None,
    rtol: float = 1e-6,
    atol: float = 1e-12,
) -> GelationResult:
    """Coupled GDL + CaCO₃ + Ca²⁺ + alginate solver for homogeneous
    *internal-release* ionic gelation.

    Unlike :func:`solve_ionic_ca_gelation`, which takes a Dirichlet
    Ca²⁺ bath at ``r = R`` (shrinking-core), this solver models the
    canonical GDL/CaCO₃ internal-release route (Draget 1997). CaCO₃
    powder is dispersed uniformly in the alginate phase; GDL hydrolyses
    slowly to gluconic acid, releasing H⁺ that dissolves CaCO₃ and
    liberates Ca²⁺ in situ.

    State variables (3 spatially-uniform scalars + 3 radial fields):

    .. code-block:: text

        L(t)      GDL                                  [mol/m³]
        A(t)      gluconic acid ≈ free H⁺              [mol/m³]
        S(t)      remaining CaCO₃                      [mol/m³]
        C(r,t)    Ca²⁺                                  [mol/m³]
        G(r,t)    free guluronate                       [mol/m³]
        X(r,t)    egg-box crosslink density             [mol/m³]

    Rate equations:

    .. code-block:: text

        dL/dt = −k_hyd · L                          (GDL → gluconic acid)
        dA/dt = +k_hyd · L − k_diss · S · A          (acid produced, consumed)
        dS/dt = −k_diss · S · A / 2                  (2 H⁺ per CaCO₃)
        dC/dt = D·∇²C + k_diss·S·A/2 − 2·k_bind·C·G  (diffusion + Ca source/sink)
        dG/dt = −k_bind · C · G
        dX/dt = +½ · k_bind · C · G

    Boundary conditions:

    - ``∂C/∂r(0) = 0`` (spherical symmetry)
    - ``∂C/∂r(R) = 0`` (**no-flux** — internal release, droplet sealed
      in oil phase; all Ca²⁺ from CaCO₃ stays inside).

    The scalar/field split exploits the homogeneity assumption: GDL
    and CaCO₃ particles are dispersed uniformly at ``t = 0`` and react
    slowly compared to mixing, so ``L``, ``A``, ``S`` stay
    approximately uniform. Ca²⁺ is produced uniformly too, but free
    Ca²⁺ then diffuses in response to binding gradients — hence ``C``,
    ``G``, ``X`` are radial fields. This is a Phase 3 refinement of
    the v8.0-rc1 lumped-parameter ``effective_bath_concentration``
    approximation; it preserves spatial information (useful for large
    beads where G gradients still matter) while honouring the
    GDL-limited kinetic time scale.

    Parameters
    ----------
    params : SimulationParameters
        Uses ``formulation.c_alginate`` (fall back to
        ``formulation.c_agarose`` if 0) and ``formulation.t_crosslink``.
    props : MaterialProperties
        Uses ``f_guluronate``, ``D_Ca``, ``k_bind_Ca``.
    R_droplet : float
        Droplet radius [m].
    C_CaCO3_0 : float
        Initial CaCO₃ loading [mol/m³]. Default 20 mM matches the
        ``gdl_caco3_internal`` reagent-library profile.
    L_GDL_0 : float, optional
        Initial GDL loading [mol/m³]. Defaults to ``2 × C_CaCO3_0``
        (stoichiometric 2:1 H⁺:Ca).
    k_hyd : float
        GDL hydrolysis rate [1/s]. 1.5e-4 s⁻¹ at 25 °C from Draget
        et al. 1997 *Int. J. Biol. Macromol.* 21:47.
    k_diss : float
        CaCO₃ dissolution rate constant [m³/(mol·s)], pseudo-first
        order in [H⁺]. 1e-2 matches small (<1 µm) dispersed particles
        per Plummer et al. 1978 / Pokrovsky & Schott 2002.

    Returns
    -------
    GelationResult
        Same schema as :func:`solve_ionic_ca_gelation`, with manifest
        tagged ``L2.Gelation.IonicCaInternalRelease`` and additional
        diagnostics (``L_final``, ``A_final``, ``S_final``).
    """
    if R_droplet <= 0:
        raise ValueError(f"R_droplet must be > 0, got {R_droplet}")
    if C_CaCO3_0 < 0:
        raise ValueError(f"C_CaCO3_0 must be >= 0, got {C_CaCO3_0}")
    if n_r < 8:
        raise ValueError(f"n_r must be >= 8 for stability, got {n_r}")

    t_end = float(time) if time is not None else params.formulation.t_crosslink
    if t_end < 0:
        raise ValueError(f"time must be >= 0, got {t_end}")

    if L_GDL_0 is None:
        L_GDL_0 = 2.0 * C_CaCO3_0  # stoichiometric default

    c_alg = float(
        params.formulation.c_alginate
        if params.formulation.c_alginate > 0.0
        else params.formulation.c_agarose
    )
    f_G = float(props.f_guluronate)
    D = float(props.D_Ca)
    k_bind = float(props.k_bind_Ca)

    G0 = _guluronate_concentration_mol_m3(c_alg, f_G)
    if G0 <= 0 or C_CaCO3_0 <= 0:
        return _zero_alginate_result(R_droplet, n_r, t_end)

    # Radial grid (node-centred, uniform)
    dr = R_droplet / n_r
    r = np.linspace(0.5 * dr, R_droplet - 0.5 * dr, n_r)
    r_face = np.linspace(dr, R_droplet, n_r)  # outer face of each cell
    r_in_arr = np.concatenate([[0.0], r_face[:-1]])
    V_shells = (r_face ** 3 - r_in_arr ** 3) / 3.0  # × 4π omitted (cancels)
    V_total = float(V_shells.sum())

    # State layout: [L, A, S, C(n_r), G(n_r), X(n_r)]
    n_state = 3 + 3 * n_r
    y0 = np.zeros(n_state)
    y0[0] = L_GDL_0
    y0[1] = 0.0
    y0[2] = C_CaCO3_0
    # C(r,0) = 0 already
    y0[3 + n_r: 3 + 2 * n_r] = G0
    # X(r,0) = 0 already

    def rhs(_t, y):
        L = max(float(y[0]), 0.0)
        A = max(float(y[1]), 0.0)
        S = max(float(y[2]), 0.0)
        C = np.maximum(y[3: 3 + n_r], 0.0)
        G = np.maximum(y[3 + n_r: 3 + 2 * n_r], 0.0)

        # Scalar ODEs
        rate_hyd = k_hyd * L                       # GDL → acid [mol/(m³·s)]
        rate_H_loss = k_diss * S * A               # H⁺ consumed [mol/(m³·s)]
        rate_Ca_src = 0.5 * rate_H_loss            # Ca²⁺ produced [mol/(m³·s)]
        dL = -rate_hyd
        dA = +rate_hyd - rate_H_loss
        dS = -rate_Ca_src

        # Spherical Laplacian on node-centred grid, no-flux at both ends
        # Interior face gradients (n_r - 1 of them)
        grad_interior = (C[1:] - C[:-1]) / dr      # at r_face[0..n_r-2]

        dC = np.zeros(n_r)
        rxn_bind = -2.0 * k_bind * C * G
        for i in range(n_r):
            # outer face of cell i
            if i == n_r - 1:
                flux_out = 0.0                      # no-flux BC at r=R
            else:
                flux_out = D * (r_face[i] ** 2) * grad_interior[i]
            # inner face of cell i
            if i == 0:
                flux_in = 0.0                       # symmetry BC at r=0
            else:
                flux_in = D * (r_face[i - 1] ** 2) * grad_interior[i - 1]
            dC[i] = (flux_out - flux_in) / V_shells[i] + rate_Ca_src + rxn_bind[i]

        dG = -k_bind * C * G
        dX = 0.5 * k_bind * C * G

        out = np.empty_like(y)
        out[0] = dL
        out[1] = dA
        out[2] = dS
        out[3: 3 + n_r] = dC
        out[3 + n_r: 3 + 2 * n_r] = dG
        out[3 + 2 * n_r:] = dX
        return out

    sol = solve_ivp(
        rhs, (0.0, t_end), y0,
        method="BDF", rtol=rtol, atol=atol,
    )

    if not sol.success:
        logger.warning("Internal-release gelation solver failed: %s", sol.message)

    y_final = sol.y[:, -1]
    L_f = float(max(y_final[0], 0.0))
    A_f = float(max(y_final[1], 0.0))
    S_f = float(max(y_final[2], 0.0))
    C_f = np.maximum(y_final[3: 3 + n_r], 0.0)
    G_f = np.maximum(y_final[3 + n_r: 3 + 2 * n_r], 0.0)
    X_f = np.maximum(y_final[3 + 2 * n_r:], 0.0)

    G_mean_final = float((G_f * V_shells).sum() / V_total)
    X_mean_final = float((X_f * V_shells).sum() / V_total)
    conversion_mean = float(1.0 - G_mean_final / max(G0, 1e-300))

    # Homogeneity metric: coefficient of variation of X across r
    X_std = float(np.std(X_f)) if X_f.size > 0 else 0.0
    X_cov = X_std / max(X_mean_final, 1e-300)

    # Pore / porosity mirror the shrinking-core solver's estimates
    fraction_ungelled = G_mean_final / max(G0, 1e-300)
    pore_size_mean = max(dr * max(fraction_ungelled, 0.01), 1.0e-9)
    pore_size_std = 0.3 * pore_size_mean
    phi_p = c_alg / 1600.0
    porosity = max(0.1, min(0.98, 1.0 - phi_p * 3.0))
    chord_samples = np.array([
        pore_size_mean * 0.7, pore_size_mean, pore_size_mean * 1.5,
    ])

    manifest = ModelManifest(
        model_name="L2.Gelation.IonicCaInternalRelease",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=[
            "1D spherical with no-flux outer boundary",
            "GDL, gluconic acid, CaCO₃ spatially uniform",
            "Fickian Ca²⁺ diffusion, second-order Ca²⁺+2G binding",
            f"k_hyd={k_hyd:.2e} /s (Draget 1997 25 °C)",
            f"k_diss={k_diss:.2e} m³/(mol·s) (pseudo-first-order in H⁺)",
            f"D_Ca={D:.2e} m²/s",
            f"k_bind={k_bind:.2e} M⁻²·s⁻¹",
            f"f_guluronate={f_G}",
            f"C_CaCO3_0={C_CaCO3_0:.2f} mol/m³",
            f"L_GDL_0={L_GDL_0:.2f} mol/m³",
            "full dissociation of gluconic acid (pKa 3.86) → H⁺ = [acid]",
            "stoichiometry 2 H⁺ : 1 Ca²⁺ : 1 CaCO₃",
        ],
        diagnostics={
            "conversion_mean": conversion_mean,
            "X_mean_final": X_mean_final,
            "X_cov": X_cov,               # homogeneity metric
            "G_mean_final": G_mean_final,
            "L_final": L_f,
            "A_final": A_f,
            "S_final": S_f,
            "solver_success": bool(sol.success),
            "n_steps": int(sol.t.size),
        },
    )

    return GelationResult(
        r_grid=r,
        phi_field=X_f / max(G0, 1e-300),
        pore_size_mean=pore_size_mean,
        pore_size_std=pore_size_std,
        pore_size_distribution=chord_samples,
        porosity=porosity,
        alpha_final=conversion_mean,
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


def _zero_alginate_result(
    R_droplet: float, n_r: int, t_end: float
) -> GelationResult:
    """Empty-gel fallback for zero-alginate inputs."""
    r = np.linspace(0.0, R_droplet, n_r)
    return GelationResult(
        r_grid=r,
        phi_field=np.zeros(n_r),
        pore_size_mean=R_droplet,  # "all pore"
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
            model_name="L2.Gelation.IonicCaShrinkingCore",
            evidence_tier=ModelEvidenceTier.UNSUPPORTED,
            assumptions=["zero alginate input -> no gel"],
            diagnostics={"G0": 0.0},
        ),
    )
