"""Chromatography orchestrators: breakthrough and gradient elution.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase C + E.

Phase C: run_breakthrough — single-component Langmuir LRM breakthrough.
Phase E: run_gradient_elution — multi-component gradient chromatography with
         competitive Langmuir, process metrics (yield, purity, resolution).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import numpy as np
from scipy.integrate import solve_ivp

from ..datatypes import ModelEvidenceTier, ModelManifest
from .hydrodynamics import ColumnGeometry
from .isotherms.langmuir import LangmuirIsotherm
from .isotherms.competitive_langmuir import CompetitiveLangmuirIsotherm
from .transport.lumped_rate import solve_lrm
from .detection.uv import compute_uv_signal, apply_detector_broadening
from .gradient import GradientProgram

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


# ─── Phase C: Breakthrough Result ────────────────────────────────────────────

@dataclass
class BreakthroughResult:
    """Complete breakthrough simulation result.

    Attributes:
        time: Time array [s].
        uv_signal: UV absorbance signal [mAU].
        C_outlet: Outlet concentration [mol/m^3].
        dbc_5pct: Dynamic binding capacity at 5% breakthrough [mol/m^3 column].
        dbc_10pct: Dynamic binding capacity at 10% breakthrough [mol/m^3 column].
        dbc_50pct: Dynamic binding capacity at 50% breakthrough [mol/m^3 column].
        pressure_drop: Column pressure drop [Pa].
        mass_balance_error: Relative mass balance error [-].
    """

    time: np.ndarray
    uv_signal: np.ndarray
    C_outlet: np.ndarray
    dbc_5pct: float
    dbc_10pct: float
    dbc_50pct: float
    pressure_drop: float
    mass_balance_error: float
    # v6.1 (Node 5): evidence provenance — populated by run_breakthrough.
    # Inherits the weakest tier from the upstream FMC manifest, then applies
    # M3-specific gates (mass-balance > 5% downgrades to QUALITATIVE_TREND).
    model_manifest: Optional[ModelManifest] = None


def _compute_dbc(
    time: np.ndarray,
    C_outlet: np.ndarray,
    C_feed: float,
    flow_rate: float,
    column_volume: float,
    threshold_fraction: float,
) -> float:
    """Compute dynamic binding capacity at a given breakthrough threshold.

    DBC is defined as the mass of protein loaded per unit column volume
    at the time when the outlet concentration first reaches
    threshold_fraction * C_feed.

    Args:
        time: Time array [s].
        C_outlet: Outlet concentration [mol/m^3].
        C_feed: Feed concentration [mol/m^3].
        flow_rate: Volumetric flow rate [m^3/s].
        column_volume: Total bed volume [m^3].
        threshold_fraction: Breakthrough threshold (e.g., 0.05 for 5%).

    Returns:
        DBC [mol/m^3] (of column volume).
        Returns NaN if breakthrough is never reached.
    """
    C_threshold = threshold_fraction * C_feed

    # Find first time C_outlet >= C_threshold
    above = np.where(C_outlet >= C_threshold)[0]
    if len(above) == 0:
        # Breakthrough not reached — DBC is the total loaded
        t_bt = time[-1]
    else:
        idx = above[0]
        # Linear interpolation for more precise breakthrough time
        if idx > 0:
            C_lo, C_hi = C_outlet[idx - 1], C_outlet[idx]
            t_lo, t_hi = time[idx - 1], time[idx]
            if C_hi > C_lo:
                frac = (C_threshold - C_lo) / (C_hi - C_lo)
                t_bt = t_lo + frac * (t_hi - t_lo)
            else:
                t_bt = time[idx]
        else:
            t_bt = time[0]

    # Mass loaded up to breakthrough time (feed minus what passed through)
    mask = time <= t_bt
    if np.sum(mask) < 2:
        return 0.0

    t_bt_arr = time[mask]
    C_out_bt = C_outlet[mask]

    # Net mass captured = Q * integral(C_feed - C_outlet) dt
    mass_captured = flow_rate * float(np.trapezoid(C_feed - C_out_bt, t_bt_arr))
    mass_captured = max(mass_captured, 0.0)

    # DBC = mass_captured / column_volume [mol/m^3]
    return mass_captured / column_volume


def run_breakthrough(
    column: ColumnGeometry,
    microsphere=None,
    C_feed: float = 1.0,
    flow_rate: float = 1e-8,
    feed_duration: float = 600.0,
    total_time: float = 1200.0,
    extinction_coeff: float = 36000.0,
    sigma_detector: float = 1.0,
    isotherm: LangmuirIsotherm | None = None,
    n_z: int = 50,
    D_molecular: float = 7e-11,
    k_ads: float = 100.0,
    fmc=None,
    process_state: dict | None = None,
) -> BreakthroughResult:
    """Run a single-component breakthrough simulation.

    If a FunctionalMicrosphere is provided, its mechanical properties
    (G_DN, E_star) and particle geometry override column defaults.

    Args:
        column: Column geometry and hydraulic parameters.
        microsphere: Optional FunctionalMicrosphere from Module 2.
        C_feed: Feed concentration [mol/m^3].
        flow_rate: Volumetric flow rate [m^3/s].
        feed_duration: Duration of feed/loading step [s].
        total_time: Total simulation time [s].
        extinction_coeff: Molar extinction coefficient [1/(M*cm)].
        sigma_detector: Detector broadening sigma [s].
        isotherm: Langmuir isotherm (default parameters if None).
        n_z: Number of axial finite-volume cells.
        D_molecular: Molecular diffusivity [m^2/s].
        k_ads: Adsorption rate constant [1/s].

    Returns:
        BreakthroughResult with DBC values, UV signal, and mass balance.
    """
    # ── Override column with microsphere properties if available ──
    if microsphere is not None:
        m1 = microsphere.m1_contract
        column = ColumnGeometry(
            diameter=column.diameter,
            bed_height=column.bed_height,
            particle_diameter=m1.bead_d50,
            bed_porosity=column.bed_porosity,
            particle_porosity=m1.porosity,
            G_DN=microsphere.G_DN_updated or m1.G_DN,
            E_star=microsphere.E_star_updated or m1.E_star,
        )

    # ── Validate flow rate ──
    flow_warnings = column.validate_flow_rate(flow_rate)
    for w in flow_warnings:
        logger.warning(w)

    # ── Isotherm selection (v5.9.5 H3: auto-route from FMC) ──
    if isotherm is None and fmc is not None:
        from .isotherms.adapter import select_isotherm_from_fmc
        isotherm = select_isotherm_from_fmc(fmc, process_state)
        logger.info("Auto-selected isotherm from FMC: %s", type(isotherm).__name__)
    elif isotherm is None:
        isotherm = LangmuirIsotherm()

    # ── Pressure drop ──
    pressure_drop = column.pressure_drop(flow_rate)

    # ── Solve LRM ──
    lrm_result = solve_lrm(
        column=column,
        isotherm=isotherm,
        C_feed=C_feed,
        feed_duration=feed_duration,
        flow_rate=flow_rate,
        total_time=total_time,
        n_z=n_z,
        D_molecular=D_molecular,
        k_ads=k_ads,
    )

    # ── UV detection ──
    uv_raw = compute_uv_signal(
        lrm_result.C_outlet,
        extinction_coeff=extinction_coeff,
    )
    uv_signal = apply_detector_broadening(uv_raw, lrm_result.time, sigma_detector)

    # ── Dynamic binding capacity ──
    V_col = column.bed_volume
    dbc_5 = _compute_dbc(
        lrm_result.time, lrm_result.C_outlet, C_feed, flow_rate, V_col, 0.05
    )
    dbc_10 = _compute_dbc(
        lrm_result.time, lrm_result.C_outlet, C_feed, flow_rate, V_col, 0.10
    )
    dbc_50 = _compute_dbc(
        lrm_result.time, lrm_result.C_outlet, C_feed, flow_rate, V_col, 0.50
    )

    manifest = _build_m3_chrom_manifest(
        model_basename="M3.breakthrough.LRM",
        isotherm=isotherm,
        fmc=fmc,
        worst_mass_balance_error=lrm_result.mass_balance_error,
        diagnostics_extra={
            "dbc_5pct": float(dbc_5),
            "dbc_10pct": float(dbc_10),
            "dbc_50pct": float(dbc_50),
            "pressure_drop_Pa": float(pressure_drop),
        },
    )

    return BreakthroughResult(
        time=lrm_result.time,
        uv_signal=uv_signal,
        C_outlet=lrm_result.C_outlet,
        dbc_5pct=dbc_5,
        dbc_10pct=dbc_10,
        dbc_50pct=dbc_50,
        pressure_drop=pressure_drop,
        mass_balance_error=lrm_result.mass_balance_error,
        model_manifest=manifest,
    )


# ─── Node 5: M3 manifest builder (shared by breakthrough + gradient) ────────


# Mass-balance gates (consensus plan §5):
#   <= 2%  -> no downgrade
#   2-5%  -> caution flag in diagnostics
#   > 5%  -> tier capped at QUALITATIVE_TREND (output not decision-grade)
_M3_MB_CAUTION = 0.02
_M3_MB_BLOCKER = 0.05


def _build_m3_chrom_manifest(
    model_basename: str,
    isotherm,
    fmc,
    worst_mass_balance_error: float,
    diagnostics_extra: dict | None = None,
) -> ModelManifest:
    """Build a chromatography-result ModelManifest.

    Tier inheritance:
      * Start from the upstream FMC manifest tier when available (M3 cannot
        be stronger than its inputs). Fall back to SEMI_QUANTITATIVE when
        the caller passed no FMC (CLI breakthrough with default isotherm).
      * Apply the mass-balance gate: a balance error worse than 5% caps the
        tier at QUALITATIVE_TREND because DBC/yield/purity numbers cannot be
        defended quantitatively when the solver has lost or fabricated mass.

    Diagnostics carry the mass balance, isotherm class, and any extras
    supplied by the caller (DBC values, peak resolution, etc.) so RunReport
    consumers can roll the M3 evidence into a pipeline-wide tier.
    """
    _ORDER = list(ModelEvidenceTier)

    # Inherit upstream tier from FMC manifest, else default to SEMI.
    fmc_manifest = getattr(fmc, "model_manifest", None) if fmc is not None else None
    if fmc_manifest is not None:
        tier = fmc_manifest.evidence_tier
        upstream_assumptions = list(fmc_manifest.assumptions)
        upstream_calref = fmc_manifest.calibration_ref
    else:
        tier = ModelEvidenceTier.SEMI_QUANTITATIVE
        upstream_assumptions = []
        upstream_calref = ""

    # Mass-balance gate.
    mb_status = "ok"
    if worst_mass_balance_error > _M3_MB_BLOCKER:
        mb_status = "blocker"
        # Cap tier at QUALITATIVE_TREND: take the weaker of (current tier,
        # QUALITATIVE_TREND) using the canonical ordering where larger
        # _ORDER index == weaker tier.
        capped_idx = max(
            _ORDER.index(tier),
            _ORDER.index(ModelEvidenceTier.QUALITATIVE_TREND),
        )
        tier = _ORDER[capped_idx]
    elif worst_mass_balance_error > _M3_MB_CAUTION:
        mb_status = "caution"

    diagnostics: dict = {
        "isotherm_class": type(isotherm).__name__ if isotherm is not None else "None",
        "mass_balance_error": float(worst_mass_balance_error),
        "mass_balance_status": mb_status,
        "fmc_provided": fmc is not None,
    }
    if diagnostics_extra:
        diagnostics.update(diagnostics_extra)

    assumptions = upstream_assumptions + [
        "M3 inherits the weakest tier from the upstream FunctionalMediaContract; "
        "isotherm parameters are calibrated only when the FMC manifest is calibrated.",
    ]
    if mb_status == "blocker":
        assumptions.append(
            f"Mass balance error {worst_mass_balance_error:.1%} exceeds {_M3_MB_BLOCKER:.0%}; "
            "DBC, yield, and purity outputs are not decision-grade."
        )

    return ModelManifest(
        model_name=model_basename,
        evidence_tier=tier,
        valid_domain={},
        calibration_ref=upstream_calref,
        assumptions=assumptions,
        diagnostics=diagnostics,
    )


# ─── Phase E: Gradient Elution ───────────────────────────────────────────────

@dataclass
class PeakInfo:
    """Chromatographic peak statistics for one component.

    Attributes:
        component_idx: Component index [-].
        retention_time: Peak apex time [s].
        peak_height: Peak apex concentration [mol/m^3].
        peak_area: Integrated peak area [mol*s/m^3].
        peak_width_half: Peak width at half height [s].
        yield_fraction: Fraction of injected mass recovered in peak [-].
        purity: Fraction of peak area from this component [-].
    """

    component_idx: int
    retention_time: float
    peak_height: float
    peak_area: float
    peak_width_half: float
    yield_fraction: float
    purity: float


@dataclass
class GradientElutionResult:
    """Result of a multi-component gradient elution simulation.

    Attributes:
        time: Time array [s], shape (N_t,).
        C_outlet: Outlet concentration per component [mol/m^3], shape (n_comp, N_t).
        gradient_profile: Gradient value (salt / pH) vs time, shape (N_t,).
        uv_signal: Summed UV absorbance [mAU], shape (N_t,).
        peaks: Peak table — one PeakInfo per component.
        resolution: Peak resolution between consecutive component pairs, shape (n_comp-1,).
        overall_yield: Overall mass recovery fraction [-].
        pressure_drop: Column pressure drop [Pa].
        mass_balance_errors: Mass balance error per component [-], shape (n_comp,).
    """

    time: np.ndarray
    C_outlet: np.ndarray                # shape (n_comp, N_t)
    gradient_profile: np.ndarray        # shape (N_t,)
    uv_signal: np.ndarray               # shape (N_t,)
    peaks: list[PeakInfo]
    resolution: np.ndarray              # shape (n_comp - 1,) or empty
    overall_yield: float
    pressure_drop: float
    mass_balance_errors: np.ndarray     # shape (n_comp,)
    gradient_affects_binding: bool = False  # True only for SMA/IMAC isotherms
    # v6.1 (Node 5): evidence provenance — populated by run_gradient_elution.
    # The mass-balance gate uses the worst per-component error.
    model_manifest: Optional[ModelManifest] = None


def _find_peak(
    time: np.ndarray,
    C_outlet_component: np.ndarray,
) -> tuple[float, float, float, float]:
    """Find peak apex and width for a single component outlet profile.

    Args:
        time: Time array [s].
        C_outlet_component: Outlet concentration [mol/m^3] for one component.

    Returns:
        Tuple of (retention_time, peak_height, peak_area, peak_width_half).
        Returns (nan, 0, 0, nan) if no discernible peak exists.
    """
    C = np.maximum(C_outlet_component, 0.0)

    if C.max() < 1e-20:
        return float("nan"), 0.0, 0.0, float("nan")

    # Apex
    apex_idx = int(np.argmax(C))
    t_apex = float(time[apex_idx])
    h_apex = float(C[apex_idx])

    # Area by trapezoidal integration
    area = float(np.trapezoid(C, time))

    # Width at half-height
    half_h = h_apex / 2.0
    above_half = C >= half_h

    # Find left and right edges of the half-height region
    if above_half.sum() < 2:
        width = 0.0
    else:
        indices = np.where(above_half)[0]
        i_left = indices[0]
        i_right = indices[-1]

        # Linear interpolation at left edge
        if i_left > 0 and C[i_left - 1] < half_h:
            frac_l = (half_h - C[i_left - 1]) / (C[i_left] - C[i_left - 1])
            t_left = time[i_left - 1] + frac_l * (time[i_left] - time[i_left - 1])
        else:
            t_left = time[i_left]

        # Linear interpolation at right edge
        if i_right < len(C) - 1 and C[i_right + 1] < half_h:
            frac_r = (half_h - C[i_right]) / (C[i_right + 1] - C[i_right])
            t_right = time[i_right] + frac_r * (time[i_right + 1] - time[i_right])
        else:
            t_right = time[i_right]

        width = float(t_right - t_left)

    return t_apex, h_apex, area, width


def _compute_resolution(
    t1: float, w1: float,
    t2: float, w2: float,
) -> float:
    """Compute chromatographic resolution between two peaks.

    Rs = 2 * |t2 - t1| / (w1 + w2)

    where w1, w2 are peak widths at half-height (FWHM).

    Args:
        t1, t2: Peak retention times [s].
        w1, w2: Peak widths at half height [s] (FWHM).

    Returns:
        Resolution Rs [-].  Returns 0 if widths are undefined.
    """
    if w1 + w2 <= 0 or np.isnan(t1) or np.isnan(t2):
        return 0.0
    return 2.0 * abs(t2 - t1) / (w1 + w2)


def _build_gradient_lrm_rhs(
    n_z: int,
    n_comp: int,
    dz: float,
    u: float,
    D_ax: np.ndarray,
    eps_b: float,
    eps_p: float,
    R_p: float,
    k_f: np.ndarray,
    k_ads: np.ndarray,
    isotherm: CompetitiveLangmuirIsotherm,
    C_feed: np.ndarray,
    gradient: GradientProgram,
    total_time: float,
    feed_duration: float,
    equilibrium_adapter=None,
    gradient_field: str = "salt_concentration",
):
    """Build the ODE RHS for multi-component gradient LRM.

    State vector layout (length = 3 * n_comp * n_z):
        y[0         : n_comp*n_z]          = C    (bulk, n_comp x n_z, row-major)
        y[n_comp*n_z: 2*n_comp*n_z]        = Cp   (pore, n_comp x n_z)
        y[2*n_comp*n_z: 3*n_comp*n_z]      = q    (bound, n_comp x n_z)

    v6.0 H6: When equilibrium_adapter is provided, it is updated with the
    gradient value at each time step to enable time-varying equilibrium.

    Args:
        feed_duration: Duration of protein loading phase [s].  After this time
            the inlet protein concentration is set to zero (elution/wash phase).
        equilibrium_adapter: Optional EquilibriumAdapter wrapping a gradient-sensitive
            isotherm. When provided, used instead of isotherm for equilibrium_loading().
        gradient_field: ProcessState field to update from gradient value.
    """
    nc = n_comp
    nz = n_z
    block = nc * nz

    # Precompute mass transfer coefficients
    mt_coeff = (3.0 / R_p) * k_f  # shape (n_comp,)

    # v6.0 H6: Determine equilibrium dispatch mode
    _use_adapter = equilibrium_adapter is not None

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        # Reshape state blocks: each (nc, nz)
        C  = y[:block].reshape(nc, nz).copy()
        Cp = y[block:2*block].reshape(nc, nz).copy()
        q  = y[2*block:3*block].reshape(nc, nz).copy()

        np.maximum(C,  0.0, out=C)
        np.maximum(Cp, 0.0, out=Cp)
        np.maximum(q,  0.0, out=q)

        # Gradient value at current time (salt, pH, imidazole, sugar, etc.)
        grad_val = gradient.value_at_time(t)

        # Feed phase switching: protein is present only during the load phase.
        # After feed_duration the inlet switches to protein-free buffer so that
        # bound protein is eluted by the rising gradient rather than being
        # continuously replaced by fresh feed.
        if t <= feed_duration:
            C_in = C_feed.copy()   # loading phase: protein in feed
        else:
            C_in = np.zeros(nc)    # elution/wash phase: no protein

        # Isotherm: q_eq for all cells, shape (nc, nz)
        # v6.0 H6: When adapter is present, update process_state with gradient
        # and use adapter for equilibrium dispatch (SMA, HIC, IMAC, lectin, pH).
        if _use_adapter:
            equilibrium_adapter._state[gradient_field] = grad_val
            q_eq = equilibrium_adapter.equilibrium_loading(Cp)
        else:
            q_eq = isotherm.equilibrium_loading(Cp)  # (nc, nz)

        dCdt  = np.zeros((nc, nz))
        dCpdt = np.zeros((nc, nz))
        dqdt  = np.zeros((nc, nz))

        for i in range(nc):
            # Advection: first-order upwind
            dCdz = np.empty(nz)
            dCdz[0]  = (C[i, 0] - C_in[i]) / dz
            dCdz[1:] = (C[i, 1:] - C[i, :-1]) / dz

            # Dispersion: central differences
            d2Cdz2 = np.empty(nz)
            d2Cdz2[0]    = (C[i, 1] - 2.0 * C[i, 0] + C_in[i]) / dz**2
            d2Cdz2[1:-1] = (C[i, 2:] - 2.0 * C[i, 1:-1] + C[i, :-2]) / dz**2
            d2Cdz2[-1]   = (C[i, -2] - C[i, -1]) / dz**2

            # Adsorption kinetics
            dqdt_i = k_ads[i] * (q_eq[i] - q[i])

            # Film mass transfer
            film = mt_coeff[i] * (C[i] - Cp[i])

            # PDEs
            dCdt[i]  = (-u * dCdz + D_ax[i] * d2Cdz2 - (1.0 - eps_b) * film) / eps_b
            dCpdt[i] = (film - (1.0 - eps_p) * dqdt_i) / eps_p
            dqdt[i]  = dqdt_i

        return np.concatenate([dCdt.ravel(), dCpdt.ravel(), dqdt.ravel()])

    return rhs


def run_gradient_elution(
    column: ColumnGeometry,
    C_feed: np.ndarray,
    gradient: GradientProgram,
    flow_rate: float,
    total_time: float,
    feed_duration: float | None = None,
    isotherm: CompetitiveLangmuirIsotherm | None = None,
    n_z: int = 50,
    D_molecular: float | np.ndarray = 7e-11,
    k_ads: float | np.ndarray = 100.0,
    extinction_coeffs: np.ndarray | None = None,
    sigma_detector: float = 1.0,
    mu: float = 1e-3,
    rho: float = 1000.0,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    process_state=None,
    gradient_field: str | None = None,
) -> GradientElutionResult:
    """Simulate multi-component gradient elution chromatography.

    Uses the Lumped Rate Model (LRM) with competitive Langmuir equilibrium
    and method-of-lines + BDF time integration.

    Args:
        column: Packed column geometry.
        C_feed: Feed concentrations per component [mol/m^3], shape (n_comp,).
        gradient: GradientProgram controlling the elution gradient
            (e.g., salt or pH profile).
        flow_rate: Volumetric flow rate [m^3/s].
        total_time: Total simulation time [s].
        feed_duration: Duration of the protein loading phase [s].  After this
            time the inlet is switched to protein-free buffer (elution phase).
            Defaults to the start time of the first gradient segment, i.e.
            ``gradient.segments[0][0]`` if > 0, otherwise ``total_time / 3``.
        isotherm: CompetitiveLangmuirIsotherm. If None, a 2-component
            default is used.
        n_z: Number of axial finite-volume cells.
        D_molecular: Molecular diffusivity [m^2/s]. Scalar or per-component
            array of shape (n_comp,).
        k_ads: Adsorption rate constant [1/s]. Scalar or per-component.
        extinction_coeffs: Molar extinction coefficients for UV [1/(M*cm)].
            Shape (n_comp,). Default 36000 for all components.
        sigma_detector: UV detector broadening sigma [s].
        mu: Dynamic viscosity [Pa.s].
        rho: Fluid density [kg/m^3].
        rtol: ODE solver relative tolerance.
        atol: ODE solver absolute tolerance.
        process_state: Optional ProcessState or dict with initial conditions
            for gradient-sensitive isotherms (SMA, HIC, IMAC, ProteinA, lectin).
            When provided with a gradient-sensitive isotherm, an EquilibriumAdapter
            is created and updated at each time step.
        gradient_field: ProcessState field to update from gradient value.
            Auto-detected from isotherm.gradient_field if not specified.
            Examples: "salt_concentration" (IEX/HIC), "imidazole" (IMAC),
            "ph" (Protein A), "sugar_competitor" (lectin).

    Returns:
        GradientElutionResult with outlet profiles, peak table, metrics, and
        ``gradient_affects_binding`` flag indicating whether the isotherm
        equilibrium depends on the gradient value.

    Raises:
        RuntimeError: If the ODE solver fails.
    """
    from .transport.lumped_rate import _wilson_geankoplis_kf, _axial_dispersion

    C_feed = np.asarray(C_feed, dtype=float)
    n_comp = len(C_feed)

    if isotherm is None:
        isotherm = CompetitiveLangmuirIsotherm(
            q_max=np.full(n_comp, 100.0),
            K_L=np.full(n_comp, 1e3),
        )

    # Determine feed_duration: protein inlet is active only during load phase.
    # Default: use the gradient's first segment start time if it is > 0
    # (i.e. the gradient only begins after an isocratic load step), otherwise
    # fall back to one-third of the total simulation time.
    if feed_duration is None:
        first_segment_time = float(gradient.segments[0][0])
        if first_segment_time > 0.0:
            feed_duration = first_segment_time
        else:
            feed_duration = total_time / 3.0

    # Detect whether the isotherm equilibrium is sensitive to the gradient
    # (e.g. SMA uses salt, IMAC uses imidazole, HIC uses salt, ProteinA uses pH).
    # Competitive Langmuir is not gradient-sensitive.
    gradient_affects_binding: bool = bool(
        getattr(isotherm, "gradient_sensitive", False)
    )

    # v6.0 H6: Build EquilibriumAdapter when isotherm is gradient-sensitive
    _adapter = None
    _grad_field = gradient_field
    if gradient_affects_binding:
        from .isotherms.adapter import EquilibriumAdapter
        _adapter = EquilibriumAdapter(isotherm, process_state)
        # Auto-detect gradient_field from isotherm if not explicitly provided
        if _grad_field is None:
            _grad_field = getattr(isotherm, "gradient_field", "salt_concentration")
        logger.info(
            "Gradient-aware elution: %s, field=%s",
            type(isotherm).__name__, _grad_field,
        )
    else:
        if _grad_field is None:
            _grad_field = "salt_concentration"

    if isotherm.n_components != n_comp:
        raise ValueError(
            f"isotherm has {isotherm.n_components} components but "
            f"C_feed has {n_comp}."
        )

    # ── Derived quantities ──
    u = column.superficial_velocity(flow_rate)
    eps_b = column.bed_porosity
    eps_p = column.particle_porosity
    R_p = column.particle_radius
    dp = column.particle_diameter
    L = column.bed_height

    # Per-component arrays
    D_mol = np.broadcast_to(np.asarray(D_molecular, dtype=float), (n_comp,)).copy()
    k_ads_arr = np.broadcast_to(np.asarray(k_ads, dtype=float), (n_comp,)).copy()

    # Mass transfer and dispersion coefficients
    k_f = np.array([
        _wilson_geankoplis_kf(u, dp, eps_b, D_mol[i], mu, rho)
        for i in range(n_comp)
    ])
    D_ax = np.array([
        _axial_dispersion(u, dp, eps_b)
        for _ in range(n_comp)
    ])

    # Grid
    dz = L / n_z

    # ── Build RHS (v6.0 H6: gradient-aware when adapter available) ──
    rhs = _build_gradient_lrm_rhs(
        n_z=n_z,
        n_comp=n_comp,
        dz=dz,
        u=u,
        D_ax=D_ax,
        eps_b=eps_b,
        eps_p=eps_p,
        R_p=R_p,
        k_f=k_f,
        k_ads=k_ads_arr,
        isotherm=isotherm,
        C_feed=C_feed,
        gradient=gradient,
        total_time=total_time,
        feed_duration=feed_duration,
        equilibrium_adapter=_adapter,
        gradient_field=_grad_field,
    )

    # ── Initial conditions: column empty ──
    block = n_comp * n_z
    y0 = np.zeros(3 * block)

    # ── Time evaluation points ──
    n_eval = max(300, int(total_time / 0.5))
    t_eval = np.linspace(0.0, total_time, n_eval)

    logger.info(
        "Gradient LRM: %d components, n_z=%d, u=%.3e m/s, total_time=%.0f s",
        n_comp, n_z, u, total_time,
    )

    # ── Solve ──
    # BDF stays here: this is the gradient-elution path, where the binding
    # equilibrium changes with time. LSODA gets stuck oscillating between
    # stiff/non-stiff modes when the Jacobian shifts each step (v9.1.1 #2).
    sol = solve_ivp(
        rhs,
        t_span=(0.0, total_time),
        y0=y0,
        method="BDF",
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        max_step=total_time / 50.0,
    )

    if not sol.success:
        raise RuntimeError(f"Gradient LRM solver failed: {sol.message}")

    time = sol.t
    N_t = len(time)

    # ── Extract outlet profiles ──
    y_out = sol.y  # shape (3*block, N_t)
    C_outlet = np.zeros((n_comp, N_t))
    q_avg = np.zeros((n_comp, N_t))

    for i in range(n_comp):
        # Bulk concentration in all cells for component i: rows i*n_z..(i+1)*n_z-1
        C_all_i = y_out[i * n_z:(i + 1) * n_z, :]  # (n_z, N_t)
        C_all_i = np.maximum(C_all_i, 0.0)
        C_outlet[i] = C_all_i[-1, :]   # last cell = outlet

        q_all_i = y_out[2 * block + i * n_z: 2 * block + (i + 1) * n_z, :]
        q_all_i = np.maximum(q_all_i, 0.0)
        q_avg[i] = np.mean(q_all_i, axis=0)

    # ── Gradient profile at output times ──
    # value_at_time returns float | ndarray; time is always ndarray here so the
    # result is ndarray, but the union confuses mypy at the GradientElutionResult
    # call site below. asarray narrows it without a runtime cost.
    gradient_profile = np.asarray(gradient.value_at_time(time))

    # ── UV signal (summed over components) ──
    if extinction_coeffs is None:
        extinction_coeffs = np.full(n_comp, 36000.0)
    extinction_coeffs = np.asarray(extinction_coeffs, dtype=float)

    uv_total = np.zeros(N_t)
    for i in range(n_comp):
        uv_total += compute_uv_signal(C_outlet[i], extinction_coeff=extinction_coeffs[i])
    uv_signal = apply_detector_broadening(uv_total, time, sigma_detector)

    # ── Peak analysis ──
    peaks: list[PeakInfo] = []
    total_areas = np.array([float(np.trapezoid(C_outlet[i], time)) for i in range(n_comp)])

    for i in range(n_comp):
        t_apex, h_apex, area_i, width_i = _find_peak(time, C_outlet[i])

        # Yield: recovered mass / injected mass
        # Injected mass approximated as C_feed * flow_rate * total_time
        mass_in = float(C_feed[i]) * flow_rate * total_time
        mass_out = area_i * flow_rate
        yield_i = min(mass_out / mass_in, 1.0) if mass_in > 0 else 0.0

        # Purity: fraction of this component's area in summed area
        sum_all_areas = total_areas.sum()
        purity_i = area_i / sum_all_areas if sum_all_areas > 0 else 0.0

        peaks.append(PeakInfo(
            component_idx=i,
            retention_time=t_apex,
            peak_height=h_apex,
            peak_area=area_i,
            peak_width_half=width_i,
            yield_fraction=yield_i,
            purity=purity_i,
        ))

    # ── Resolution between consecutive peaks ──
    resolution = np.zeros(max(n_comp - 1, 0))
    for j in range(n_comp - 1):
        pk1, pk2 = peaks[j], peaks[j + 1]
        resolution[j] = _compute_resolution(
            pk1.retention_time, pk1.peak_width_half,
            pk2.retention_time, pk2.peak_width_half,
        )

    # ── Overall yield (average across components) ──
    overall_yield = float(np.mean([p.yield_fraction for p in peaks]))

    # ── Mass balance errors ──
    V_col = column.bed_volume
    mass_balance_errors = np.zeros(n_comp)
    for i in range(n_comp):
        mass_in = float(C_feed[i]) * flow_rate * total_time
        if mass_in > 0:
            mass_out = float(np.trapezoid(C_outlet[i] * flow_rate, time))
            # Last-time-step bound mass
            q_bound_i = q_avg[i, -1] * (1.0 - eps_p) * (1.0 - eps_b) * V_col
            mass_balance_errors[i] = abs(mass_in - mass_out - q_bound_i) / mass_in

    logger.info(
        "Gradient elution complete: %d peaks, Rs=%s, overall_yield=%.1f%%",
        len(peaks),
        [f"{r:.2f}" for r in resolution],
        overall_yield * 100,
    )

    # Manifest: gradient elution has no FMC argument in this signature, so the
    # caller upstream (UI/CLI) is responsible for providing one if available.
    # We still build a manifest so the result is uniformly evidence-tagged;
    # mass-balance gate uses the worst per-component error.
    worst_mb = float(np.max(mass_balance_errors)) if len(mass_balance_errors) > 0 else 0.0
    manifest = _build_m3_chrom_manifest(
        model_basename="M3.gradient_elution.LRM",
        isotherm=isotherm,
        fmc=None,  # Not currently passed to gradient API; caller supplies if known.
        worst_mass_balance_error=worst_mb,
        diagnostics_extra={
            "n_components": int(n_comp),
            "overall_yield": float(overall_yield),
            "gradient_affects_binding": bool(gradient_affects_binding),
            "n_peaks": int(len(peaks)),
        },
    )

    return GradientElutionResult(
        time=time,
        C_outlet=C_outlet,
        gradient_profile=gradient_profile,
        uv_signal=uv_signal,
        peaks=peaks,
        resolution=resolution,
        overall_yield=overall_yield,
        pressure_drop=column.pressure_drop(flow_rate, mu),
        mass_balance_errors=mass_balance_errors,
        gradient_affects_binding=gradient_affects_binding,
        model_manifest=manifest,
    )
