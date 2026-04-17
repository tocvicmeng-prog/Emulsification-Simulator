"""Objective function definitions for multi-objective optimization.

Maps simulation FullResult → objective vector for BoTorch.
All objectives are MINIMIZED (lower is better).

Supports two fixed target sets (legacy / stirred-vessel) plus — as of
Node F3-a (v8.0 Phase 1) — a user-supplied ``TargetSpec`` with per-dimension
tolerances for inverse-design runs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from ..datatypes import FullResult, ModelEvidenceTier


# ─── Legacy rotor-stator targets ─────────────────────────────────────────

TARGET_D32 = 2.0e-6      # 2 µm
TARGET_PORE = 80.0e-9    # 80 nm (centre of 60-100 nm range)
TARGET_LOG_G = 4.0        # log10(G_DN) = 4 → 10 kPa

# ─── Stirred-vessel targets ──────────────────────────────────────────────
# d_mode (volume-weighted modal diameter) is used instead of d32 for
# stirred-vessel mode because: (a) downstream gelation/crosslinking acts
# on individual droplets, so the most-probable size matters more than
# the surface-area-averaged d32; (b) d32 is dominated by fines in broad
# distributions, which is misleading for the stirred-vessel regime where
# the target is 75-150 um microspheres.

TARGET_D_MODE = 100.0e-6  # 100 µm modal peak
TARGET_PORE_SV = 100.0e-9 # 100 nm (centre of 80-120 nm range)

# ─── Constraint thresholds ───────────────────────────────────────────────

MAX_SPAN = 2.0
MIN_G_DN = 1e3            # 1 kPa minimum


def compute_objectives(result: FullResult, mode: str | None = None) -> np.ndarray:
    """Compute 3 objective values from a full pipeline result.

    Legacy mode:
        f1 = |d32 - 2 µm| / 2 µm
        f2 = |pore - 80 nm| / 80 nm
        f3 = |log10(G_DN) - 4.0|

    Stirred-vessel mode:
        f1 = |d_mode - 100 µm| / 100 µm
        f2 = |pore - 100 nm| / 100 nm
        f3 = |log10(G_DN) - 4.0|

    All minimised.

    If ``mode`` is None, it is inferred from the result's parameters.
    """
    if mode is None:
        mode = getattr(result.parameters.emulsification, 'mode', 'rotor_stator_legacy')

    pore = result.gelation.pore_size_mean
    G_DN = max(result.mechanical.G_DN, 1.0)

    if mode == "stirred_vessel":
        d_mode = getattr(result.emulsification, 'd_mode', 0.0)
        f1 = abs(d_mode - TARGET_D_MODE) / TARGET_D_MODE
        f2 = abs(pore - TARGET_PORE_SV) / TARGET_PORE_SV
    else:
        d32 = result.emulsification.d32
        f1 = abs(d32 - TARGET_D32) / TARGET_D32
        f2 = abs(pore - TARGET_PORE) / TARGET_PORE

    f3 = abs(np.log10(G_DN) - TARGET_LOG_G)

    return np.array([f1, f2, f3])


# ── v6.1: Trust-aware optimization ───────────────────────────────────────

# Penalty values added to ALL objectives when evidence tier is weak.
#
# Calibration to engine REF_POINT (engine.py REF_POINT = 5.0 for all 3
# objectives). The Pareto filter at engine.py:338 excludes any candidate whose
# max objective exceeds REF_POINT, so the penalty must push QUALITATIVE_TREND
# and UNSUPPORTED candidates *above* 5.0 in at least one objective for them
# to drop off the front. Base objectives are dimensionless ratios bounded
# roughly [0, 3]; with the values below:
#
#   VALIDATED_QUANTITATIVE: max objective ≈ base − 0.05    -> kept (preferred)
#   CALIBRATED_LOCAL:       max objective ≈ base           -> kept
#   SEMI_QUANTITATIVE:      max objective ≈ base + 0.1     -> kept (mild deprioritization)
#   QUALITATIVE_TREND:      max objective ≈ base + 5.5     -> EXCLUDED (>REF_POINT)
#   UNSUPPORTED:            max objective ≈ base + 10.0    -> EXCLUDED (>>REF_POINT)
#
# Node 6 (consensus plan §5): bumped QUALITATIVE_TREND from 1.0 -> 5.5 so
# the documented "effectively excluded from Pareto" behaviour actually holds.
_TRUST_PENALTIES = {
    ModelEvidenceTier.VALIDATED_QUANTITATIVE: -0.05,
    ModelEvidenceTier.CALIBRATED_LOCAL: 0.0,
    ModelEvidenceTier.SEMI_QUANTITATIVE: 0.1,
    ModelEvidenceTier.QUALITATIVE_TREND: 5.5,
    ModelEvidenceTier.UNSUPPORTED: 10.0,
}


def trust_penalty_for_tier(tier: ModelEvidenceTier) -> float:
    """Public lookup for the per-tier optimizer penalty (Node 6).

    Returns the additive penalty applied to each objective component for
    candidates whose weakest evidence tier is `tier`. Defaults to the
    SEMI_QUANTITATIVE penalty for unknown tiers (defensive).
    """
    return _TRUST_PENALTIES.get(tier, _TRUST_PENALTIES[ModelEvidenceTier.SEMI_QUANTITATIVE])


def compute_trust_penalty(result: FullResult) -> float:
    """Return a scalar penalty based on the weakest evidence tier in the run.

    Added to each objective component, so candidates with weaker evidence
    are deprioritized in the Pareto front. A run without a RunReport is
    treated as SEMI_QUANTITATIVE (the v6.1 default for unmanaged outputs).
    """
    rr = getattr(result, "run_report", None)
    if rr is None:
        return _TRUST_PENALTIES[ModelEvidenceTier.SEMI_QUANTITATIVE]
    return trust_penalty_for_tier(rr.compute_min_tier())


def compute_objectives_trust_aware(
    result: FullResult, mode: str | None = None,
) -> np.ndarray:
    """Compute trust-penalized objectives (v6.1).

    Same as compute_objectives but adds trust penalty to each component.
    QUALITATIVE_TREND and UNSUPPORTED candidates land above the engine
    REF_POINT and are dropped from the Pareto front by the feasible_mask.
    """
    base = compute_objectives(result, mode=mode)
    penalty = compute_trust_penalty(result)
    return base + penalty


# ═════════════════════════════════════════════════════════════════════════
# Node F3-a (v8.0 Phase 1): Inverse-design target-spec objectives
# ═════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class TargetSpec:
    """User-supplied target specification for inverse-design BO.

    Each target field pairs a target value with a tolerance. Objectives
    are computed as tolerance-normalised absolute distances so every
    dimension lands on a comparable [0, ~few] scale regardless of
    physical units. Setting a target or its tolerance to ``None``
    removes that dimension from the objective vector.

    Attributes:
        d32_target: Target Sauter mean diameter [m]. Use None to
            target d_mode instead (stirred-vessel mode).
        d32_tol: Tolerance [m]; used as the scale for
            ``|d32 - target| / tol``.
        pore_target: Target pore size [m].
        pore_tol: Pore tolerance [m].
        G_DN_target: Target shear modulus [Pa]. Distance uses
            ``|log10(G_DN) - log10(target)| / log10_tol``.
        G_DN_log10_tol: log10-scale tolerance for G_DN.
        Kav_target: Target distribution coefficient [-] (optional;
            requires result.m3 block).
        Kav_tol: Tolerance on Kav.
    """
    d32_target: Optional[float] = None
    d32_tol: Optional[float] = None
    pore_target: Optional[float] = None
    pore_tol: Optional[float] = None
    G_DN_target: Optional[float] = None
    G_DN_log10_tol: Optional[float] = None
    Kav_target: Optional[float] = None
    Kav_tol: Optional[float] = None

    def active_dims(self) -> list[str]:
        """Return the list of dimensions that will contribute to the
        objective vector (i.e. both target and tol set)."""
        dims: list[str] = []
        if self.d32_target is not None and self.d32_tol is not None:
            dims.append("d32")
        if self.pore_target is not None and self.pore_tol is not None:
            dims.append("pore")
        if self.G_DN_target is not None and self.G_DN_log10_tol is not None:
            dims.append("G_DN")
        if self.Kav_target is not None and self.Kav_tol is not None:
            dims.append("Kav")
        return dims

    def validate(self) -> None:
        """Raise ValueError if the spec has no active dimension or any
        tolerance is non-positive."""
        dims = self.active_dims()
        if not dims:
            raise ValueError(
                "TargetSpec has no active dimension; at least one of "
                "(d32, pore, G_DN, Kav) must be specified with its tolerance."
            )
        for name, tol in (
            ("d32_tol", self.d32_tol),
            ("pore_tol", self.pore_tol),
            ("G_DN_log10_tol", self.G_DN_log10_tol),
            ("Kav_tol", self.Kav_tol),
        ):
            if tol is not None and tol <= 0:
                raise ValueError(f"{name} must be > 0, got {tol}")


def compute_inverse_design_objectives(
    result: FullResult,
    target: TargetSpec,
    trust_aware: bool = True,
    mode: str | None = None,
) -> np.ndarray:
    """Compute user-target-normalised objectives for inverse design.

    Each active dimension contributes ``|value - target| / tol`` (or the
    log10 variant for G_DN). Trust penalty (Node 6) is added to every
    dimension when ``trust_aware=True`` so weak-evidence candidates
    land above the engine REF_POINT and drop off the Pareto front.

    Parameters
    ----------
    result : FullResult
        Pipeline result whose emulsification/gelation/crosslinking/
        mechanical blocks provide the observables.
    target : TargetSpec
        User-supplied target specification. Must have at least one
        active dimension; ``target.validate()`` is called.
    trust_aware : bool, default True
        When True, each objective component has the Node 6 trust
        penalty added so the Pareto filter can drop weak-evidence
        candidates.
    mode : str, optional
        Pipeline mode ("rotor_stator_legacy" | "stirred_vessel"). When
        ``d32_target`` is set the d32 attribute is used; when it is
        None the function substitutes d_mode for stirred-vessel mode,
        matching the convention from the fixed-target functions above.
    """
    target.validate()
    if mode is None:
        mode = getattr(result.parameters.emulsification, "mode", "rotor_stator_legacy")

    components: list[float] = []

    if target.d32_target is not None and target.d32_tol is not None:
        if mode == "stirred_vessel":
            d_val = getattr(result.emulsification, "d_mode", None)
            if d_val is None:
                d_val = result.emulsification.d32
        else:
            d_val = result.emulsification.d32
        components.append(abs(d_val - target.d32_target) / target.d32_tol)

    if target.pore_target is not None and target.pore_tol is not None:
        pore = result.gelation.pore_size_mean
        components.append(abs(pore - target.pore_target) / target.pore_tol)

    if target.G_DN_target is not None and target.G_DN_log10_tol is not None:
        G_DN = max(result.mechanical.G_DN, 1.0)
        components.append(
            abs(np.log10(G_DN) - np.log10(max(target.G_DN_target, 1.0)))
            / target.G_DN_log10_tol
        )

    if target.Kav_target is not None and target.Kav_tol is not None:
        # Kav lives on the M3 result block if present. When absent,
        # emit +inf so the candidate is dropped by the feasibility
        # filter — a user asking for Kav-target on a run that doesn't
        # produce M3 is a configuration error.
        m3 = getattr(result, "m3", None)
        if m3 is None:
            components.append(float("inf"))
        else:
            kav = getattr(m3, "Kav", None)
            if kav is None:
                components.append(float("inf"))
            else:
                components.append(abs(kav - target.Kav_target) / target.Kav_tol)

    obj = np.array(components, dtype=float)
    if trust_aware:
        obj = obj + compute_trust_penalty(result)
    return obj


def check_constraints(result: FullResult, mode: str | None = None) -> tuple[bool, list[str]]:
    """Check optimisation constraints. Returns (feasible, violations).

    Constraints checked:
      1. Span < 2.0 (polydispersity limit)
      2. G_DN > 1 kPa (minimum mechanical integrity)
      3. d_mode in [75, 150] µm for stirred_vessel mode
      4. RPM <= 20000 (validated calibration range for L1 PBE)
      5. Crosslinker/NH2 ratio >= 0.02 (minimum stoichiometry to avoid
         severely crosslinker-limited regime)
      6. Total polymer (c_agarose + c_chitosan) <= 120 kg/m³ (pumpability
         limit, equivalent to 12% w/v)
      7. Pore size >= mesh size xi_final (inaccessible pore structure check)
    """
    if mode is None:
        mode = getattr(result.parameters.emulsification, 'mode', 'rotor_stator_legacy')

    params = result.parameters
    violations = []

    # ── Original constraints ──────────────────────────────────────────────

    if result.emulsification.span > MAX_SPAN:
        violations.append(f"span={result.emulsification.span:.2f} > {MAX_SPAN}")

    if result.mechanical.G_DN < MIN_G_DN:
        violations.append(f"G_DN={result.mechanical.G_DN:.0f} Pa < {MIN_G_DN} Pa")

    # Stirred-vessel: check d_mode is in ideal range [75, 150] µm
    if mode == "stirred_vessel":
        d_mode = getattr(result.emulsification, 'd_mode', 0.0)
        if d_mode < 75e-6 or d_mode > 150e-6:
            violations.append(
                f"d_mode={d_mode*1e6:.1f} µm outside ideal range [75, 150] µm"
            )

    # ── Additional feasibility constraints (F8b) ──────────────────────────

    # 4. RPM monotonicity guard: nonphysical regime until L1 fully recalibrated
    if params.emulsification.rpm > 20000:
        violations.append("RPM exceeds validated calibration range (max 20000)")

    # 5. Crosslinker stoichiometry check: reject if crosslinker/NH2 < 0.02
    from ..level3_crosslinking.solver import available_amine_concentration
    # Use MaterialProperties defaults for DDA and M_GlcN when not available on params
    _DDA = 0.90       # MaterialProperties default: degree of deacetylation
    _M_GlcN = 161.16  # MaterialProperties default: glucosamine molar mass [g/mol]
    NH2 = available_amine_concentration(
        params.formulation.c_chitosan, _DDA, _M_GlcN
    )
    if NH2 > 0 and params.formulation.c_genipin / NH2 < 0.02:
        violations.append("Crosslinker/NH2 ratio < 0.02 (severely limited)")

    # 6. Total polymer concentration cap: > 120 kg/m³ is extremely viscous
    total_polymer = params.formulation.c_agarose + params.formulation.c_chitosan
    if total_polymer > 120.0:
        violations.append(
            f"Total polymer {total_polymer:.0f} kg/m3 exceeds pumpability limit (120)"
        )

    # 7. Pore-mesh size consistency: pore must be accessible through crosslinked mesh
    if result.gelation.pore_size_mean < result.crosslinking.xi_final:
        violations.append("Pore size smaller than crosslinked mesh size (inaccessible)")

    return len(violations) == 0, violations


# ─── Parameter bounds (legacy — rotor-stator) ────────────────────────────

PARAM_NAMES = [
    "RPM", "c_span80", "agarose_frac", "T_oil",
    "cooling_rate", "c_genipin", "t_crosslink",
]

PARAM_BOUNDS = np.array([
    [3000.0,  25000.0],    # RPM
    [5.0,     50.0],       # c_span80 [kg/m³]
    [0.6,     0.9],        # agarose fraction
    [333.15,  368.15],     # T_oil [K] (60-95°C)
    [0.033,   0.333],      # cooling_rate [K/s] (2-20°C/min)
    [0.5,     10.0],       # c_genipin [mol/m³]
    [3600.0,  172800.0],   # t_crosslink [s] (1-48 h)
])

# Indices for log-scale parameters
LOG_SCALE_INDICES = [0, 1, 4, 5, 6]  # RPM, span80, cooling, genipin, t_xlink

# ─── Parameter bounds (stirred-vessel) ───────────────────────────────────

PARAM_NAMES_SV = [
    "RPM", "c_span80", "agarose_frac", "T_oil",
    "cooling_rate", "c_genipin", "t_crosslink",
]

PARAM_BOUNDS_SV = np.array([
    [800.0,   2000.0],     # RPM (Stirrer A: pitched blade)
    [5.0,     30.0],       # c_span80 [kg/m³] (1.5% v/v default)
    [0.6,     0.9],        # agarose fraction
    [333.15,  358.15],     # T_oil [K] (60-85°C, limited by heating)
    [0.005,   0.05],       # cooling_rate [K/s] (0.3-3°C/min, natural cooling)
    [0.5,     10.0],       # c_genipin [mol/m³]
    [3600.0,  172800.0],   # t_crosslink [s] (1-48 h)
])

PARAM_BOUNDS_SV_RS = np.array([
    [800.0,   9000.0],     # RPM (Stirrer B: small rotor-stator)
    [5.0,     30.0],       # c_span80 [kg/m³]
    [0.6,     0.9],        # agarose fraction
    [333.15,  358.15],     # T_oil [K]
    [0.005,   0.05],       # cooling_rate [K/s]
    [0.5,     10.0],       # c_genipin [mol/m³]
    [3600.0,  172800.0],   # t_crosslink [s]
])

LOG_SCALE_INDICES_SV = [0, 1, 4, 5, 6]


def get_param_bounds(mode: str = "rotor_stator_legacy",
                     stirrer_type: str = "pitched_blade") -> np.ndarray:
    """Return parameter bounds appropriate for the given mode and stirrer."""
    if mode == "stirred_vessel":
        if stirrer_type == "rotor_stator_small":
            return PARAM_BOUNDS_SV_RS.copy()
        return PARAM_BOUNDS_SV.copy()
    return PARAM_BOUNDS.copy()
