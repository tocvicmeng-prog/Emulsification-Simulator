"""Node F1-a Phase 2c (v8.0-beta): Alginate gelant / Ca²⁺-source library.

Alginate ionic gelation is driven by Ca²⁺ availability, not by a
covalent crosslinker. Two practical source modes are supported here:

  1. **External bath** — droplets are extruded or pre-formed and then
     dropped into a CaCl₂ bath. Ca²⁺ diffuses in from the boundary,
     egg-box junctions propagate inward as a shrinking front. This is
     the default mode solved by `level2_gelation.ionic_ca`.

  2. **Internal release** — a sparingly-soluble Ca source (e.g.,
     CaCO₃) is dispersed into the alginate phase before emulsification.
     A slow-dissolving acidifier (glucono-δ-lactone, GDL) then
     hydrolyses to gluconic acid, dropping the pH and releasing Ca²⁺
     in situ. This produces a far more homogeneous gel than diffusion
     from the outside and is the standard route for "structured"
     alginate foods/biomaterials.

References
----------
External CaCl₂ bath:
    Kuo & Ma (2001) *Biomaterials* 22:511-521 — Ca²⁺ diffusion
        coefficients and bath-concentration scaling in 1-3 mm beads.
    Draget et al. (1997) *Carbohydr. Polym.* 32:79 — 1-500 mM CaCl₂
        screen, gel homogeneity vs bath strength.

Internal GDL/CaCO₃:
    Draget, Skjåk-Bræk, Smidsrød (1997) *Int. J. Biol. Macromol.*
        21:47 — canonical paper on GDL + CaCO₃ for homogeneous
        alginate gels; GDL hydrolysis rate at 25 °C.
    Liu et al. (2003) *Food Hydrocoll.* 17:661 — CaCO₃ particle size
        vs gel strength; stoichiometry 2 COO⁻ : 1 Ca²⁺ as usual.

Coupling to the simulator
-------------------------
The L2 ionic-Ca solver accepts a single scalar `C_Ca_bath` [mol/m³].
For external mode this is literally the bath concentration. For
internal mode the simulator uses an *effective* time-averaged Ca²⁺
concentration ≈ (C_Ca_source · (1 − exp(−k_release · t_end))) — a
coarse lumped-parameter approximation that holds when GDL hydrolysis
is slow compared to Ca²⁺ binding. A future Phase 3 solver variant
could solve coupled GDL + CaCO₃ + alginate, but Phase 2c keeps the
interface minimal and exposes the two profiles as presets only.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class AlginateGelantProfile:
    """Profile for a Ca²⁺ source for alginate ionic gelation.

    Phase 2c stores only the parameters needed to build a
    ``SimulationParameters.formulation.c_Ca_bath`` value plus a
    literature-backed default time scale. Downstream solvers
    consume these via the existing ``solve_ionic_ca_gelation`` API.
    """

    name: str
    cas: str
    mode: str  # "external_bath" | "internal_release"
    # External-mode parameter
    C_Ca_bath: float  # [mol/m³] bath Ca²⁺ concentration (external mode)
    # Internal-mode parameters (ignored if mode="external_bath")
    C_Ca_source: float  # [mol/m³] total Ca²⁺ available from dispersed source
    k_release: float  # [1/s] effective first-order Ca²⁺ release rate
    # Process defaults
    T_default: float  # [K] recommended bath / process temperature
    t_default: float  # [s] recommended gelation time
    suitability: int  # 1-10 score for microsphere production
    notes: str


# ═══════════════════════════════════════════════════════════════════════
#  ALGINATE GELANT LIBRARY
# ═══════════════════════════════════════════════════════════════════════

GELANTS_ALGINATE: dict[str, AlginateGelantProfile] = {
    # ── 1. CaCl₂ external bath (baseline, shrinking-core) ──────────────
    "cacl2_external": AlginateGelantProfile(
        name="CaCl₂ external bath",
        cas="10043-52-4",
        mode="external_bath",
        # Typical industrial range 50-200 mM; 100 mM is the standard
        # default in most alginate microsphere recipes (Draget 1997).
        C_Ca_bath=100.0,
        C_Ca_source=0.0,
        k_release=0.0,
        T_default=298.15,  # 25 °C — no thermal activation needed
        t_default=1800.0,  # 30 min typical for ≤ 500 µm beads
        suitability=8,
        notes=(
            "Classical route: pre-formed alginate droplets (via extrusion, "
            "emulsification, or electrospray) fall into a stirred CaCl₂ "
            "bath and gel from the outside in. Shrinking-core kinetics "
            "— gel front depth ∝ √(D·t). Fast and robust, but produces "
            "an inhomogeneous gel (denser outer shell, softer core) for "
            "beads > 500 µm. Cost: ~USD 0.01/g CaCl₂·2H₂O."
        ),
    ),
    # ── 2. GDL + CaCO₃ internal release (homogeneous) ──────────────────
    "gdl_caco3_internal": AlginateGelantProfile(
        name="GDL + CaCO₃ internal release",
        cas="90-80-2 (GDL) / 471-34-1 (CaCO₃)",
        mode="internal_release",
        C_Ca_bath=0.0,  # no external bath — Ca released from within
        # Typical loading 10-40 mM CaCO₃ equivalent in the alginate phase.
        # 20 mM is a reasonable default at ~1 % alginate.
        C_Ca_source=20.0,
        # GDL hydrolysis k ≈ 1.5e-4 s⁻¹ at 25 °C (Draget 1997);
        # CaCO₃ dissolution tracks pH drop roughly linearly, so the
        # effective Ca²⁺ release rate is GDL-limited at ~1-2e-4 s⁻¹.
        k_release=1.5e-4,
        T_default=298.15,  # 25 °C
        # GDL hydrolysis is slow; gelation plateau reached at ~4-6 h.
        t_default=14400.0,  # 4 h
        suitability=7,
        notes=(
            "In-situ gelation: CaCO₃ powder dispersed in the alginate "
            "phase before emulsification, GDL added just before droplet "
            "formation. GDL hydrolyses to gluconic acid (t½ ≈ 80 min at "
            "25 °C), pH drops from ~8 to ~4, CaCO₃ solubilises, Ca²⁺ is "
            "released homogeneously throughout the droplet. Result: far "
            "more uniform crosslink density than CaCl₂ bath, at the cost "
            "of much longer gelation times (hours) and strict stoichiometric "
            "balance (GDL/CaCO₃ molar ratio ~2:1). Current simulator "
            "handles internal release via the effective-bath lumped "
            "approximation C_eff = C_Ca_source·(1 − exp(−k_release·t)); "
            "full coupled GDL/CaCO₃/alginate solver = Phase 3 follow-up. "
            "Cost: ~USD 0.02/g GDL, ~USD 0.01/g CaCO₃."
        ),
    ),
}


def effective_bath_concentration(
    profile: AlginateGelantProfile, t_end: float
) -> float:
    """Lumped-parameter effective bath concentration for internal-release mode.

    For ``mode="external_bath"`` returns ``profile.C_Ca_bath`` unchanged.
    For ``mode="internal_release"`` returns
    ``C_Ca_source · (1 − exp(−k_release · t_end))`` — a monotone,
    saturating approximation to the time-averaged Ca²⁺ availability.

    Parameters
    ----------
    profile : AlginateGelantProfile
    t_end : float
        Simulation end time [s]. Used only for internal mode.

    Returns
    -------
    float
        Effective Ca²⁺ bath concentration [mol/m³] to feed into the
        ionic-Ca L2 solver.
    """
    if profile.mode == "external_bath":
        return float(profile.C_Ca_bath)
    if profile.mode == "internal_release":
        import math
        if t_end <= 0.0 or profile.k_release <= 0.0:
            return 0.0
        return float(
            profile.C_Ca_source * (1.0 - math.exp(-profile.k_release * t_end))
        )
    raise ValueError(
        f"Unknown alginate gelant mode: {profile.mode!r}. "
        f"Expected 'external_bath' or 'internal_release'."
    )
