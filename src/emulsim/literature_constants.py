"""Literature-reported values for simulation constants.

Each constant has a 'literature' preset with value, source, and confidence,
plus guidance on when custom calibration is recommended.

References are from peer-reviewed publications only.
"""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class LiteratureValue:
    """A constant value from published literature."""
    value: float
    unit: str
    source: str
    doi: str
    confidence: str  # "high", "medium", "low"
    notes: str = ""
    system_match: str = ""  # how well this matches our specific system


# ─── Constant 1: K_L (Langmuir adsorption constant for Span-80) ──────────

K_L = LiteratureValue(
    value=0.75,
    unit="m³/mol",
    source="Santini et al. (2007) Colloids Surf. A 298:12-21",
    doi="10.1016/j.colsurfa.2006.12.004",
    confidence="medium",
    notes=(
        "Fitted from pendant-drop tensiometry of Span-80 at paraffin oil/water "
        "interface at 25°C. Gamma_inf = 3.5e-6 mol/m². Our system operates at "
        "90°C where adsorption may differ. K_L calibrated to give σ ≈ 5 mN/m "
        "at 2% w/v Span-80, consistent with Opawale & Burgess (1998)."
    ),
    system_match="Partial — same oil/surfactant, different T and aqueous phase",
)

GAMMA_INF = LiteratureValue(
    value=3.5e-6,
    unit="mol/m²",
    source="Santini et al. (2007) Colloids Surf. A 298:12-21",
    doi="10.1016/j.colsurfa.2006.12.004",
    confidence="medium",
    notes=(
        "Maximum surface excess concentration of Span-80. Consistent across "
        "multiple studies (Peltonen & Yliruusi 2000; Aveyard et al. 2003). "
        "Range in literature: 2.5-5.0 × 10⁻⁶ mol/m²."
    ),
    system_match="Good — surfactant property, weakly system-dependent",
)

# ─── Constant 2: Chitosan intrinsic viscosity ─────────────────────────────

ETA_INTR_CHIT = LiteratureValue(
    value=800.0,
    unit="mL/g",
    source="Rinaudo et al. (1993) Int. J. Biol. Macromol. 15:281-285",
    doi="10.1016/0141-8130(93)90027-J",
    confidence="low",
    notes=(
        "Typical intrinsic viscosity for high-MW (300 kDa), high-DD (90%) "
        "chitosan in 0.2M AcOH/0.1M NaOAc at 30°C. Mark-Houwink: "
        "[η] = K·M^a with K=16.8×10⁻³ mL/g, a=0.81 for DD=91% "
        "(Roberts & Domszy 1982, DOI:10.1016/0141-8130(82)90096-7). "
        "Value is HIGHLY dependent on MW, DD, solvent, and temperature. "
        "At 90°C (emulsification), viscosity will be significantly lower."
    ),
    system_match="Poor — measured at 30°C, not 90°C; MW unknown for user's chitosan",
)

# ─── Constant 3: Breakage viscous correction C3 ──────────────────────────

BREAKAGE_C3 = LiteratureValue(
    value=0.2,
    unit="-",
    source="Alopaeus et al. (2002) Chem. Eng. Sci. 57:1815-1825",
    doi="10.1016/S0009-2509(02)00067-2",
    confidence="low",
    notes=(
        "Alopaeus viscous correction constant. Range 0.1-0.7 in the original "
        "paper, fitted to toluene/water and kerosene/water systems. No published "
        "values for paraffin oil/polymer solution systems. The value is strongly "
        "system-dependent — breakage dynamics differ for viscous biopolymer "
        "dispersed phases vs. simple organic solvents."
    ),
    system_match="Poor — different oil/dispersed phase system",
)

# ─── Constant 4: Pore structure empirical coefficients ────────────────────

PORE_A = LiteratureValue(
    value=600e-9,
    unit="m",
    source="Pernodet et al. (1997) Electrophoresis 18:55-58",
    doi="10.1002/elps.1150180111",
    confidence="medium",
    notes=(
        "Pre-factor in d_pore = A·c^α·(dT/dt)^β. Derived from AFM "
        "measurements: 2% agarose ≈ 300 nm, 4% ≈ 150 nm. The power law "
        "exponent γ = 0.5-0.75 from Pernodet (AFM) and Narayanan (2006). "
        "Chen et al. (2017) J. Sep. Sci. 40:4467 showed cooling rate effect "
        "with a two-stage relationship (transition at 6°C/min). Our coefficient "
        "of -0.2 for cooling is conservative."
    ),
    system_match="Good for agarose-only; unknown effect of chitosan addition",
)

PORE_ALPHA = LiteratureValue(
    value=-0.7,
    unit="-",
    source="Pernodet et al. (1997) Electrophoresis 18:55-58; De Gennes theory",
    doi="10.1002/elps.1150180111",
    confidence="medium",
    notes=(
        "Concentration exponent. Ogston model predicts -0.5 (random fibers), "
        "De Gennes scaling gives -0.75 (flexible chains). AFM data from "
        "Pernodet gives γ ≈ 0.59. Our value of -0.7 is within the theoretical "
        "bounds and consistent with the semi-flexible agarose helix bundles."
    ),
    system_match="Good for agarose-only gels",
)

PORE_BETA = LiteratureValue(
    value=-0.2,
    unit="-",
    source="Chen et al. (2017) J. Sep. Sci. 40:4467-4474; Aymard et al. (2001)",
    doi="10.1002/jssc.201700546",
    confidence="medium",
    notes=(
        "Cooling rate exponent. Chen et al. showed a power law relationship "
        "below 6°C/min, plateauing above. Aymard et al. (2001) demonstrated "
        "that thermal history affects agarose gel structure and mechanics. "
        "The exponent -0.2 is a conservative estimate; the actual relationship "
        "may be more complex (two-stage)."
    ),
    system_match="Good for agarose gels in the 2-10°C/min range",
)

# ─── Constant 5: IPN coupling coefficient ─────────────────────────────────

ETA_COUPLING = LiteratureValue(
    value=-0.15,
    unit="-",
    source="Estimated from Gong (2010) Soft Matter 6:2583; Nakajima et al. (2009)",
    doi="10.1039/B911622D",
    confidence="low",
    notes=(
        "No published value exists for agarose-chitosan IPN coupling. The value "
        "-0.15 is estimated from DN gel literature: true DN gels (Gong group) "
        "show synergistic effects (η > 0) due to sacrificial bonds, but our "
        "sequential IPN without sacrificial bonds should show antagonistic "
        "coupling (η < 0) from mutual swelling constraint. "
        "DeepakKumar et al. (2021, PMC7960570) showed IPN hydrogels with "
        "enhanced mechanical properties, but specific coupling coefficients "
        "are not reported in standard form."
    ),
    system_match="Poor — no published data for this specific system",
)


# ─── Summary table for all constants ──────────────────────────────────────

ALL_CONSTANTS = {
    "K_L": K_L,
    "Gamma_inf": GAMMA_INF,
    "eta_intr_chit": ETA_INTR_CHIT,
    "breakage_C3": BREAKAGE_C3,
    "pore_A": PORE_A,
    "pore_alpha": PORE_ALPHA,
    "pore_beta": PORE_BETA,
    "eta_coupling": ETA_COUPLING,
}


def summary_table() -> str:
    """Generate a markdown summary table of all literature constants."""
    lines = [
        "| Constant | Value | Unit | Confidence | System Match | Source |",
        "|----------|-------|------|------------|--------------|--------|",
    ]
    for name, lv in ALL_CONSTANTS.items():
        val_str = f"{lv.value:.2e}" if abs(lv.value) < 0.01 or abs(lv.value) > 1000 else f"{lv.value}"
        lines.append(
            f"| {name} | {val_str} | {lv.unit} | {lv.confidence} | {lv.system_match[:30]}... | {lv.source[:40]}... |"
        )
    return "\n".join(lines)
