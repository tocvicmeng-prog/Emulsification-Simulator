"""Node F1-c Phase 1 (v8.2-alpha): PLGA grade registry.

Four canonical PLGA variants used in solvent-evaporation
microsphere fabrication, parameterised from the literature anchors
in :file:`docs/f1c_plga_protocol.md` §5:

- **PLGA 50:50** — most common drug-delivery grade; fastest
  hydrolytic degradation; lower T_g.
- **PLGA 75:25** — intermediate lactide-rich formulation.
- **PLGA 85:15** — slow-release grade with higher T_g.
- **PLA** — pure poly(L-lactide); slowest degradation; highest T_g.

References
----------
Wang & Schwendeman (1999) *J. Pharm. Sci.* 88:1090 — DCM effective
    diffusivity in PLGA/DCM matrix (~10⁻¹⁰ to 10⁻⁹ m²/s).
Park et al. (1998) *Biomaterials* 19:745 — PLGA glass transition
    and modulus vs L:G ratio and M_n.
Freitas et al. (2005) *Int. J. Pharm.* 282:1 — microsphere
    fabrication review; Henry-law DCM partitioning and typical
    process phi_PLGA_0 loadings.

All four presets ship in Phase 1 (the grade-switching surface is
trivially implemented as a registry-plus-apply_preset pattern —
no need to defer three of them to a Phase 3 release).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class PLGAGradeProfile:
    """Parameter bundle for a single PLGA grade.

    Consumed by :func:`apply_preset` below to patch a
    ``MaterialProperties`` instance before the L2 solvent-evaporation
    solver runs.
    """

    name: str
    L_fraction: float   # lactide fraction (0.5 = PLGA 50:50)
    M_n: float          # [g/mol] number-average molecular weight
    T_g_C: float        # [°C] glass transition temperature
    D_DCM: float        # [m²/s] DCM effective diffusivity in gel
    phi_DCM_eq: float   # [-] DCM fraction in water sink at equilibrium
    G_glassy: float     # [Pa] glassy-state shear modulus
    n_plga: float       # [-] Gibson-Ashby exponent (porous-solid scaling)
    phi_PLGA_0_typical: float  # [-] typical initial polymer loading
    notes: str


PLGA_50_50 = PLGAGradeProfile(
    name="PLGA 50:50",
    L_fraction=0.5,
    M_n=30_000.0,
    T_g_C=45.0,
    D_DCM=1.0e-9,
    phi_DCM_eq=0.005,
    G_glassy=7.0e8,
    n_plga=2.0,
    phi_PLGA_0_typical=0.10,
    notes=(
        "Most common drug-delivery PLGA grade. Fastest hydrolytic "
        "degradation (weeks). Low T_g (~45 °C) means process must "
        "stay cool (< 30 °C bath) to avoid particle aggregation. "
        "Typical 10 wt% loading in DCM. Cost: ~USD 200-500/g (research)."
    ),
)


PLGA_75_25 = PLGAGradeProfile(
    name="PLGA 75:25",
    L_fraction=0.75,
    M_n=50_000.0,
    T_g_C=50.0,
    D_DCM=8.0e-10,
    phi_DCM_eq=0.005,
    G_glassy=9.0e8,
    n_plga=2.0,
    phi_PLGA_0_typical=0.10,
    notes=(
        "Intermediate lactide-rich formulation. Slower degradation "
        "(months). Slightly higher T_g and modulus. Common for "
        "controlled-release depot injectables with month-scale duration."
    ),
)


PLGA_85_15 = PLGAGradeProfile(
    name="PLGA 85:15",
    L_fraction=0.85,
    M_n=60_000.0,
    T_g_C=55.0,
    D_DCM=6.0e-10,
    phi_DCM_eq=0.005,
    G_glassy=1.0e9,
    n_plga=2.0,
    phi_PLGA_0_typical=0.08,
    notes=(
        "Slow-release grade. Degrades over 6-12 months. Higher crystallinity "
        "than 50:50 / 75:25; stiffer and harder. Lower DCM diffusivity "
        "means longer evaporation / extraction times."
    ),
)


PLA = PLGAGradeProfile(
    name="PLA (PLLA)",
    L_fraction=1.0,
    M_n=60_000.0,
    T_g_C=60.0,
    D_DCM=5.0e-10,
    phi_DCM_eq=0.005,
    G_glassy=1.2e9,
    n_plga=2.0,
    phi_PLGA_0_typical=0.08,
    notes=(
        "Pure poly(L-lactide). Semi-crystalline when annealed above "
        "T_g. Degradation in vivo ≥ 1 year. Used for structural "
        "applications (orthopaedic fixation, tissue scaffolds) more "
        "than drug delivery."
    ),
)


PLGA_GRADE_PRESETS: dict[str, PLGAGradeProfile] = {
    "50_50": PLGA_50_50,
    "75_25": PLGA_75_25,
    "85_15": PLGA_85_15,
    "pla": PLA,
}


def apply_preset(props, grade_name: str) -> None:
    """Patch a ``MaterialProperties`` instance with a PLGA grade preset.

    Raises
    ------
    KeyError
        If ``grade_name`` is not in :data:`PLGA_GRADE_PRESETS`.
    """
    if grade_name not in PLGA_GRADE_PRESETS:
        raise KeyError(
            f"Unknown PLGA grade {grade_name!r}. "
            f"Available: {sorted(PLGA_GRADE_PRESETS)}"
        )
    g = PLGA_GRADE_PRESETS[grade_name]
    props.D_DCM_plga = g.D_DCM
    props.phi_DCM_eq = g.phi_DCM_eq
    props.G_glassy_plga = g.G_glassy
    props.n_plga_modulus = g.n_plga
