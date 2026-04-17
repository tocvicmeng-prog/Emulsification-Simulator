"""Node F1-b Phase 1 (v8.1-alpha): cellulose NIPS solvent-system defaults.

Phase 1 ships NaOH/urea (Zhang Lab Wuhan) as the default system — widest
literature base, room-temperature processing, patent-clear. NMMO,
EMIM-Ac, and DMAc/LiCl presets are stubbed for F1-b Phase 3.

References
----------
Cuissinat & Navard (2006) *Macromol. Symp.* 244:1 — dissolution in
    NaOH/urea, NMMO, IL.
Xu et al. (2010) *Biomacromolecules* 11:1724 — cellulose/NaOH-urea
    transport constants.
Lindman et al. (2010) *Phys. Chem. Chem. Phys.* 12:4369 — Flory-Huggins
    χ parameters for cellulose solvents.
Zhang et al. (2020) *Cellulose* 27:1071 — modulus vs φ_cell scaling
    for regenerated cellulose hydrogels.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CelluloseSolventPreset:
    """Parameter bundle for a single cellulose solvent / non-solvent pair.

    Consumed by ``apply_preset(mat_props, preset_name)`` below to patch
    a ``MaterialProperties`` instance before the L2 NIPS solver runs.
    """

    name: str
    solvent: str
    nonsolvent: str
    N_p: float             # [-] degree of polymerisation of cellulose
    chi_PS: float          # [-] χ polymer-solvent
    chi_PN: float          # [-] χ polymer-nonsolvent
    chi_SN: float          # [-] χ solvent-nonsolvent
    D_solvent: float       # [m²/s] solvent self-diffusion in gel
    D_nonsolvent: float    # [m²/s] non-solvent self-diffusion in gel
    kappa_CH: float        # [J·m⁻¹] Cahn-Hilliard gradient energy coef
    K_cell: float          # [Pa] modulus prefactor
    alpha_cell: float      # [-] modulus exponent
    phi_cellulose_0_typical: float  # [-] typical initial cellulose vol. frac.
    T_process: float       # [K] recommended processing temperature
    notes: str


# Canonical NaOH/urea preset (F1-b Phase 1 default). Values from §5 of
# docs/f1b_cellulose_nips_protocol.md.
NAOH_UREA = CelluloseSolventPreset(
    name="NaOH/urea",
    solvent="7 wt% NaOH + 12 wt% urea (aq.)",
    nonsolvent="water",
    N_p=370.0,
    chi_PS=0.45,
    chi_PN=0.85,
    chi_SN=0.30,
    D_solvent=5.0e-11,
    D_nonsolvent=1.0e-10,
    kappa_CH=1.0e-17,
    K_cell=5.0e5,
    alpha_cell=2.25,
    phi_cellulose_0_typical=0.05,  # ~5 wt% cellulose in solution
    T_process=298.15,
    notes=(
        "Zhang Lab Wuhan system — cellulose dissolves in pre-cooled 7 wt% "
        "NaOH + 12 wt% urea after a brief freeze-thaw. Room-temperature "
        "processing, aqueous non-solvent (water). Patent-clear for most "
        "jurisdictions. Typical M_n ≈ 60 kDa gives N_p ≈ 370."
    ),
)


# NMMO (Lyocell) — industrial Lenzing process. Higher M_n typical for
# pulp feedstock, slower solvent diffusivity at 90 °C, slightly lower
# χ_PS (better cellulose solubility) — all per Cuissinat & Navard 2006
# and Xu et al. 2010.
NMMO = CelluloseSolventPreset(
    name="NMMO",
    solvent="N-methylmorpholine N-oxide (~80 wt% aq.)",
    nonsolvent="water",
    N_p=500.0,                   # ~81 kDa M_n typical for Lyocell pulp
    chi_PS=0.40,
    chi_PN=0.80,
    chi_SN=0.20,
    D_solvent=2.0e-11,           # viscous at 90 °C
    D_nonsolvent=8.0e-11,
    kappa_CH=1.0e-17,
    K_cell=8.0e5,                # higher modulus from higher M_n
    alpha_cell=2.25,
    phi_cellulose_0_typical=0.08,  # 5-15 wt% industrial range
    T_process=363.15,             # 90 °C
    notes=(
        "Lenzing Lyocell system. Dissolves high-M_n wood-pulp cellulose "
        "in ~80 wt% NMMO/water at 90-120 °C. Non-solvent bath is water. "
        "Significant Lenzing patent overhang; prefer NaOH/urea for "
        "patent-clear R&D. Higher polymer fraction and higher modulus "
        "possible, at the cost of elevated T and moisture-control."
    ),
)


# EMIM-Ac (1-ethyl-3-methylimidazolium acetate) — research-scale ionic
# liquid solvent. Swatloski et al. 2002 / Zhu et al. 2006 — excellent
# cellulose dissolution with the lowest χ_PS of the common systems.
EMIM_AC = CelluloseSolventPreset(
    name="EMIM-Ac",
    solvent="1-ethyl-3-methylimidazolium acetate",
    nonsolvent="water",
    N_p=370.0,
    chi_PS=0.38,
    chi_PN=0.80,
    chi_SN=0.25,
    D_solvent=1.0e-11,            # viscous IL
    D_nonsolvent=5.0e-11,
    kappa_CH=1.0e-17,
    K_cell=4.0e5,
    alpha_cell=2.25,
    phi_cellulose_0_typical=0.06,
    T_process=353.15,             # 80 °C
    notes=(
        "Imidazolium-acetate ionic liquid (Swatloski 2002). Dissolves "
        "cellulose at 80-100 °C without derivatisation. Expensive "
        "(~USD 100/g), moisture-sensitive, but highly solvating — "
        "lowest χ_PS among the systems here. Water is a strong "
        "non-solvent: recovery of EMIM-Ac from the coagulation bath is "
        "the major process bottleneck and cost driver."
    ),
)


# DMAc/LiCl (N,N-dimethylacetamide + lithium chloride) — classic
# cellulose solvent for analytical chromatography (GPC, etc.). Less
# common for microsphere fabrication because of the ternary complexity
# and the water non-solvent quench step. Included here for completeness.
DMAC_LICL = CelluloseSolventPreset(
    name="DMAc/LiCl",
    solvent="N,N-dimethylacetamide + 8 wt% LiCl",
    nonsolvent="water",
    N_p=450.0,
    chi_PS=0.42,
    chi_PN=0.82,
    chi_SN=0.28,
    D_solvent=3.0e-11,
    D_nonsolvent=7.0e-11,
    kappa_CH=1.0e-17,
    K_cell=6.0e5,
    alpha_cell=2.25,
    phi_cellulose_0_typical=0.05,
    T_process=333.15,             # 60 °C activation
    notes=(
        "McCormick system (DMAc/LiCl, pre-activated with hot DMAc or "
        "polar swelling). Widely used for cellulose GPC/MALLS analytics "
        "but rare for microsphere fabrication — LiCl residuals are hard "
        "to wash out and may perturb downstream chromatography. "
        "Quaternary dissolution mechanism: Li⁺ coordinates to DMAc, "
        "Cl⁻ H-bonds to cellulose -OH; water quench precipitates "
        "cellulose with moderate χ_PN."
    ),
)


# Phase 3 registry — all four systems wired.
CELLULOSE_SOLVENT_PRESETS: dict[str, CelluloseSolventPreset] = {
    "naoh_urea": NAOH_UREA,
    "nmmo": NMMO,
    "emim_ac": EMIM_AC,
    "dmac_licl": DMAC_LICL,
}


def apply_preset(props, preset_name: str) -> None:
    """Patch a ``MaterialProperties`` instance with a solvent preset.

    Mutates ``props`` in place — the usual pattern is to build default
    ``MaterialProperties`` then call ``apply_preset(props, "naoh_urea")``
    before handing the props to the L2 NIPS solver.

    Raises
    ------
    KeyError
        If ``preset_name`` is not in :data:`CELLULOSE_SOLVENT_PRESETS`.
    """
    if preset_name not in CELLULOSE_SOLVENT_PRESETS:
        raise KeyError(
            f"Unknown cellulose solvent preset {preset_name!r}. "
            f"Available: {sorted(CELLULOSE_SOLVENT_PRESETS)}"
        )
    p = CELLULOSE_SOLVENT_PRESETS[preset_name]
    props.N_p_cellulose = p.N_p
    props.chi_PS_cellulose = p.chi_PS
    props.chi_PN_cellulose = p.chi_PN
    props.chi_SN_cellulose = p.chi_SN
    props.D_solvent_cellulose = p.D_solvent
    props.D_nonsolvent_cellulose = p.D_nonsolvent
    props.kappa_CH_cellulose = p.kappa_CH
    props.K_cell_modulus = p.K_cell
    props.alpha_cell_modulus = p.alpha_cell
