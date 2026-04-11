"""Reagent profiles for Module 2 functionalization chemistry.

Phase B: Minimal Module 2 — 2 Workflows.
Architecture: module2_module3_final_implementation_plan.md, Phase B.

Only 4 reagent profiles needed for Phase B:
  - genipin_secondary: Secondary amine crosslinking (genipin)
  - glutaraldehyde_secondary: Secondary amine crosslinking (glutaraldehyde)
  - ech_activation: Hydroxyl activation with epichlorohydrin
  - dvs_activation: Hydroxyl activation with divinyl sulfone

Phase 2 expansion (Module 2 extension) adds 10 further profiles:
  Ligand coupling (4): deae_coupling, ida_coupling, phenyl_coupling, sp_coupling
  Protein coupling (2): protein_a_coupling, protein_g_coupling
  Quenching (4): ethanolamine_quench, mercaptoethanol_quench, nabh4_quench,
                 acetic_anhydride_quench

Follows the CrosslinkerProfile pattern from reagent_library.py.
Literature values sourced from the Scientific Advisor report.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .acs import ACSSiteType


@dataclass
class ReagentProfile:
    """Profile for a Module 2 functionalization reagent.

    Maps directly to the reaction engine inputs in reactions.py.
    Follows the CrosslinkerProfile pattern from reagent_library.py.

    Attributes:
        name: Human-readable reagent name.
        cas: CAS registry number.
        reaction_type: Chemistry category ("crosslinking", "activation", "blocking").
        target_acs: ACS site type consumed by this reagent.
        product_acs: ACS site type produced (None for crosslinking).
        k_forward: Forward rate constant at reference T [m^3/(mol*s)].
        E_a: Activation energy [J/mol].
        stoichiometry: Moles reagent consumed per mole ACS consumed [-].
        hydrolysis_rate: First-order hydrolysis rate constant [1/s].
        ph_optimum: Optimal pH for reaction [-].
        temperature_default: Default reaction temperature [K].
        time_default: Default reaction time [s].
        notes: Literature references and notes.

        # ── Extended identity (Phase 2) ──
        reagent_identity: Chemical name of actual reactive form.
        installed_ligand: Functional group after coupling.
        functional_mode: One of "crosslinker", "activator", "iex_ligand",
            "affinity_ligand", "hic_ligand", "imac_chelator", "quencher".
        chemistry_class: One of "epoxide_amine", "epoxide_thiol", "vs_amine",
            "vs_thiol", "aldehyde_amine", "reduction", "acetylation",
            "amine_covalent".

        # ── Validity windows ──
        ph_min: Minimum valid pH [-].
        ph_max: Maximum valid pH [-].
        temperature_min: Minimum valid temperature [K].
        temperature_max: Maximum valid temperature [K].

        # ── Macromolecule fields ──
        ligand_mw: Ligand molecular weight [Da].
        ligand_r_h: Hydrodynamic radius [m].
        is_macromolecule: True if ligand is a macromolecule (protein, etc.).
        activity_retention: Fraction of activity retained after coupling [0, 1].
        activity_retention_uncertainty: Uncertainty on activity_retention [0, 1].
        max_surface_density: Steric jamming limit [mol/m^2].

        # ── Metadata ──
        confidence_tier: "semi_quantitative" or "ranking_only".
        calibration_source: Source or rationale for kinetic parameters.
        hazard_class: GHS / lab hazard descriptor string.
    """
    name: str
    cas: str
    reaction_type: str
    target_acs: ACSSiteType
    product_acs: Optional[ACSSiteType]
    k_forward: float         # [m^3/(mol*s)] at reference T
    E_a: float               # [J/mol]
    stoichiometry: float     # [-]
    hydrolysis_rate: float = 0.0   # [1/s]
    ph_optimum: float = 7.0
    temperature_default: float = 298.15  # [K]
    time_default: float = 3600.0         # [s]
    notes: str = ""

    # ── Extended identity (Phase 2) ──
    reagent_identity: str = ""         # Chemical name of actual reactive form
    installed_ligand: str = ""         # Functional group after coupling
    functional_mode: str = ""          # "crosslinker", "activator", "iex_ligand", "affinity_ligand", "hic_ligand", "imac_chelator", "quencher"
    chemistry_class: str = ""          # "epoxide_amine", "epoxide_thiol", "vs_amine", "vs_thiol", "aldehyde_amine", "reduction", "acetylation"

    # ── Validity windows ──
    ph_min: float = 0.0
    ph_max: float = 14.0
    temperature_min: float = 273.15    # [K]
    temperature_max: float = 373.15    # [K]

    # ── Macromolecule fields ──
    ligand_mw: float = 0.0            # [Da]
    ligand_r_h: float = 0.5e-9        # [m] hydrodynamic radius
    is_macromolecule: bool = False
    activity_retention: float = 1.0    # [0,1]
    activity_retention_uncertainty: float = 0.0
    max_surface_density: float = 0.0   # [mol/m^2] steric jamming limit

    # ── Metadata ──
    confidence_tier: str = "semi_quantitative"
    calibration_source: str = ""
    hazard_class: str = ""


# ─── Phase B Reagent Library (4 profiles) ─────────────────────────────

REAGENT_PROFILES: dict[str, ReagentProfile] = {

    # ── 1. Genipin secondary crosslinking ─────────────────────────────
    # Butler et al. (2003): k ~ 0.002 L/(mol*s) = 2e-6 m^3/(mol*s) at 37 degC
    # E_a ~ 45-55 kJ/mol.  1 genipin bridges 2 NH2 groups (stoich = 0.5).
    "genipin_secondary": ReagentProfile(
        name="Genipin (secondary crosslinking)",
        cas="6902-77-8",
        reaction_type="crosslinking",
        target_acs=ACSSiteType.AMINE_PRIMARY,
        product_acs=None,
        k_forward=2e-6,          # [m^3/(mol*s)] at 37 degC
        E_a=45000.0,             # [J/mol]
        stoichiometry=0.5,       # 1 genipin per 2 NH2
        hydrolysis_rate=0.0,
        ph_optimum=7.4,
        temperature_default=310.15,  # 37 degC
        time_default=14400.0,        # 4 h
        functional_mode="crosslinker",
        chemistry_class="amine_covalent",
        ph_min=6.0,
        ph_max=9.0,
        notes=(
            "Secondary crosslinking after primary L3 genipin cure. "
            "Butler et al. (2003) J. Polym. Sci. A: Polym. Chem. 41:3941. "
            "Low cytotoxicity, FDA-approved colorant."
        ),
    ),

    # ── 2. Glutaraldehyde secondary crosslinking ──────────────────────
    # Migneault et al. (2004): k ~ 0.01 L/(mol*s) = 1e-5 m^3/(mol*s)
    # E_a ~ 35-40 kJ/mol.  1 glutaraldehyde bridges 2 NH2 (Schiff base).
    "glutaraldehyde_secondary": ReagentProfile(
        name="Glutaraldehyde (secondary crosslinking)",
        cas="111-30-8",
        reaction_type="crosslinking",
        target_acs=ACSSiteType.AMINE_PRIMARY,
        product_acs=None,
        k_forward=1e-5,          # [m^3/(mol*s)]
        E_a=40000.0,             # [J/mol]
        stoichiometry=0.5,       # 1 GA per 2 NH2 (Schiff base bridge)
        hydrolysis_rate=0.0,
        ph_optimum=7.0,
        temperature_default=298.15,  # 25 degC
        time_default=3600.0,         # 1 h
        functional_mode="crosslinker",
        chemistry_class="amine_covalent",
        ph_min=6.0,
        ph_max=8.0,
        hazard_class="toxic",
        notes=(
            "Schiff base crosslinking of primary amines. Fast kinetics. "
            "Migneault et al. (2004) BioTechniques 37:790."
        ),
    ),

    # ── 3. Epichlorohydrin (ECH) activation ───────────────────────────
    # Sundberg & Porath (1974): ECH reacts with OH under alkaline conditions
    # to introduce epoxide groups.  Significant hydrolysis at pH > 11.
    # k ~ 1.5e-5 m^3/(mol*s), E_a ~ 60 kJ/mol.
    # Hydrolysis rate ~ 1e-4 /s at pH 12, 25 degC.
    "ech_activation": ReagentProfile(
        name="Epichlorohydrin (OH activation)",
        cas="106-89-8",
        reaction_type="activation",
        target_acs=ACSSiteType.HYDROXYL,
        product_acs=ACSSiteType.EPOXIDE,
        k_forward=1.5e-5,        # [m^3/(mol*s)]
        E_a=60000.0,             # [J/mol]
        stoichiometry=1.0,       # 1 ECH per OH consumed
        hydrolysis_rate=1e-4,    # [1/s] at pH 12
        ph_optimum=12.0,
        temperature_default=298.15,  # 25 degC
        time_default=7200.0,         # 2 h
        functional_mode="activator",
        chemistry_class="epoxide_amine",
        ph_min=10.0,
        ph_max=13.0,
        hazard_class="toxic_carcinogen",
        notes=(
            "Alkaline activation of agarose hydroxyl groups. "
            "Introduces epoxide for subsequent ligand coupling. "
            "Sundberg & Porath (1974) J. Chromatogr. 90:87. "
            "Competing hydrolysis significant at pH > 11."
        ),
    ),

    # ── 4. Divinyl sulfone (DVS) activation ───────────────────────────
    # Porath & Fornstedt (1970): DVS reacts with OH to introduce vinyl
    # sulfone groups.  More stable than ECH epoxides (no hydrolysis).
    # k ~ 5e-6 m^3/(mol*s), E_a ~ 55 kJ/mol.
    "dvs_activation": ReagentProfile(
        name="Divinyl sulfone (OH activation)",
        cas="77-77-0",
        reaction_type="activation",
        target_acs=ACSSiteType.HYDROXYL,
        product_acs=ACSSiteType.VINYL_SULFONE,
        k_forward=5e-6,          # [m^3/(mol*s)]
        E_a=55000.0,             # [J/mol]
        stoichiometry=1.0,       # 1 DVS per OH consumed
        hydrolysis_rate=0.0,     # Negligible hydrolysis
        ph_optimum=11.0,
        temperature_default=298.15,  # 25 degC
        time_default=3600.0,         # 1 h
        functional_mode="activator",
        chemistry_class="vs_amine",
        ph_min=10.0,
        ph_max=12.0,
        hazard_class="toxic",
        notes=(
            "Alkaline activation of agarose hydroxyl groups. "
            "Introduces vinyl sulfone for nucleophilic coupling. "
            "Porath & Fornstedt (1970) J. Chromatogr. 51:479. "
            "More hydrolytically stable than ECH epoxides."
        ),
    ),

    # ─── Phase 2 Expansion — Ligand Coupling (4) ──────────────────────

    # ── 5. DEAE coupling ──────────────────────────────────────────────
    "deae_coupling": ReagentProfile(
        name="DEAE (weak anion exchange)",
        cas="100-36-7",
        reagent_identity="2-(Diethylamino)ethylamine",
        installed_ligand="DEAE",
        functional_mode="iex_ligand",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=5e-5,
        E_a=50000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-5,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        temperature_min=288.15, temperature_max=313.15,
        time_default=14400.0,
        ligand_mw=116.0,
        ligand_r_h=0.4e-9,
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from general epoxide-amine kinetics",
        notes="Weak anion exchanger; pKa ~11.5; fully charged below pH 9",
    ),

    # ── 6. IDA coupling ───────────────────────────────────────────────
    "ida_coupling": ReagentProfile(
        name="IDA (IMAC chelator)",
        cas="142-73-4",
        reagent_identity="Iminodiacetic acid",
        installed_ligand="IDA",
        functional_mode="imac_chelator",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=2e-5,
        E_a=45000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-5,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        temperature_min=288.15, temperature_max=313.15,
        time_default=21600.0,
        ligand_mw=133.0,
        ligand_r_h=0.4e-9,
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from general epoxide-amine kinetics",
        notes="IMAC chelator; requires metal charging (Ni/Co/Cu/Zn) for function — not modeled here",
    ),

    # ── 7. Phenyl coupling ────────────────────────────────────────────
    "phenyl_coupling": ReagentProfile(
        name="Phenyl (HIC ligand)",
        cas="62-53-3",
        reagent_identity="Phenylamine (aniline)",
        installed_ligand="Phenyl",
        functional_mode="hic_ligand",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=3e-5,
        E_a=48000.0,
        stoichiometry=1.0,
        hydrolysis_rate=5e-6,
        ph_optimum=9.5,
        ph_min=8.0, ph_max=11.0,
        temperature_default=298.15,
        temperature_min=288.15, temperature_max=313.15,
        time_default=14400.0,
        ligand_mw=93.0,
        ligand_r_h=0.3e-9,
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from aniline-epoxide kinetics",
        hazard_class="toxic",
        notes="HIC ligand; weak nucleophile; hydrophobic interaction depends on salt concentration",
    ),

    # ── 8. Sulfopropyl coupling ───────────────────────────────────────
    "sp_coupling": ReagentProfile(
        name="Sulfopropyl (strong cation exchange)",
        cas="1120-71-4",
        reagent_identity="1,3-Propane sultone",
        installed_ligand="Sulfopropyl",
        functional_mode="iex_ligand",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=4e-5,
        E_a=50000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-5,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        temperature_min=288.15, temperature_max=313.15,
        time_default=14400.0,
        ligand_mw=122.0,
        ligand_r_h=0.4e-9,
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from sultone ring-opening kinetics",
        hazard_class="toxic",
        notes="Strong cation exchanger; permanently charged sulfonate at all pH",
    ),

    # ─── Phase 2 Expansion — Protein Coupling (2) ─────────────────────

    # ── 9. Protein A coupling ─────────────────────────────────────────
    "protein_a_coupling": ReagentProfile(
        name="Protein A (IgG affinity)",
        cas="91932-65-9",
        reagent_identity="Recombinant Protein A (rSPA)",
        installed_ligand="Protein A",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=5e-7,
        E_a=25000.0,
        stoichiometry=1.0,
        hydrolysis_rate=0.0,
        ph_optimum=9.0,
        ph_min=7.5, ph_max=10.0,
        temperature_default=277.15,
        temperature_min=273.15, temperature_max=283.15,
        time_default=57600.0,
        ligand_mw=42000.0,
        ligand_r_h=2.5e-9,
        is_macromolecule=True,
        activity_retention=0.60,
        activity_retention_uncertainty=0.15,
        max_surface_density=2e-8,
        confidence_tier="ranking_only",
        calibration_source="Estimated; activity retention from Cytiva Protein A Sepharose literature",
        notes="Couple at 4C to preserve folding; 1 Protein A binds 2 IgG Fc; ranking_only unless calibrated",
    ),

    # ── 10. Protein G coupling ────────────────────────────────────────
    "protein_g_coupling": ReagentProfile(
        name="Protein G (IgG affinity, broad subclass)",
        cas="122441-07-8",
        reagent_identity="Recombinant Protein G (rSPG)",
        installed_ligand="Protein G",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=5e-7,
        E_a=25000.0,
        stoichiometry=1.0,
        hydrolysis_rate=0.0,
        ph_optimum=9.0,
        ph_min=7.5, ph_max=10.0,
        temperature_default=277.15,
        temperature_min=273.15, temperature_max=283.15,
        time_default=57600.0,
        ligand_mw=22000.0,
        ligand_r_h=2.0e-9,
        is_macromolecule=True,
        activity_retention=0.65,
        activity_retention_uncertainty=0.15,
        max_surface_density=3e-8,
        confidence_tier="ranking_only",
        calibration_source="Estimated; broader IgG subclass coverage than Protein A",
        notes="Couple at 4C; broader subclass binding; ranking_only unless calibrated",
    ),

    # ─── Phase 2 Expansion — Quenching (4) ────────────────────────────

    # ── 11. Ethanolamine quench ───────────────────────────────────────
    "ethanolamine_quench": ReagentProfile(
        name="Ethanolamine (epoxide quench)",
        cas="141-43-5",
        reagent_identity="Ethanolamine",
        installed_ligand="beta-hydroxyethylamine (blocked)",
        functional_mode="quencher",
        reaction_type="blocking",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=1e-3,
        E_a=30000.0,
        stoichiometry=1.0,
        hydrolysis_rate=0.0,
        ph_optimum=8.5,
        ph_min=7.0, ph_max=10.0,
        temperature_default=298.15,
        time_default=7200.0,
        confidence_tier="semi_quantitative",
        calibration_source="Standard Sepharose quenching protocol",
        notes="Standard quench for epoxide-activated media; 1M concentration typical",
    ),

    # ── 12. 2-Mercaptoethanol quench ──────────────────────────────────
    "mercaptoethanol_quench": ReagentProfile(
        name="2-Mercaptoethanol (VS quench)",
        cas="60-24-2",
        reagent_identity="2-Mercaptoethanol",
        installed_ligand="Thioether (blocked)",
        functional_mode="quencher",
        reaction_type="blocking",
        chemistry_class="vs_thiol",
        target_acs=ACSSiteType.VINYL_SULFONE,
        product_acs=None,
        k_forward=5e-3,
        E_a=25000.0,
        stoichiometry=1.0,
        hydrolysis_rate=0.0,
        ph_optimum=6.5,
        ph_min=5.0, ph_max=8.0,
        temperature_default=298.15,
        time_default=3600.0,
        confidence_tier="semi_quantitative",
        hazard_class="toxic_malodorous",
        notes="Fast thiol-VS Michael addition; handle in fume hood",
    ),

    # ── 13. Sodium borohydride quench ─────────────────────────────────
    "nabh4_quench": ReagentProfile(
        name="Sodium borohydride (aldehyde quench)",
        cas="16940-66-2",
        reagent_identity="Sodium borohydride (NaBH4)",
        installed_ligand="Alcohol -CH2OH (blocked)",
        functional_mode="quencher",
        reaction_type="blocking",
        chemistry_class="reduction",
        target_acs=ACSSiteType.ALDEHYDE,
        product_acs=None,
        k_forward=1e-1,
        E_a=15000.0,
        stoichiometry=1.0,
        hydrolysis_rate=0.0,
        ph_optimum=7.5,
        ph_min=6.0, ph_max=9.0,
        temperature_default=298.15,
        time_default=1800.0,
        confidence_tier="semi_quantitative",
        hazard_class="flammable_corrosive",
        notes="Also reduces Schiff base (imine) linkages to stable secondary amines",
    ),

    # ── 14. Acetic anhydride quench ───────────────────────────────────
    "acetic_anhydride_quench": ReagentProfile(
        name="Acetic anhydride (amine quench)",
        cas="108-24-7",
        reagent_identity="Acetic anhydride",
        installed_ligand="Acetamide (blocked)",
        functional_mode="quencher",
        reaction_type="blocking",
        chemistry_class="acetylation",
        target_acs=ACSSiteType.AMINE_PRIMARY,
        product_acs=None,
        k_forward=5e-3,
        E_a=25000.0,
        stoichiometry=1.0,
        hydrolysis_rate=0.0,
        ph_optimum=7.5,
        ph_min=6.0, ph_max=9.0,
        temperature_default=298.15,
        time_default=3600.0,
        confidence_tier="semi_quantitative",
        hazard_class="corrosive_flammable",
        notes="Caps free amines to reduce non-specific binding",
    ),
}
