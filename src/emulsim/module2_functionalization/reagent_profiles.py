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

    # ── Spacer arm support (Phase 1 multiplier model) ──
    spacer_key: str = ""                    # Reference to spacer profile key (empty = direct)
    spacer_length_angstrom: float = 0.0     # Spacer length [angstrom]
    spacer_activity_multiplier: float = 1.0 # Multiplier on activity_retention [>=1.0]

    # ── Charge type for IEX (audit F14) ──
    charge_type: str = ""                   # "anion", "cation", ""

    # ── IMAC metal state (audit F5) ──
    metal_ion: str = ""                     # "Ni2+", "Co2+", "Cu2+", "Zn2+", ""
    metal_loaded_fraction: float = 1.0      # [0,1] assumed metal loading

    # ── Binding model hint for M3 (audit F15) ──
    binding_model_hint: str = ""            # "charge_exchange", "metal_chelation",
                                            # "salt_promoted", "fc_affinity",
                                            # "gst_glutathione", "near_irreversible",
                                            # "mixed_mode", ""

    # ── v5.8 fields (audit F2/F3/F4/F5/F10/F15) ──
    profile_role: str = "final_ligand"      # "native", "activated", "spacer_intermediate",
                                             # "heterobifunctional_intermediate", "final_ligand",
                                             # "spacer_metadata", "quencher"
    m3_support_level: str = "mapped_estimated"  # "mapped_quantitative", "mapped_estimated",
                                                 # "not_mapped", "requires_user_calibration"
    distal_group_yield: float = 1.0         # Fraction consumed sites producing distal group [0,1]
    maleimide_decay_rate: float = 0.0       # [1/s] first-order hydrolysis of immobilized maleimide
    buffer_incompatibilities: str = ""      # Comma-separated: "Tris,glycine,DTT"
    requires_reduced_thiol: bool = False    # True for maleimide-thiol coupling
    thiol_accessibility_fraction: float = 1.0  # Fraction of protein Cys accessible [0,1]

    # ── v5.9.1-5.9.4 fields ──
    metal_association_constant: float = 0.0  # [m3/mol] metal-chelator association constant
    reduction_efficiency: float = 0.95       # [0,1] protein pretreatment efficiency
    regulatory_limit_ppm: float = 0.0        # [ppm] for washing compliance check
    pKa_nucleophile: float = 0.0             # pH scaling (0 = disabled)

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
        charge_type="anion",
        binding_model_hint="charge_exchange",
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
        metal_ion="Ni2+",
        metal_loaded_fraction=1.0,
        binding_model_hint="metal_chelation",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from general epoxide-amine kinetics",
        notes="IMAC chelator; assumes fully Ni2+-loaded; no metal leaching modeled",
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
        binding_model_hint="salt_promoted",
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
        charge_type="cation",
        binding_model_hint="charge_exchange",
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
        binding_model_hint="fc_affinity",
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
        binding_model_hint="fc_affinity",
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

    # ═══════════════════════════════════════════════════════════════════
    # v5.7 Expansion — 6 new ligand coupling profiles (WN-2)
    # ═══════════════════════════════════════════════════════════════════

    "q_coupling": ReagentProfile(
        name="Q (strong anion exchange)",
        cas="3033-77-0",
        reagent_identity="Glycidyltrimethylammonium chloride",
        installed_ligand="Q",
        functional_mode="iex_ligand",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=6e-5,
        E_a=50000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-5,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        time_default=14400.0,
        ligand_mw=152.0,
        ligand_r_h=0.4e-9,
        charge_type="anion",
        binding_model_hint="charge_exchange",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from epoxide-amine kinetics",
        notes="Strong anion exchanger; permanent positive charge at all pH",
    ),

    "cm_coupling": ReagentProfile(
        name="CM (weak cation exchange)",
        cas="79-11-8",
        reagent_identity="Chloroacetic acid + amino spacer",
        installed_ligand="CM-like carboxymethyl",
        functional_mode="iex_ligand",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=3e-5,
        E_a=45000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-5,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        time_default=21600.0,
        ligand_mw=94.0,
        ligand_r_h=0.4e-9,
        charge_type="cation",
        binding_model_hint="charge_exchange",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated; CM-like via amino-carboxyl ligand coupling",
        notes="Weak cation exchanger; charged above pH ~4 (carboxyl pKa)",
    ),

    "nta_coupling": ReagentProfile(
        name="NTA (IMAC chelator, His-tag)",
        cas="139-13-9",
        reagent_identity="Nitrilotriacetic acid",
        installed_ligand="NTA",
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
        time_default=21600.0,
        ligand_mw=191.0,
        ligand_r_h=0.5e-9,
        metal_ion="Ni2+",
        metal_loaded_fraction=1.0,
        binding_model_hint="metal_chelation",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from epoxide-amine kinetics",
        notes="Tetradentate IMAC chelator; assumes fully Ni2+-loaded; higher specificity than IDA",
    ),

    "butyl_coupling": ReagentProfile(
        name="Butyl (HIC, mild)",
        cas="109-73-9",
        reagent_identity="n-Butylamine",
        installed_ligand="Butyl",
        functional_mode="hic_ligand",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=4e-5,
        E_a=48000.0,
        stoichiometry=1.0,
        hydrolysis_rate=5e-6,
        ph_optimum=9.5,
        ph_min=8.0, ph_max=11.0,
        temperature_default=298.15,
        time_default=14400.0,
        ligand_mw=73.0,
        ligand_r_h=0.3e-9,
        binding_model_hint="salt_promoted",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from alkylamine-epoxide kinetics",
        notes="Mild HIC ligand; lower hydrophobicity than Phenyl; gentler elution conditions",
    ),

    "glutathione_coupling": ReagentProfile(
        name="Glutathione (GST-tag affinity)",
        cas="70-18-8",
        reagent_identity="Glutathione reduced (GSH)",
        installed_ligand="Glutathione",
        functional_mode="gst_affinity",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=3e-5,
        E_a=45000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-5,
        ph_optimum=9.5,
        ph_min=8.0, ph_max=11.0,
        temperature_default=298.15,
        time_default=14400.0,
        ligand_mw=307.0,
        ligand_r_h=0.5e-9,
        activity_retention=0.80,
        activity_retention_uncertainty=0.15,
        spacer_key="aha_spacer",
        spacer_activity_multiplier=1.15,
        binding_model_hint="gst_glutathione",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated; orientation-sensitive coupling",
        notes="Orientation-sensitive; activity_retention=0.80 reflects random coupling penalty",
    ),

    "heparin_coupling": ReagentProfile(
        name="Heparin (affinity + IEX)",
        cas="9005-49-6",
        reagent_identity="Heparin sodium (porcine intestinal)",
        installed_ligand="Heparin",
        functional_mode="heparin_affinity",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=1e-5,
        E_a=40000.0,
        stoichiometry=1.0,
        hydrolysis_rate=5e-6,
        ph_optimum=9.5,
        ph_min=8.0, ph_max=11.0,
        temperature_default=298.15,
        time_default=28800.0,
        ligand_mw=14000.0,
        ligand_r_h=3.0e-9,
        is_macromolecule=True,
        spacer_key="dadpa_spacer",
        spacer_activity_multiplier=1.22,
        binding_model_hint="mixed_mode",
        confidence_tier="ranking_only",
        calibration_source="Estimated; macromolecular polysaccharide",
        notes="Mixed affinity + cation exchange; q_max highly target-dependent; macromolecule",
    ),

    # ═══════════════════════════════════════════════════════════════════
    # v5.7 Expansion — 2 new protein coupling profiles (WN-3)
    # ═══════════════════════════════════════════════════════════════════

    "protein_ag_coupling": ReagentProfile(
        name="Protein A/G Fusion (broadest IgG)",
        cas="N/A",
        reagent_identity="Recombinant Protein A/G",
        installed_ligand="Protein A/G",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=5e-7,
        E_a=25000.0,
        stoichiometry=1.0,
        ph_optimum=9.0,
        ph_min=7.5, ph_max=10.0,
        temperature_default=277.15,
        temperature_min=273.15, temperature_max=283.15,
        time_default=57600.0,
        ligand_mw=51000.0,
        ligand_r_h=2.8e-9,
        is_macromolecule=True,
        activity_retention=0.55,
        activity_retention_uncertainty=0.15,
        max_surface_density=2e-8,
        spacer_key="dadpa_spacer",
        spacer_activity_multiplier=1.22,
        binding_model_hint="fc_affinity",
        confidence_tier="ranking_only",
        calibration_source="Estimated; broadest IgG coverage via A+G domains",
        notes="Fusion of Protein A + G binding domains; ranking_only unless calibrated",
    ),

    "streptavidin_coupling": ReagentProfile(
        name="Streptavidin (biotin-tag)",
        cas="9013-20-1",
        reagent_identity="Streptavidin (S. avidinii)",
        installed_ligand="Streptavidin",
        functional_mode="biotin_affinity",
        reaction_type="protein_coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=4e-7,
        E_a=25000.0,
        stoichiometry=1.0,
        ph_optimum=9.0,
        ph_min=7.5, ph_max=10.0,
        temperature_default=277.15,
        temperature_min=273.15, temperature_max=283.15,
        time_default=57600.0,
        ligand_mw=53000.0,
        ligand_r_h=2.8e-9,
        is_macromolecule=True,
        activity_retention=0.70,
        activity_retention_uncertainty=0.10,
        max_surface_density=3e-8,
        spacer_key="dadpa_spacer",
        spacer_activity_multiplier=1.22,
        binding_model_hint="near_irreversible",
        confidence_tier="ranking_only",
        calibration_source="Estimated; Kd ~10^-15 M streptavidin-biotin",
        notes="Near-irreversible binding; effective stoichiometry capped at 2.5 (not theoretical 4)",
    ),

    # ═══════════════════════════════════════════════════════════════════
    # v5.7 Expansion — 3 spacer arm metadata profiles (WN-3)
    # These are NOT executable as modification steps. They serve as
    # lookup records for spacer_key on coupling profiles.
    # ═══════════════════════════════════════════════════════════════════

    "dadpa_spacer": ReagentProfile(
        name="DADPA spacer (13 A, EAH-standard)",
        cas="56-18-8",
        reagent_identity="Diaminodipropylamine",
        installed_ligand="DADPA spacer arm",
        functional_mode="spacer",
        reaction_type="spacer",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=0.0,
        E_a=0.0,
        stoichiometry=1.0,
        spacer_length_angstrom=13.0,
        spacer_activity_multiplier=1.22,
        confidence_tier="semi_quantitative",
        calibration_source="EAH-Sepharose standard; multiplier from Protein A literature",
        hazard_class="irritant",
        notes="9-atom amine spacer; industry standard for protein ligand immobilization",
    ),

    "aha_spacer": ReagentProfile(
        name="AHA spacer (10 A, NHS-standard)",
        cas="60-32-2",
        reagent_identity="6-Aminohexanoic acid",
        installed_ligand="AHA spacer arm",
        functional_mode="spacer",
        reaction_type="spacer",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=0.0,
        E_a=0.0,
        stoichiometry=1.0,
        spacer_length_angstrom=10.0,
        spacer_activity_multiplier=1.15,
        confidence_tier="semi_quantitative",
        calibration_source="NHS-Sepharose HP standard; multiplier from affinity literature",
        hazard_class="low_hazard",
        notes="7-atom acid spacer; provides -COOH distal group (EDC/NHS path not modeled in Phase 1)",
    ),

    "dah_spacer": ReagentProfile(
        name="DAH spacer (9 A, AH-standard)",
        cas="124-09-4",
        reagent_identity="1,6-Diaminohexane",
        installed_ligand="DAH spacer arm",
        functional_mode="spacer",
        reaction_type="spacer",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=0.0,
        E_a=0.0,
        stoichiometry=1.0,
        spacer_length_angstrom=9.0,
        spacer_activity_multiplier=1.08,
        confidence_tier="semi_quantitative",
        calibration_source="AH-Sepharose 4B standard; shorter than DADPA",
        hazard_class="irritant",
        notes="6-atom diamine spacer; simpler than DADPA but more hydrophobic",
    ),

    # ═══════════════════════════════════════════════════════════════════
    # v5.8 Phase 1 — 6 deferred profiles (WN-0)
    # ═══════════════════════════════════════════════════════════════════

    "protein_l_coupling": ReagentProfile(
        name="Protein L (kappa light chain)",
        cas="N/A",
        reagent_identity="Recombinant Protein L (Peptostreptococcus magnus)",
        installed_ligand="Protein L",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=5e-7,
        E_a=25000.0,
        stoichiometry=1.0,
        ph_optimum=9.0,
        ph_min=7.5, ph_max=10.0,
        temperature_default=277.15,
        temperature_min=273.15, temperature_max=283.15,
        time_default=57600.0,
        ligand_mw=36000.0,
        ligand_r_h=2.3e-9,
        is_macromolecule=True,
        activity_retention=0.55,
        activity_retention_uncertainty=0.15,
        max_surface_density=3e-8,
        spacer_key="dadpa_spacer",
        spacer_activity_multiplier=1.20,
        binding_model_hint="kappa_light_chain_affinity",
        profile_role="final_ligand",
        m3_support_level="requires_user_calibration",
        confidence_tier="ranking_only",
        calibration_source="Estimated; Fab/scFv purification via kappa light chain",
        notes="Binds kappa light chains; useful for Fab, scFv, single-domain antibodies",
    ),

    "con_a_coupling": ReagentProfile(
        name="Concanavalin A (lectin, mannose/glucose)",
        cas="11028-71-0",
        reagent_identity="Concanavalin A (Canavalia ensiformis)",
        installed_ligand="Con A",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=2e-7,
        E_a=25000.0,
        stoichiometry=1.0,
        ph_optimum=8.5,
        ph_min=7.0, ph_max=9.5,
        temperature_default=277.15,
        temperature_min=273.15, temperature_max=283.15,
        time_default=57600.0,
        ligand_mw=104000.0,
        ligand_r_h=4.0e-9,
        is_macromolecule=True,
        activity_retention=0.40,
        activity_retention_uncertainty=0.20,
        max_surface_density=1e-8,
        spacer_key="peg600_spacer",
        spacer_activity_multiplier=1.35,
        binding_model_hint="lectin_mannose_glucose",
        profile_role="final_ligand",
        m3_support_level="requires_user_calibration",
        confidence_tier="ranking_only",
        calibration_source="Estimated; lectin affinity literature",
        notes="Tetramer at pH>7; requires Ca2+/Mn2+ cofactors; sugar-competition elution",
    ),

    "octyl_coupling": ReagentProfile(
        name="Octyl (HIC, strong)",
        cas="111-86-4",
        reagent_identity="n-Octylamine",
        installed_ligand="Octyl",
        functional_mode="hic_ligand",
        reaction_type="coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=4e-5,
        E_a=48000.0,
        stoichiometry=1.0,
        hydrolysis_rate=5e-6,
        ph_optimum=9.5,
        ph_min=8.0, ph_max=11.0,
        temperature_default=298.15,
        time_default=14400.0,
        ligand_mw=129.0,
        ligand_r_h=0.4e-9,
        binding_model_hint="salt_promoted",
        profile_role="final_ligand",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated from alkylamine-epoxide kinetics",
        notes="Very hydrophobic HIC; risk of irreversible binding at high density",
    ),

    "wga_coupling": ReagentProfile(
        name="WGA (wheat germ agglutinin, GlcNAc/sialic)",
        cas="9001-31-2",
        reagent_identity="Wheat Germ Agglutinin",
        installed_ligand="WGA",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=4e-7,
        E_a=25000.0,
        stoichiometry=1.0,
        ph_optimum=8.5,
        ph_min=7.0, ph_max=9.5,
        temperature_default=277.15,
        temperature_min=273.15, temperature_max=283.15,
        time_default=57600.0,
        ligand_mw=36000.0,
        ligand_r_h=2.3e-9,
        is_macromolecule=True,
        activity_retention=0.50,
        activity_retention_uncertainty=0.15,
        max_surface_density=3e-8,
        spacer_key="dadpa_spacer",
        spacer_activity_multiplier=1.20,
        binding_model_hint="lectin_glcnac_sialic",
        profile_role="final_ligand",
        m3_support_level="requires_user_calibration",
        confidence_tier="ranking_only",
        calibration_source="Estimated; lectin affinity literature",
        notes="Dimer 2x18kDa; binds GlcNAc and sialic acid; sugar-competition elution",
    ),

    "peg600_spacer": ReagentProfile(
        name="PEG-diamine Mn 600 spacer (35 A)",
        cas="929-59-9",
        reagent_identity="Bis(aminopropyl) PEG Mn 600",
        installed_ligand="PEG-diamine spacer arm",
        functional_mode="spacer",
        reaction_type="spacer",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=None,
        k_forward=0.0,
        E_a=0.0,
        stoichiometry=1.0,
        spacer_length_angstrom=35.0,
        spacer_activity_multiplier=1.35,
        ligand_mw=600.0,
        profile_role="spacer_metadata",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated; PEG-diamine literature",
        notes="Long hydrophilic PEG spacer; critical for large proteins (>50 kDa)",
    ),

    "bdge_activation": ReagentProfile(
        name="BDGE activation (18 A spacer epoxide)",
        cas="2425-79-8",
        reagent_identity="1,4-Butanediol diglycidyl ether",
        installed_ligand="Epoxide (long-arm)",
        functional_mode="activator",
        reaction_type="activation",
        chemistry_class="epoxide_amine",
        target_acs=ACSSiteType.HYDROXYL,
        product_acs=ACSSiteType.EPOXIDE,
        k_forward=1.2e-5,
        E_a=60000.0,
        stoichiometry=1.0,
        hydrolysis_rate=5e-5,
        ph_optimum=11.5,
        ph_min=10.0, ph_max=13.0,
        temperature_default=298.15,
        time_default=14400.0,
        spacer_length_angstrom=18.0,
        profile_role="activated",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated; bis-epoxide activation literature",
        hazard_class="irritant",
        notes="Bis-epoxide; creates long-arm epoxide (18A vs ECH ~5A)",
    ),

    # ═══════════════════════════════════════════════════════════════════
    # v5.8 Phase 2 — SPACER_ARM executable profiles (WN-5)
    # ═══════════════════════════════════════════════════════════════════

    # ── Diamine spacer-arm profiles (EPOXIDE → AMINE_DISTAL) ──

    "eda_spacer_arm": ReagentProfile(
        name="EDA spacer arm (3 A, shortest)",
        cas="107-15-3",
        reagent_identity="Ethylenediamine",
        installed_ligand="EDA spacer (distal -NH2)",
        functional_mode="spacer",
        reaction_type="spacer_arm",
        chemistry_class="epoxide_amine_spacer",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=ACSSiteType.AMINE_DISTAL,
        k_forward=8e-5,
        E_a=45000.0,
        stoichiometry=1.0,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        time_default=14400.0,
        spacer_length_angstrom=3.0,
        distal_group_yield=0.60,
        profile_role="spacer_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated; high bridging risk for short diamine",
        notes="Shortest diamine; 60% distal yield (40% bridging)",
    ),

    "dadpa_spacer_arm": ReagentProfile(
        name="DADPA spacer arm (13 A)",
        cas="56-18-8",
        reagent_identity="Diaminodipropylamine",
        installed_ligand="DADPA spacer (distal -NH2)",
        functional_mode="spacer",
        reaction_type="spacer_arm",
        chemistry_class="epoxide_amine_spacer",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=ACSSiteType.AMINE_DISTAL,
        k_forward=5e-5,
        E_a=45000.0,
        stoichiometry=1.0,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        time_default=14400.0,
        spacer_length_angstrom=13.0,
        distal_group_yield=0.80,
        profile_role="spacer_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="EAH-Sepharose standard; internal amine reduces bridging",
        notes="Industry standard protein spacer; 80% distal yield",
    ),

    "dah_spacer_arm": ReagentProfile(
        name="DAH spacer arm (9 A)",
        cas="124-09-4",
        reagent_identity="1,6-Diaminohexane",
        installed_ligand="DAH spacer (distal -NH2)",
        functional_mode="spacer",
        reaction_type="spacer_arm",
        chemistry_class="epoxide_amine_spacer",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=ACSSiteType.AMINE_DISTAL,
        k_forward=6e-5,
        E_a=45000.0,
        stoichiometry=1.0,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        time_default=14400.0,
        spacer_length_angstrom=9.0,
        distal_group_yield=0.70,
        profile_role="spacer_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="AH-Sepharose standard",
        hazard_class="irritant",
        notes="Simple diamine spacer; 70% distal yield (more hydrophobic than DADPA)",
    ),

    "peg600_spacer_arm": ReagentProfile(
        name="PEG-diamine Mn 600 spacer arm (35 A)",
        cas="929-59-9",
        reagent_identity="Bis(aminopropyl) PEG Mn 600",
        installed_ligand="PEG spacer (distal -NH2)",
        functional_mode="spacer",
        reaction_type="spacer_arm",
        chemistry_class="epoxide_amine_spacer",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=ACSSiteType.AMINE_DISTAL,
        k_forward=2e-5,
        E_a=45000.0,
        stoichiometry=1.0,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        time_default=21600.0,
        spacer_length_angstrom=35.0,
        distal_group_yield=0.90,
        ligand_mw=600.0,
        profile_role="spacer_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="PEG-diamine literature; PEG flexibility favors monoattachment",
        notes="Long hydrophilic PEG spacer; 90% distal yield; critical for large proteins",
    ),

    # ── SM(PEG)n heterobifunctional crosslinkers (AMINE_DISTAL → MALEIMIDE) ──

    "sm_peg2": ReagentProfile(
        name="SM(PEG)2 (NHS-PEG2-Maleimide, 18 A)",
        cas="1334179-85-1",
        reagent_identity="Succinimidyl-[(N-maleimidopropionamido)-diethyleneglycol] ester",
        installed_ligand="Maleimide (PEG2 arm)",
        functional_mode="heterobifunctional_crosslinker",
        reaction_type="heterobifunctional",
        chemistry_class="nhs_amine",
        target_acs=ACSSiteType.AMINE_DISTAL,
        product_acs=ACSSiteType.MALEIMIDE,
        k_forward=1e-3,
        E_a=35000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-3,
        ph_optimum=7.4,
        ph_min=7.0, ph_max=7.5,
        temperature_default=298.15,
        time_default=1800.0,
        spacer_length_angstrom=18.0,
        ligand_mw=425.0,
        distal_group_yield=0.85,
        maleimide_decay_rate=1e-5,
        buffer_incompatibilities="Tris,glycine,primary_amine_buffers",
        profile_role="heterobifunctional_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Thermo Fisher SM(PEG)n protocols; NHS hydrolysis from Hermanson",
        notes="NHS half-life ~10 min at pH 7.4; react quickly after dissolving",
    ),

    "sm_peg4": ReagentProfile(
        name="SM(PEG)4 (NHS-PEG4-Maleimide, 32 A)",
        cas="1229578-42-6",
        reagent_identity="Succinimidyl-[(N-maleimidopropionamido)-tetraethyleneglycol] ester",
        installed_ligand="Maleimide (PEG4 arm)",
        functional_mode="heterobifunctional_crosslinker",
        reaction_type="heterobifunctional",
        chemistry_class="nhs_amine",
        target_acs=ACSSiteType.AMINE_DISTAL,
        product_acs=ACSSiteType.MALEIMIDE,
        k_forward=1e-3,
        E_a=35000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-3,
        ph_optimum=7.4,
        ph_min=7.0, ph_max=7.5,
        temperature_default=298.15,
        time_default=1800.0,
        spacer_length_angstrom=32.0,
        ligand_mw=513.0,
        distal_group_yield=0.85,
        maleimide_decay_rate=1e-5,
        buffer_incompatibilities="Tris,glycine,primary_amine_buffers",
        profile_role="heterobifunctional_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Thermo Fisher SM(PEG)n protocols",
        notes="Most commonly used SM(PEG)n variant; good balance of length and solubility",
    ),

    "sm_peg12": ReagentProfile(
        name="SM(PEG)12 (NHS-PEG12-Maleimide, 60 A)",
        cas="1334179-86-2",
        reagent_identity="Succinimidyl-[(N-maleimidopropionamido)-dodecaethyleneglycol] ester",
        installed_ligand="Maleimide (PEG12 arm)",
        functional_mode="heterobifunctional_crosslinker",
        reaction_type="heterobifunctional",
        chemistry_class="nhs_amine",
        target_acs=ACSSiteType.AMINE_DISTAL,
        product_acs=ACSSiteType.MALEIMIDE,
        k_forward=1e-3,
        E_a=35000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-3,
        ph_optimum=7.4,
        ph_min=7.0, ph_max=7.5,
        temperature_default=298.15,
        time_default=1800.0,
        spacer_length_angstrom=60.0,
        ligand_mw=865.0,
        distal_group_yield=0.85,
        maleimide_decay_rate=1e-5,
        buffer_incompatibilities="Tris,glycine,primary_amine_buffers",
        profile_role="heterobifunctional_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Thermo Fisher SM(PEG)n protocols",
        notes="Long PEG arm for maximum protein conformational freedom",
    ),

    "sm_peg24": ReagentProfile(
        name="SM(PEG)24 (NHS-PEG24-Maleimide, 95 A)",
        cas="1334179-87-3",
        reagent_identity="Succinimidyl-[(N-maleimidopropionamido)-tetracosaethyleneglycol] ester",
        installed_ligand="Maleimide (PEG24 arm)",
        functional_mode="heterobifunctional_crosslinker",
        reaction_type="heterobifunctional",
        chemistry_class="nhs_amine",
        target_acs=ACSSiteType.AMINE_DISTAL,
        product_acs=ACSSiteType.MALEIMIDE,
        k_forward=1e-3,
        E_a=35000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-3,
        ph_optimum=7.4,
        ph_min=7.0, ph_max=7.5,
        temperature_default=298.15,
        time_default=1800.0,
        spacer_length_angstrom=95.0,
        ligand_mw=1393.0,
        distal_group_yield=0.85,
        maleimide_decay_rate=1e-5,
        buffer_incompatibilities="Tris,glycine,primary_amine_buffers",
        profile_role="heterobifunctional_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Thermo Fisher SM(PEG)n protocols",
        notes="Longest SM(PEG)n variant; ~95 A; may introduce excess flexibility",
    ),

    # ── Protein-Cys coupling profiles (MALEIMIDE → thioether) ──

    "protein_a_cys_coupling": ReagentProfile(
        name="Protein A-Cys (oriented, maleimide-thiol)",
        cas="91932-65-9",
        reagent_identity="Recombinant Protein A with engineered C-terminal Cys",
        installed_ligand="Protein A (oriented)",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="maleimide_thiol",
        target_acs=ACSSiteType.MALEIMIDE,
        product_acs=None,
        k_forward=5e-2,
        E_a=20000.0,
        stoichiometry=1.0,
        ph_optimum=7.0,
        ph_min=6.5, ph_max=7.5,
        temperature_default=298.15,
        time_default=7200.0,
        ligand_mw=42000.0,
        ligand_r_h=2.5e-9,
        is_macromolecule=True,
        activity_retention=0.80,
        activity_retention_uncertainty=0.10,
        max_surface_density=2e-8,
        requires_reduced_thiol=True,
        thiol_accessibility_fraction=0.90,
        binding_model_hint="fc_affinity",
        buffer_incompatibilities="DTT,beta-mercaptoethanol,TCEP,free_thiols",
        profile_role="final_ligand",
        m3_support_level="mapped_estimated",
        confidence_tier="ranking_only",
        calibration_source="Site-specific conjugation literature; improved activity over random",
        notes="Oriented immobilization via C-terminal Cys; ~80% activity (vs 60% random)",
    ),

    "protein_g_cys_coupling": ReagentProfile(
        name="Protein G-Cys (oriented, maleimide-thiol)",
        cas="122441-07-8",
        reagent_identity="Recombinant Protein G with engineered C-terminal Cys",
        installed_ligand="Protein G (oriented)",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="maleimide_thiol",
        target_acs=ACSSiteType.MALEIMIDE,
        product_acs=None,
        k_forward=5e-2,
        E_a=20000.0,
        stoichiometry=1.0,
        ph_optimum=7.0,
        ph_min=6.5, ph_max=7.5,
        temperature_default=298.15,
        time_default=7200.0,
        ligand_mw=22000.0,
        ligand_r_h=2.0e-9,
        is_macromolecule=True,
        activity_retention=0.85,
        activity_retention_uncertainty=0.10,
        max_surface_density=3e-8,
        requires_reduced_thiol=True,
        thiol_accessibility_fraction=0.90,
        binding_model_hint="fc_affinity",
        buffer_incompatibilities="DTT,beta-mercaptoethanol,TCEP,free_thiols",
        profile_role="final_ligand",
        m3_support_level="mapped_estimated",
        confidence_tier="ranking_only",
        calibration_source="Site-specific conjugation literature",
        notes="Oriented Protein G via Cys; broader subclass than Protein A",
    ),

    "generic_cys_protein_coupling": ReagentProfile(
        name="Generic Cys-protein (oriented, maleimide-thiol)",
        cas="N/A",
        reagent_identity="User-supplied Cys-tagged protein",
        installed_ligand="Cys-protein (oriented)",
        functional_mode="affinity_ligand",
        reaction_type="protein_coupling",
        chemistry_class="maleimide_thiol",
        target_acs=ACSSiteType.MALEIMIDE,
        product_acs=None,
        k_forward=5e-2,
        E_a=20000.0,
        stoichiometry=1.0,
        ph_optimum=7.0,
        ph_min=6.5, ph_max=7.5,
        temperature_default=298.15,
        time_default=7200.0,
        ligand_mw=50000.0,
        ligand_r_h=2.5e-9,
        is_macromolecule=True,
        activity_retention=0.70,
        activity_retention_uncertainty=0.15,
        max_surface_density=2e-8,
        requires_reduced_thiol=True,
        thiol_accessibility_fraction=0.80,
        binding_model_hint="fc_affinity",
        buffer_incompatibilities="DTT,beta-mercaptoethanol,TCEP,free_thiols",
        profile_role="final_ligand",
        m3_support_level="requires_user_calibration",
        confidence_tier="ranking_only",
        calibration_source="Generic; user must provide target-specific parameters",
        notes="Generic maleimide-thiol coupling; user must specify protein MW and activity",
    ),

    # ═══════════════════════════════════════════════════════════════════
    # v5.9.1 — IMAC Metal Charging (5 profiles)
    # ═══════════════════════════════════════════════════════════════════

    "nickel_charging": ReagentProfile(
        name="Nickel(II) charging (Ni-NTA/IDA)",
        cas="7786-81-4", reagent_identity="Nickel(II) sulfate",
        installed_ligand="Ni2+-loaded chelator",
        functional_mode="metal_charging", reaction_type="metal_charging",
        chemistry_class="metal_chelation",
        target_acs=ACSSiteType.EPOXIDE, product_acs=None,
        k_forward=0.0, E_a=0.0, stoichiometry=1.0,
        metal_ion="Ni2+", metal_association_constant=3e11,
        ph_optimum=7.0, ph_min=5.0, ph_max=8.0,
        temperature_default=298.15, time_default=1800.0,
        profile_role="native", m3_support_level="mapped_estimated",
        confidence_tier="semi_quantitative",
        calibration_source="log K(NTA-Ni)=11.5; Martell & Smith Critical Stability Constants",
        notes="Standard Ni2+ charging for His-tag IMAC; 50 mM NiSO4 typical",
    ),
    "cobalt_charging": ReagentProfile(
        name="Cobalt(II) charging (Co-NTA/IDA)",
        cas="10026-24-1", reagent_identity="Cobalt(II) chloride",
        installed_ligand="Co2+-loaded chelator",
        functional_mode="metal_charging", reaction_type="metal_charging",
        chemistry_class="metal_chelation",
        target_acs=ACSSiteType.EPOXIDE, product_acs=None,
        k_forward=0.0, E_a=0.0, stoichiometry=1.0,
        metal_ion="Co2+", metal_association_constant=1e10,
        ph_optimum=7.0, ph_min=5.0, ph_max=8.0,
        temperature_default=298.15, time_default=1800.0,
        profile_role="native", m3_support_level="mapped_estimated",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated; Co2+ has higher specificity than Ni2+",
        notes="Higher specificity, lower capacity than Ni2+",
    ),
    "copper_charging": ReagentProfile(
        name="Copper(II) charging (Cu-IDA)",
        cas="7758-99-8", reagent_identity="Copper(II) sulfate",
        installed_ligand="Cu2+-loaded chelator",
        functional_mode="metal_charging", reaction_type="metal_charging",
        chemistry_class="metal_chelation",
        target_acs=ACSSiteType.EPOXIDE, product_acs=None,
        k_forward=0.0, E_a=0.0, stoichiometry=1.0,
        metal_ion="Cu2+", metal_association_constant=1e13,
        ph_optimum=7.0, ph_min=4.0, ph_max=8.0,
        temperature_default=298.15, time_default=1800.0,
        profile_role="native", m3_support_level="mapped_estimated",
        confidence_tier="semi_quantitative",
        calibration_source="log K(IDA-Cu)=10.6; strongest divalent binding",
        notes="Highest affinity but lowest specificity; more non-specific binding",
    ),
    "zinc_charging": ReagentProfile(
        name="Zinc(II) charging (Zn-NTA/IDA)",
        cas="7446-20-0", reagent_identity="Zinc(II) sulfate",
        installed_ligand="Zn2+-loaded chelator",
        functional_mode="metal_charging", reaction_type="metal_charging",
        chemistry_class="metal_chelation",
        target_acs=ACSSiteType.EPOXIDE, product_acs=None,
        k_forward=0.0, E_a=0.0, stoichiometry=1.0,
        metal_ion="Zn2+", metal_association_constant=1e10,
        ph_optimum=7.0, ph_min=5.0, ph_max=8.0,
        temperature_default=298.15, time_default=1800.0,
        profile_role="native", m3_support_level="mapped_estimated",
        confidence_tier="semi_quantitative",
        calibration_source="Estimated; Zn2+ used for some phosphopeptide IMAC",
    ),
    "edta_stripping": ReagentProfile(
        name="EDTA metal stripping",
        cas="60-00-4", reagent_identity="EDTA disodium salt",
        installed_ligand="Metal-free chelator",
        functional_mode="metal_charging", reaction_type="metal_stripping",
        chemistry_class="metal_chelation",
        target_acs=ACSSiteType.EPOXIDE, product_acs=None,
        k_forward=0.0, E_a=0.0, stoichiometry=1.0,
        ph_optimum=7.5, ph_min=6.0, ph_max=9.0,
        temperature_default=298.15, time_default=1800.0,
        profile_role="native", m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Standard EDTA stripping; 50 mM strips >99% Ni/Co/Cu",
        notes="Strips loaded metal for regeneration or metal switching",
    ),

    # ═══════════════════════════════════════════════════════════════════
    # v5.9.2 — Protein Pretreatment (2 profiles)
    # ═══════════════════════════════════════════════════════════════════

    "tcep_reduction": ReagentProfile(
        name="TCEP disulfide reduction",
        cas="51805-45-9", reagent_identity="TCEP-HCl",
        installed_ligand="Reduced protein (free -SH)",
        functional_mode="protein_pretreatment", reaction_type="protein_pretreatment",
        chemistry_class="reduction",
        target_acs=ACSSiteType.MALEIMIDE, product_acs=None,
        k_forward=0.01, E_a=20000.0, stoichiometry=1.0,
        ph_optimum=7.0, ph_min=6.0, ph_max=8.0,
        temperature_default=298.15, time_default=1800.0,
        reduction_efficiency=0.95,
        activity_retention=0.95,
        profile_role="native", m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Hermanson Bioconjugate Techniques; TCEP is maleimide-compatible",
        notes="Preferred over DTT; does not interfere with maleimide coupling",
    ),
    "dtt_reduction": ReagentProfile(
        name="DTT disulfide reduction",
        cas="3483-12-3", reagent_identity="Dithiothreitol",
        installed_ligand="Reduced protein (free -SH)",
        functional_mode="protein_pretreatment", reaction_type="protein_pretreatment",
        chemistry_class="reduction",
        target_acs=ACSSiteType.MALEIMIDE, product_acs=None,
        k_forward=0.005, E_a=20000.0, stoichiometry=1.0,
        ph_optimum=7.5, ph_min=6.5, ph_max=8.5,
        temperature_default=298.15, time_default=1800.0,
        reduction_efficiency=0.90,
        activity_retention=0.90,
        buffer_incompatibilities="maleimide,free_thiols_in_coupling_buffer",
        profile_role="native", m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Standard DTT reduction; must remove excess before maleimide",
        hazard_class="irritant",
        notes="Must remove DTT before maleimide coupling (desalt or dialysis)",
    ),

    # ═══════════════════════════════════════════════════════════════════
    # v5.9.3 — EDC/NHS Chemistry (2 profiles)
    # ═══════════════════════════════════════════════════════════════════

    "aha_carboxyl_spacer_arm": ReagentProfile(
        name="AHA spacer arm (EPOXIDE -> CARBOXYL_DISTAL)",
        cas="60-32-2",
        reagent_identity="6-Aminohexanoic acid",
        installed_ligand="AHA spacer (distal -COOH)",
        functional_mode="spacer",
        reaction_type="spacer_arm",
        chemistry_class="epoxide_amine_spacer",
        target_acs=ACSSiteType.EPOXIDE,
        product_acs=ACSSiteType.CARBOXYL_DISTAL,
        k_forward=3e-5,
        E_a=45000.0,
        stoichiometry=1.0,
        ph_optimum=10.5,
        ph_min=9.0, ph_max=12.0,
        temperature_default=298.15,
        time_default=14400.0,
        spacer_length_angstrom=10.0,
        distal_group_yield=0.85,
        profile_role="spacer_intermediate",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="NHS-Sepharose HP standard; AHA provides -COOH distal",
        notes="Creates carboxyl terminus for EDC/NHS activation pathway",
    ),
    "edc_nhs_activation": ReagentProfile(
        name="EDC/NHS carboxyl activation",
        cas="25952-53-8",
        reagent_identity="EDC (1-Ethyl-3-(3-dimethylaminopropyl)carbodiimide) + NHS",
        installed_ligand="NHS ester (amine-reactive)",
        functional_mode="activator",
        reaction_type="activation",
        chemistry_class="edc_nhs",
        target_acs=ACSSiteType.CARBOXYL_DISTAL,
        product_acs=ACSSiteType.NHS_ESTER,
        k_forward=0.1,
        E_a=30000.0,
        stoichiometry=1.0,
        hydrolysis_rate=1e-4,
        ph_optimum=5.5,
        ph_min=4.5, ph_max=6.5,
        temperature_default=298.15,
        time_default=900.0,
        buffer_incompatibilities="Tris,glycine,primary_amine_buffers",
        profile_role="activated",
        m3_support_level="not_mapped",
        confidence_tier="ranking_only",
        calibration_source="Pseudo-single-step EDC/NHS; ranking_only (F5)",
        notes="Two-step mechanism simplified to single activation. NHS ester half-life ~2h at pH 7.5",
    ),

    # ═══════════════════════════════════════════════════════════════════
    # v5.9.4 — Washing (1 profile)
    # ═══════════════════════════════════════════════════════════════════

    "wash_buffer": ReagentProfile(
        name="Wash buffer (advisory residual removal)",
        cas="N/A",
        reagent_identity="Phosphate buffer pH 7.4",
        installed_ligand="N/A",
        functional_mode="washing",
        reaction_type="washing",
        chemistry_class="diffusion_out",
        target_acs=ACSSiteType.EPOXIDE, product_acs=None,
        k_forward=0.0, E_a=0.0, stoichiometry=1.0,
        ph_optimum=7.4,
        temperature_default=298.15,
        time_default=3600.0,
        regulatory_limit_ppm=1.0,
        profile_role="native",
        m3_support_level="not_mapped",
        confidence_tier="semi_quantitative",
        calibration_source="Advisory diffusion-out screening model",
        notes="Advisory only; does not claim GMP pass/fail without validated residual assays",
    ),
}
