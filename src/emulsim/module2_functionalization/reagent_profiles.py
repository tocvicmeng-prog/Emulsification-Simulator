"""Reagent profiles for Module 2 functionalization chemistry.

Phase B: Minimal Module 2 — 2 Workflows.
Architecture: module2_module3_final_implementation_plan.md, Phase B.

Only 4 reagent profiles needed for Phase B:
  - genipin_secondary: Secondary amine crosslinking (genipin)
  - glutaraldehyde_secondary: Secondary amine crosslinking (glutaraldehyde)
  - ech_activation: Hydroxyl activation with epichlorohydrin
  - dvs_activation: Hydroxyl activation with divinyl sulfone

Follows the CrosslinkerProfile pattern from reagent_library.py.
Literature values sourced from the Scientific Advisor report.
"""

from __future__ import annotations

from dataclasses import dataclass
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
        notes=(
            "Alkaline activation of agarose hydroxyl groups. "
            "Introduces vinyl sulfone for nucleophilic coupling. "
            "Porath & Fornstedt (1970) J. Chromatogr. 51:479. "
            "More hydrolytically stable than ECH epoxides."
        ),
    ),
}
