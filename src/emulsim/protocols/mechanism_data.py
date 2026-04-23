"""Mechanism data for reagent detail pages.

Each MechanismDescriptor provides a step-by-step mechanistic account of the
chemistry performed by a given reagent, keyed by the same string identifiers
used in reagent_library.CROSSLINKERS and
module2_functionalization.reagent_profiles.REAGENT_PROFILES.

Do NOT import streamlit here — this is pure data.

References
----------
Genipin:
    Butler et al. (2003) J. Polym. Sci. A 41:3941-3953
    Mi et al. (2005) Biomaterials 26:5983-5990
    Muzzarelli (2009) Carbohydr. Polym. 77:1-9
Glutaraldehyde / NaBH4:
    Migneault et al. (2004) BioTechniques 37:790-802
    Jentoft & Dearborn (1979) J. Biol. Chem. 254:4359-4365
EDC/NHS:
    Hermanson (2013) Bioconjugate Techniques, 3rd ed., Academic Press
    Grabarek & Gergely (1990) Anal. Biochem. 185:131-135
PEGDA / UV:
    Lin & Anseth (2009) Pharm. Res. 26:631-643
    Nguyen & West (2002) Biomaterials 23:4307-4314
TPP:
    Calvo et al. (1997) J. Appl. Polym. Sci. 63:125-132
Epichlorohydrin (ECH):
    Sundberg & Porath (1974) J. Chromatogr. 90:87-98
    Zhao et al. (2020) Carbohydr. Polym. 227:115352
Divinyl Sulfone (DVS):
    Porath et al. (1975) J. Chromatogr. 103:49-62
    Mateo et al. (2006) Enzyme Microb. Technol. 39:274-280
Citric Acid:
    Demitri et al. (2008) J. Appl. Polym. Sci. 110:2453-2460
    Reddy & Yang (2010) Food Chem. 118:702-711
BDGE:
    Gustavsson & Larsson (1996) J. Chromatogr. A 734:231-240
IDA/NTA chelation:
    Porath et al. (1975) Nature 258:598-599
    Sulkowski (1985) Trends Biotechnol. 3:1-7
Sultone ring-opening:
    Helferich & Schafer (1929) Liebigs Ann. 450:219
Nickel charging:
    Martell & Smith (1974-1989) Critical Stability Constants, Plenum Press
"""

from __future__ import annotations

from dataclasses import dataclass, field


# ═══════════════════════════════════════════════════════════════════════════
#  DATACLASSES
# ═══════════════════════════════════════════════════════════════════════════


@dataclass
class ReactionStep:
    """A single elementary or formal step in a reaction mechanism.

    Attributes
    ----------
    step_number:
        Sequential index starting at 1.
    description:
        Human-readable description, e.g. "Nucleophilic attack of -NH2 on
        aldehyde C=O".
    bond_formed:
        Bond(s) created in this step, e.g. "C-N (Schiff base, imine)".
    bond_broken:
        Bond(s) cleaved in this step, e.g. "C=O (aldehyde)".
    equation_latex:
        LaTeX-formatted equation for this step only (raw string).
    conditions:
        Reaction conditions, e.g. "pH 7-8, 25 C".
    notes:
        Optional additional commentary (safety, side reactions, etc.).
    """

    step_number: int
    description: str
    bond_formed: str
    bond_broken: str
    equation_latex: str
    conditions: str
    notes: str = ""


@dataclass
class MechanismDescriptor:
    """Full mechanistic description for one reagent.

    Attributes
    ----------
    reagent_key:
        Matches the key in CROSSLINKERS or REAGENT_PROFILES.
    display_name:
        Human-readable name shown in the UI.
    overall_equation_latex:
        Full balanced equation (raw LaTeX string).  Empty string means
        "Mechanism details not yet authored".
    mechanism_type:
        Category string — e.g. "nucleophilic_addition", "radical_chain",
        "ionic_interaction", "schiff_base", "carbodiimide_activation",
        "epoxide_ring_opening", "michael_addition", "ester_bond_formation",
        "metal_chelation", "reductive_amination", "sultone_ring_opening",
        "unknown".
    steps:
        Ordered list of ReactionStep objects.
    byproducts:
        Small-molecule byproducts released, e.g. ["H2O", "urea"].
    reversibility:
        One of "irreversible", "reversible_acid", "reversible_thermal",
        "reversible_ionic", "equilibrium".
    critical_notes:
        Safety, compatibility, or timing constraints shown prominently in
        the UI.
    svg_filename:
        Optional SVG file name in static/ folder, e.g.
        "genipin_mechanism.svg".  Empty string = no figure.
    buffer_system:
        Default recommended buffer, e.g. "PBS pH 7.4".
    quench_recommendation:
        How to quench after the reaction is complete.
    wash_protocol:
        Recommended post-reaction wash steps.
    storage_recommendation:
        Reagent and/or bead storage conditions.
    """

    reagent_key: str
    display_name: str
    overall_equation_latex: str
    mechanism_type: str
    steps: list[ReactionStep] = field(default_factory=list)
    byproducts: list[str] = field(default_factory=list)
    reversibility: str = "irreversible"
    critical_notes: list[str] = field(default_factory=list)
    svg_filename: str = ""
    buffer_system: str = "PBS pH 7.4"
    quench_recommendation: str = ""
    wash_protocol: str = ""
    storage_recommendation: str = ""


# ═══════════════════════════════════════════════════════════════════════════
#  MECHANISM_REGISTRY — manually authored entries
# ═══════════════════════════════════════════════════════════════════════════

MECHANISM_REGISTRY: dict[str, MechanismDescriptor] = {}

# ───────────────────────────────────────────────────────────────────────────
# 1. Genipin  (reagent_library.CROSSLINKERS["genipin"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["genipin"] = MechanismDescriptor(
    reagent_key="genipin",
    display_name="Genipin",
    overall_equation_latex=(
        r"2\,\text{Chitosan-NH}_2 + \text{Genipin} \rightarrow "
        r"\text{Chitosan-N-CH(ring)-N-Chitosan} + \text{H}_2\text{O}"
    ),
    mechanism_type="nucleophilic_addition",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Nucleophilic attack of chitosan primary amine (-NH2) on the "
                "C-3 position of the genipin iridoid ring, displacing the "
                "ester methoxy group and opening the dihydropyran ring to "
                "form a secondary amine intermediate."
            ),
            bond_formed="C3-N (secondary amine)",
            bond_broken="C3-OCH3 (methoxy ester)",
            equation_latex=(
                r"\text{Chitosan-NH}_2 + \text{Genipin(C3-OCH}_3\text{)} "
                r"\rightarrow \text{Chitosan-NH-Genipin(open)} + \text{CH}_3\text{OH}"
            ),
            conditions="pH 6.5-8.0, 25-37 C",
            notes=(
                "Rate-limiting step. kforward ~ 5e-6 m3/(mol*s) at 37 C "
                "(Butler et al. 2003). Slow: t1/2 ~ 4-24 h."
            ),
        ),
        ReactionStep(
            step_number=2,
            description=(
                "Intramolecular cyclisation of the ring-opened intermediate: "
                "the aldehyde terminus of the opened genipin attacks the "
                "alpha-carbon of the secondary amine to form a stable "
                "dihydropyridine (heterocyclic) crosslink on a single chain "
                "(monomer product)."
            ),
            bond_formed="C-N (dihydropyridine ring closure)",
            bond_broken="none (intramolecular ring formation)",
            equation_latex=(
                r"\text{Chitosan-NH-Genipin(open)} \rightarrow "
                r"\text{Chitosan-dihydropyridine} + \text{H}_2\text{O}"
            ),
            conditions="pH 6.5-8.0, spontaneous after step 1",
            notes=(
                "Produces intra-chain crosslinks. f_bridge is split between "
                "this pathway (~60%) and the intermolecular bridge (step 3)."
            ),
        ),
        ReactionStep(
            step_number=3,
            description=(
                "Intermolecular bridge formation: a second chitosan -NH2 "
                "reacts with the remaining reactive site on the genipin "
                "intermediate to form a covalent amine-bridge between two "
                "polymer chains. Subsequent auto-oxidation produces the "
                "characteristic blue pigment (oligomeric chromophore)."
            ),
            bond_formed="C-N (inter-chain amine bridge)",
            bond_broken="C=C (partial oxidation during pigment formation)",
            equation_latex=(
                r"\text{Chitosan-NH-Genipin} + \text{Chitosan'-NH}_2 "
                r"\rightarrow \text{Chitosan-N-Genipin-N-Chitosan'} "
                r"\xrightarrow{\text{O}_2} \text{blue pigment}"
            ),
            conditions="pH 6.5-8.0, 37 C, air exposure",
            notes=(
                "This step creates the elastically active crosslink modeled "
                "in Level 3 (f_bridge = 0.40). Blue coloration from pigment "
                "formation may interfere with absorbance monitoring at 280 nm "
                "if genipin is used in the functionalization workflow."
            ),
        ),
    ],
    byproducts=["methanol", "H2O", "blue_oligomeric_pigment"],
    reversibility="irreversible",
    critical_notes=[
        "Blue coloration develops over 24-72 h — may interfere with 280 nm "
        "UV absorbance monitoring of protein coupling yield.",
        "Low cytotoxicity (~10,000x less toxic than glutaraldehyde). "
        "Suitable for biocompatible applications.",
        "Slow kinetics require 24-48 h cure at 37 C for full crosslinking. "
        "Do not attempt to accelerate with heat above 50 C (agarose melts).",
        "Genipin is amine-specific and does NOT react with agarose hydroxyls "
        "under mild physiological conditions.",
    ],
    svg_filename="genipin_mechanism.svg",
    buffer_system="PBS pH 7.4",
    quench_recommendation="Dilute with excess PBS pH 7.4; 3x wash removes unreacted genipin.",
    wash_protocol="5 column volumes PBS pH 7.4, then 5 CV 0.5 M NaCl/PBS, then 5 CV PBS.",
    storage_recommendation="Store genipin powder at -20 C, desiccated. Solutions: use within 24 h.",
)

# ───────────────────────────────────────────────────────────────────────────
# 2. Glutaraldehyde  (reagent_library.CROSSLINKERS["glutaraldehyde"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["glutaraldehyde"] = MechanismDescriptor(
    reagent_key="glutaraldehyde",
    display_name="Glutaraldehyde",
    overall_equation_latex=(
        r"\text{Chitosan-NH}_2 + \text{OHC-(CH}_2\text{)}_3\text{-CHO} + "
        r"\text{H}_2\text{N-Chitosan'} \rightarrow "
        r"\text{Chitosan-N=CH-(CH}_2\text{)}_3\text{-CH=N-Chitosan'} + "
        r"2\,\text{H}_2\text{O}"
    ),
    mechanism_type="schiff_base",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Nucleophilic addition of the first chitosan primary amine "
                "(-NH2) to one aldehyde carbonyl of bifunctional "
                "glutaraldehyde (OHC-(CH2)3-CHO), forming a carbinolamine "
                "hemiaminal intermediate which rapidly dehydrates to yield "
                "a Schiff base (imine, C=N) at one end."
            ),
            bond_formed="C=N (Schiff base / imine) at aldehyde-1",
            bond_broken="C=O (aldehyde-1)",
            equation_latex=(
                r"\text{Chitosan-NH}_2 + \text{OHC-R} \rightarrow "
                r"\text{Chitosan-N=CH-R} + \text{H}_2\text{O}"
            ),
            conditions="pH 6.5-8.0, 25 C, 5-30 min",
            notes=(
                "Fast reaction (k ~ 1 L/(mol*s) at 25 C). The resulting "
                "mono-Schiff base tether is reactive at the free aldehyde end."
            ),
        ),
        ReactionStep(
            step_number=2,
            description=(
                "The free aldehyde at the other end of the glutaraldehyde "
                "bridge reacts with a second chitosan primary amine on an "
                "adjacent or distal chain, forming the second Schiff base "
                "and completing the bifunctional inter-chain crosslink. "
                "The product imine linkages (C=N) are hydrolytically "
                "reversible unless subsequently reduced with NaBH4."
            ),
            bond_formed="C=N (Schiff base / imine) at aldehyde-2",
            bond_broken="C=O (aldehyde-2)",
            equation_latex=(
                r"\text{Chitosan-N=CH-(CH}_2\text{)}_3\text{-CHO} + "
                r"\text{H}_2\text{N-Chitosan'} \rightarrow "
                r"\text{Chitosan-N=CH-(CH}_2\text{)}_3\text{-CH=N-Chitosan'} + "
                r"\text{H}_2\text{O}"
            ),
            conditions="pH 6.5-8.0, 25 C, 30-60 min",
            notes=(
                "Schiff base C=N is hydrolytically labile at pH < 5. "
                "For permanent crosslinks, reduce with NaBH4 immediately "
                "after washing (see nabh4_quench profile)."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="reversible_acid",
    critical_notes=[
        "HIGH TOXICITY: LD50 oral ~100 mg/kg (rat). Handle in fume hood with "
        "nitrile gloves. GHS Category 3 acute toxicity.",
        "Schiff base linkages are hydrolytically unstable at pH < 5. Reduce "
        "with 10 mM NaBH4 in PBS pH 7-8 for 30 min to convert C=N to stable "
        "C-NH (secondary amine).",
        "Residual glutaraldehyde leaches from beads and can modify immobilized "
        "proteins — thorough washing (>10 CV) is mandatory before use in "
        "chromatography.",
        "Produces yellow coloration at high crosslink densities.",
    ],
    buffer_system="PBS pH 7.0",
    quench_recommendation=(
        "Add 0.1 M glycine or ethanolamine in PBS pH 7.4 for 30 min to cap "
        "remaining free aldehyde groups, then wash 10 CV PBS."
    ),
    wash_protocol=(
        "5 CV 0.1 M glycine in PBS pH 7.4 (quench), then 10 CV PBS pH 7.4, "
        "then 5 CV 0.5 M NaCl/PBS, then 5 CV PBS."
    ),
    storage_recommendation=(
        "Store stock solution (25% in water) at 4 C, amber glass, up to 12 months. "
        "Working solutions: prepare fresh daily."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 3. EDC/NHS  (reagent_library.CROSSLINKERS["edc_nhs"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["edc_nhs"] = MechanismDescriptor(
    reagent_key="edc_nhs",
    display_name="EDC/NHS (carbodiimide activation)",
    overall_equation_latex=(
        r"\text{R-COOH} + \text{EDC} + \text{NHS} \rightarrow "
        r"\text{R-CO-NHS} \xrightarrow{\text{R'-NH}_2} "
        r"\text{R-CO-NH-R'} + \text{NHS} + \text{urea}"
    ),
    mechanism_type="carbodiimide_activation",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "EDC (1-ethyl-3-(3-dimethylaminopropyl)carbodiimide) activates "
                "a carboxyl group (-COOH) by forming a highly reactive "
                "O-acylisourea intermediate. This step is rapid but the "
                "O-acylisourea hydrolyzes back to -COOH with t1/2 ~ 10-20 min "
                "at pH 5, so subsequent steps must proceed quickly."
            ),
            bond_formed="C-O (O-acylisourea mixed anhydride)",
            bond_broken="O-H (carboxylate proton), rearranges C=N of carbodiimide",
            equation_latex=(
                r"\text{R-COOH} + \text{R''-N=C=N-R'''} \rightarrow "
                r"\text{R-CO-O-C(=NR'')-NHR'''}"
            ),
            conditions="pH 4.5-6.0, 4-25 C, 15-30 min",
            notes=(
                "O-acylisourea is unstable — hydrolyzes rapidly above pH 6 to "
                "regenerate -COOH and urea byproduct. Keep at pH 5 and 4 C "
                "to prolong intermediate lifetime."
            ),
        ),
        ReactionStep(
            step_number=2,
            description=(
                "NHS (N-hydroxysuccinimide) attacks the O-acylisourea to "
                "displace the urea and form a semi-stable NHS ester. The NHS "
                "ester is hydrolytically more stable (t1/2 ~ 2-4 h at pH 7.4, "
                "25 C) and can be stored or used in a second buffer step at "
                "physiological pH."
            ),
            bond_formed="CO-O-NHS ester bond",
            bond_broken="C-O (O-acylisourea) → releases urea",
            equation_latex=(
                r"\text{R-CO-O-acylisourea} + \text{NHS-OH} \rightarrow "
                r"\text{R-CO-O-NHS} + \text{urea}"
            ),
            conditions="pH 5.0-6.0, immediate addition of NHS",
            notes=(
                "NHS ester formation is the key step that allows the activated "
                "carboxyl to survive pH adjustment to 7.2-7.4 for amine coupling."
            ),
        ),
        ReactionStep(
            step_number=3,
            description=(
                "The NHS ester reacts with a primary amine (-NH2) on chitosan "
                "or the ligand of interest via aminolysis, releasing NHS and "
                "forming a stable amide bond (peptide-like C-CO-NH-). "
                "This is the final, permanent covalent linkage."
            ),
            bond_formed="C-N amide bond (zero-length crosslink)",
            bond_broken="CO-O-NHS ester",
            equation_latex=(
                r"\text{R-CO-O-NHS} + \text{R'-NH}_2 \rightarrow "
                r"\text{R-CO-NH-R'} + \text{NHS-OH}"
            ),
            conditions="pH 7.0-8.0, 4-25 C, 1-4 h",
            notes=(
                "Amide bond is permanent and hydrolytically stable across "
                "pH 2-12. NHS is released as a non-toxic byproduct. "
                "Zero-length crosslink — no spacer atom remains."
            ),
        ),
    ],
    byproducts=["urea", "NHS"],
    reversibility="irreversible",
    critical_notes=[
        "Requires carboxyl groups (-COOH) on the substrate. Pure chitosan "
        "has very few COOH groups — use carboxymethyl chitosan or add a "
        "dicarboxylic acid spacer (e.g. adipic acid) first.",
        "EDC hydrolyzes rapidly in water: t1/2 ~ 10-20 min at pH 5, 25 C. "
        "Prepare fresh solution and add immediately to reaction. Do NOT "
        "pre-mix EDC and NHS without substrate.",
        "Use pH 4.5-6.0 for activation step, then adjust to pH 7.0-8.0 for "
        "amine coupling. Two-step protocol avoids competing hydrolysis.",
        "Avoid Tris, glycine, and other amine-containing buffers — they "
        "compete with the target amine for the NHS ester.",
        "EDC is moisture-sensitive: store at -20 C desiccated, warm to RT "
        "before opening to prevent condensation.",
    ],
    buffer_system="MES pH 5.5 (activation), PBS pH 7.4 (coupling)",
    quench_recommendation=(
        "Add 20 mM ethanolamine pH 8.0 for 30 min to quench remaining "
        "NHS esters. EDC itself is quenched by water (urea byproduct)."
    ),
    wash_protocol=(
        "5 CV 20 mM ethanolamine in PBS pH 8.0 (quench NHS esters), "
        "then 10 CV PBS pH 7.4."
    ),
    storage_recommendation=(
        "Store EDC powder at -20 C in sealed desiccant container. "
        "NHS at -20 C. Working solutions: prepare fresh immediately before use."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 4. PEGDA + UV  (reagent_library.CROSSLINKERS["pegda_uv"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["pegda_uv"] = MechanismDescriptor(
    reagent_key="pegda_uv",
    display_name="PEGDA + UV photoinitiation",
    overall_equation_latex=(
        r"n\,\text{CH}_2\text{=CH-CO-O-PEG-O-CO-CH=CH}_2 "
        r"\xrightarrow{h\nu,\,\text{Irgacure}} "
        r"\text{[-CH}_2\text{-CH(CO-O-PEG-O-CO-CH}_2\text{-CH-)]}_n "
        r"\text{(crosslinked PEG network)}"
    ),
    mechanism_type="radical_chain",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "UV absorption by photoinitiator (Irgacure 2959 or LAP at "
                "365 nm / 405 nm) causes homolytic cleavage of the C-C bond "
                "alpha to the carbonyl, generating two carbon-centred radicals. "
                "Each radical is capable of initiating acrylate polymerization."
            ),
            bond_formed="radical species (C-radical)",
            bond_broken="C-C (alpha to carbonyl, Norrish type I cleavage)",
            equation_latex=(
                r"\text{PI} \xrightarrow{h\nu} 2\,\text{R}^\bullet"
            ),
            conditions="365 nm or 405 nm UV, 1-10 mW/cm2, photoinitiator 0.1-1% w/v",
            notes=(
                "Irgacure 2959 is water-soluble but weakly cytotoxic at >0.1% w/v. "
                "LAP (lithium acylphosphinate) is preferred for cell-laden gels. "
                "UV penetration depth in opaque beads is limited to ~50-100 um."
            ),
        ),
        ReactionStep(
            step_number=2,
            description=(
                "Radical addition to PEGDA acrylate double bond initiates "
                "chain-growth polymerization. Each propagation step adds one "
                "PEGDA monomer and extends the growing radical chain."
            ),
            bond_formed="C-C (polymer backbone)",
            bond_broken="C=C (acrylate vinyl)",
            equation_latex=(
                r"\text{R}^\bullet + \text{CH}_2\text{=CH-PEGDA} \rightarrow "
                r"\text{R-CH}_2\text{-}\overset{\bullet}{\text{CH}}\text{-PEGDA}"
            ),
            conditions="RT, aqueous, UV exposure",
            notes=(
                "Propagation rate kp ~ 10-100 L/(mol*s) for acrylates. "
                "Oxygen inhibition is significant in air-exposed samples — "
                "degas solutions or use inert atmosphere for high conversion."
            ),
        ),
        ReactionStep(
            step_number=3,
            description=(
                "Radical termination (combination or disproportionation) "
                "and chain-transfer events create a PEG network that "
                "interpenetrates the pre-existing agarose/chitosan matrix, "
                "forming an independent triple IPN. No covalent bonds are "
                "formed between PEG and chitosan/agarose in this step."
            ),
            bond_formed="C-C (chain termination by combination)",
            bond_broken="none (radical quenching)",
            equation_latex=(
                r"2\,\overset{\bullet}{\text{P}} \rightarrow \text{P-P} "
                r"\quad (\text{termination, independent PEG network})"
            ),
            conditions="RT, continuous UV until desired dose reached",
            notes=(
                "The resulting PEG network is physically entangled with but "
                "NOT covalently bonded to agarose or chitosan. "
                "PEG reduces non-specific protein binding — advantage for "
                "high-specificity chromatography applications."
            ),
        ),
    ],
    byproducts=["photoinitiator_fragments"],
    reversibility="irreversible",
    critical_notes=[
        "UV penetration in opaque agarose-chitosan microspheres is limited to "
        "~50-100 um from the surface. Beads > 100 um diameter will have an "
        "under-crosslinked core.",
        "Oxygen inhibits radical polymerization — degas PEGDA solution with N2 "
        "for 15 min before UV exposure.",
        "PEGDA does NOT react directly with chitosan -NH2 or agarose -OH. It "
        "forms a completely independent network. Ensure gel is fully formed "
        "before adding PEGDA.",
        "Irgacure 2959 photoinitiator fragments may remain in the bead — "
        "wash thoroughly with PBS (minimum 5 CV) before use.",
        "Photoinitiator is cytotoxic at concentrations > 0.5% w/v — use lowest "
        "effective concentration and maximize wash volume.",
    ],
    buffer_system="PBS pH 7.4 (deoxygenated)",
    quench_recommendation="Remove UV source; residual radicals quench spontaneously within seconds.",
    wash_protocol="10 CV PBS pH 7.4, then 5 CV 50% ethanol/PBS, then 10 CV PBS pH 7.4.",
    storage_recommendation=(
        "Store PEGDA at -20 C away from light. Photoinitiator: -20 C, amber vial. "
        "Crosslinked beads: 4 C in PBS + 0.02% sodium azide."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 5. TPP — sodium tripolyphosphate  (reagent_library.CROSSLINKERS["tpp"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["tpp"] = MechanismDescriptor(
    reagent_key="tpp",
    display_name="Sodium Tripolyphosphate (TPP) — ionic crosslinking",
    overall_equation_latex=(
        r"\text{Chitosan-NH}_3^+ + \tfrac{1}{5}\,\text{P}_3\text{O}_{10}^{5-} "
        r"\rightleftharpoons \text{Chitosan-NH}_3^+\cdots "
        r"\text{P}_3\text{O}_{10}^{5-}\cdots\text{NH}_3^+\text{-Chitosan}"
    ),
    mechanism_type="ionic_interaction",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Electrostatic complexation: the pentavalent tripolyphosphate "
                "anion (P3O10^5-) electrostatically bridges multiple protonated "
                "glucosamine ammonium groups (-NH3+, pKa ~6.5) on chitosan "
                "chains. At pH < 6.5, chitosan is sufficiently protonated for "
                "this interaction. Gelation is effectively instantaneous "
                "(seconds to minutes) upon mixing."
            ),
            bond_formed="electrostatic ion pair (-NH3+...OPO3-)",
            bond_broken="none (non-covalent)",
            equation_latex=(
                r"5\,\text{Chit-NH}_3^+ + \text{P}_3\text{O}_{10}^{5-} "
                r"\rightleftharpoons "
                r"[\text{Chit-NH}_3^+]_5[\text{P}_3\text{O}_{10}^{5-}]"
            ),
            conditions="pH < 6.5, RT, ionic strength < 0.3 M",
            notes=(
                "Instantaneous — t1/2 < 30 s. The crosslink density depends "
                "on the TPP:chitosan mass ratio and pH. Optimal ratio: "
                "typically 3:1 to 10:1 chitosan:TPP by mass."
            ),
        ),
    ],
    byproducts=[],
    reversibility="reversible_ionic",
    critical_notes=[
        "Ionic crosslinks are REVERSIBLE and pH-sensitive. Above pH 6.5 "
        "chitosan deprotonates (-NH2 form) and ionic bridges dissociate. "
        "Under typical chromatography buffer conditions (PBS pH 7.4, 0.15 M NaCl) "
        "TPP crosslinks WILL partially or fully dissociate.",
        "Mechanical strength is low (G ~ 0.1-1 kPa) — NOT suitable as sole "
        "crosslinker for pressure-driven column chromatography.",
        "High ionic strength competes with and weakens ionic crosslinks. "
        "Not suitable for salt gradient separations.",
        "TPP is food-grade and non-toxic (E451 additive). Safe to handle.",
    ],
    buffer_system="Acetic acid / sodium acetate pH 5.0",
    quench_recommendation="Dilute with 0.1 M NaCl to screen electrostatics; not applicable.",
    wash_protocol="3 CV deionized water, then 3 CV acetate buffer pH 5.0.",
    storage_recommendation="Store TPP powder at RT, desiccated. Aqueous solutions: stable weeks at 4 C.",
)

# ───────────────────────────────────────────────────────────────────────────
# 6. Epichlorohydrin  (reagent_library.CROSSLINKERS["epichlorohydrin"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["epichlorohydrin"] = MechanismDescriptor(
    reagent_key="epichlorohydrin",
    display_name="Epichlorohydrin (ECH)",
    overall_equation_latex=(
        r"\text{Polymer-OH} + \text{ClCH}_2\text{-CH(O)-CH}_2 "
        r"\xrightarrow{\text{NaOH}} "
        r"\text{Polymer-O-CH}_2\text{-CHOH-CH}_2\text{-O-Polymer'} + "
        r"\text{HCl(neutralized)}"
    ),
    mechanism_type="epoxide_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Under strongly alkaline conditions (0.5-2 M NaOH), "
                "a polysaccharide hydroxyl group (-OH) is deprotonated to an "
                "alkoxide (-O-), which acts as a nucleophile and opens the "
                "epoxide ring of epichlorohydrin (ECH) via SN2 attack at the "
                "less-hindered C-1 position, forming a chlorohydrin ether "
                "tether covalently bound to the polymer at one end."
            ),
            bond_formed="C-O ether (polymer-O-CH2-CHOH-CH2Cl)",
            bond_broken="C-O (epoxide ring opening at C1)",
            equation_latex=(
                r"\text{Polymer-O}^- + \text{ClCH}_2\text{CH(O)CH}_2 "
                r"\rightarrow \text{Polymer-O-CH}_2\text{-CHOH-CH}_2\text{Cl}"
            ),
            conditions="0.5-2 M NaOH, 50 C, 2-6 h",
            notes=(
                "ECH can also react with -NH2 groups at lower pH (~8-10). "
                "At pH > 12, the -OH pathway dominates (Sundberg & Porath 1974)."
            ),
        ),
        ReactionStep(
            step_number=2,
            description=(
                "Base-catalysed intramolecular HCl elimination from the "
                "chlorohydrin intermediate regenerates an epoxide group "
                "(now tethered to the polymer via an ether bridge), OR "
                "direct SN2 displacement of Cl- by a second polymer -O- "
                "forms the bis-ether crosslink in one step. "
                "Either pathway yields the stable glyceryl ether crosslink "
                "(-O-CH2-CHOH-CH2-O-) between two polymer chains."
            ),
            bond_formed="C-O ether (second polymer chain crosslink)",
            bond_broken="C-Cl (chloride departure, SN2) or O-H + C-Cl (elimination → epoxide)",
            equation_latex=(
                r"\text{Polymer-O-CH}_2\text{-CHOH-CH}_2\text{Cl} + "
                r"\text{Polymer'-O}^- \rightarrow "
                r"\text{Polymer-O-CH}_2\text{-CHOH-CH}_2\text{-O-Polymer'} + "
                r"\text{Cl}^-"
            ),
            conditions="0.5-2 M NaOH, 50 C, continuation of step 1 reaction",
            notes=(
                "The glyceryl ether bridge is extremely hydrolytically stable "
                "— resistant to acid, base, urea, and guanidine HCl. "
                "This is the industrial process for producing Sepharose CL."
            ),
        ),
    ],
    byproducts=["NaCl", "H2O"],
    reversibility="irreversible",
    critical_notes=[
        "TOXIC AND CARCINOGENIC — IARC Group 1 confirmed human carcinogen. "
        "Handle exclusively in a certified fume hood with nitrile or neoprene "
        "gloves. Avoid all skin/eye contact and inhalation.",
        "Strongly alkaline conditions (pH > 12) required for -OH crosslinking "
        "WILL disrupt agarose helix structure if gel is not fully set and "
        "stabilized first. Apply ECH only AFTER complete agarose gelation.",
        "Competing hydrolysis of ECH is significant at pH > 11. Excess ECH "
        "must be thoroughly washed out — validated wash protocols required.",
        "ECH also reacts with -NH2 at lower pH (~8-10), which may compete "
        "with the desired -OH crosslinking in chitosan-containing matrices.",
        "Validated wash protocol (e.g. Pharmacia Sepharose CL protocol) is "
        "required before using crosslinked beads in any biological assay.",
    ],
    buffer_system="0.5-2 M NaOH (reaction), PBS pH 7.4 (post-wash)",
    quench_recommendation=(
        "Dilute reaction with 10 volumes of cold water; neutralize to pH 7.4 "
        "with HCl; wash extensively."
    ),
    wash_protocol=(
        "Neutralize with 0.1 M HCl to pH 7.0, then 10 CV water, "
        "then 5 CV 0.5 M NaCl, then 10 CV PBS pH 7.4. "
        "Monitor eluate ECH by GC-headspace < 1 ppm before use."
    ),
    storage_recommendation=(
        "Store ECH in sealed bottles at 4-8 C under nitrogen. "
        "Crosslinked beads: 4 C in PBS + 20% ethanol as preservative."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 7. DVS — divinyl sulfone  (reagent_library.CROSSLINKERS["dvs"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["dvs"] = MechanismDescriptor(
    reagent_key="dvs",
    display_name="Divinyl Sulfone (DVS)",
    overall_equation_latex=(
        r"\text{Polymer-OH} + \text{CH}_2\text{=CH-SO}_2\text{-CH=CH}_2 + "
        r"\text{HO-Polymer'} \xrightarrow{\text{NaOH, pH 11}} "
        r"\text{Polymer-O-CH}_2\text{CH}_2\text{-SO}_2\text{-CH}_2\text{CH}_2"
        r"\text{-O-Polymer'}"
    ),
    mechanism_type="michael_addition",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Under alkaline conditions (0.5 M NaOH + 0.5 M Na2CO3, pH 11-12), "
                "the polysaccharide alkoxide (-O-) undergoes oxa-Michael addition "
                "to one vinyl group of divinyl sulfone (DVS). The sulfone group "
                "strongly activates both vinyl groups toward nucleophilic addition. "
                "This step couples DVS to the first polymer chain via an ether-sulfone "
                "tether, leaving one vinyl group free for the second coupling."
            ),
            bond_formed="C-O ether + vinyl sulfone addition product (C-C bond at beta-C)",
            bond_broken="C=C (one vinyl, Michael acceptor)",
            equation_latex=(
                r"\text{Polymer-O}^- + \text{CH}_2\text{=CH-SO}_2\text{-CH=CH}_2 "
                r"\rightarrow \text{Polymer-O-CH}_2\text{CH}_2\text{-SO}_2"
                r"\text{-CH=CH}_2"
            ),
            conditions="0.5 M NaOH / 0.5 M Na2CO3, pH 11-12, 25 C, 30-60 min",
            notes=(
                "DVS reacts faster than ECH at equivalent pH due to the "
                "electron-withdrawing sulfone activating both vinyl groups. "
                "Both -OH and -NH2 nucleophiles react; -OH pathway dominates "
                "at pH 11-12."
            ),
        ),
        ReactionStep(
            step_number=2,
            description=(
                "The remaining free vinyl sulfone on the DVS mono-adduct "
                "undergoes a second oxa-Michael (or aza-Michael) addition with "
                "a nucleophile on a second polymer chain (-OH or -NH2), completing "
                "the bis-ether-sulfone crosslink. The resulting bridge "
                "(-O-CH2CH2-SO2-CH2CH2-O-) is chemically inert to acids, bases, "
                "urea, and chaotropes."
            ),
            bond_formed="C-O ether (second chain) completing the crosslink bridge",
            bond_broken="C=C (second vinyl, Michael acceptor)",
            equation_latex=(
                r"\text{Polymer-O-CH}_2\text{CH}_2\text{-SO}_2\text{-CH=CH}_2 + "
                r"\text{HO-Polymer'} \xrightarrow{\text{base}} "
                r"\text{Polymer-O-CH}_2\text{CH}_2\text{-SO}_2\text{-CH}_2"
                r"\text{CH}_2\text{-O-Polymer'}"
            ),
            conditions="0.5 M NaOH / 0.5 M Na2CO3, pH 11-12, 25 C, 1-2 h total",
            notes=(
                "Crosslinked DVS agarose (Sepharose HP) shows superior rigidity and "
                "pressure resistance vs ECH-crosslinked material. "
                "f_bridge ~ 0.70 (high bridge efficiency, Porath et al. 1975)."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="irreversible",
    critical_notes=[
        "TOXIC AND LACHRYMATORY — DVS vapour is a potent eye/respiratory "
        "irritant. Handle exclusively in a certified fume hood with acid-gas "
        "cartridge respirator and chemical splash goggles.",
        "Strongly alkaline conditions (pH 11-12) are required — confirm agarose "
        "gel is fully crosslinked and stable before DVS treatment.",
        "Unreacted vinyl sulfone groups hydrolyze to -CH2CH2-SO2-OH (harmless "
        "sulfonate) within hours at physiological pH — no quench step strictly "
        "required, but washing is mandatory.",
        "DVS reacts with -NH2 groups (aza-Michael) as well as -OH. In "
        "chitosan-agarose blends this provides additional chitosan-DVS-agarose "
        "inter-network crosslinks.",
        "Industry gold standard for Sepharose CL and Sepharose HP production "
        "(Cytiva/GE Healthcare). Well-validated wash protocols available.",
    ],
    buffer_system="0.5 M NaOH + 0.5 M Na2CO3 (reaction), PBS pH 7.4 (post-wash)",
    quench_recommendation=(
        "Residual vinyl sulfone groups quench spontaneously in aqueous buffer "
        "(hydrolysis t1/2 ~ 2-4 h at pH 7.4). Wash 10 CV PBS immediately "
        "after reaction to remove excess DVS."
    ),
    wash_protocol=(
        "Neutralize pH to 7.0 with HCl, then 10 CV water, "
        "then 5 CV 0.5 M NaCl, then 10 CV PBS pH 7.4. "
        "Monitor DVS residual by Ellman's assay if regulatory compliance required."
    ),
    storage_recommendation=(
        "Store DVS at 4 C under nitrogen in amber glass. "
        "Crosslinked beads: 4 C in PBS + 20% ethanol."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 8. Citric acid  (reagent_library.CROSSLINKERS["citric_acid"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["citric_acid"] = MechanismDescriptor(
    reagent_key="citric_acid",
    display_name="Citric Acid (heat-activated ester crosslinking)",
    overall_equation_latex=(
        r"\text{Polymer-OH} + \text{HOOC-CH}_2\text{-C(OH)(COOH)-CH}_2\text{-COOH} "
        r"\xrightarrow{80\text{-}120\,^\circ\text{C}} "
        r"\text{Polymer-OOC-CH}_2\text{-C(OH)(COOH)-CH}_2\text{-COO-Polymer'} + "
        r"2\,\text{H}_2\text{O}"
    ),
    mechanism_type="ester_bond_formation",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Heat activation (80-120 C) drives Fischer esterification: "
                "a carboxyl group (-COOH) of citric acid condenses with a "
                "polysaccharide hydroxyl group (-OH) to form an ester bond "
                "(R-COO-R'), releasing water. Citric acid is trifunctional "
                "(3 COOH groups), so one or more can participate in bridge "
                "formation while remaining COOH groups stay pendant."
            ),
            bond_formed="C-O ester bond",
            bond_broken="O-H (hydroxyl proton) and C=O → C-O-C ester",
            equation_latex=(
                r"\text{Polymer-OH} + \text{HOOC-Citrate} "
                r"\xrightarrow{\Delta,\,-H_2\text{O}} "
                r"\text{Polymer-OOC-Citrate} + \text{H}_2\text{O}"
            ),
            conditions="80-120 C, dry or minimal aqueous conditions, 2-4 h",
            notes=(
                "Esterification is an equilibrium — water must be removed "
                "(heat, vacuum, or dry atmosphere) to drive conversion. "
                "At temperatures > 120 C, dehydration of citric acid itself "
                "to aconitic acid anhydride can occur."
            ),
        ),
        ReactionStep(
            step_number=2,
            description=(
                "A second esterification (or amide bond formation with -NH2 "
                "of chitosan) at a different COOH group of the same citrate "
                "molecule creates the inter-chain crosslink, bridging two "
                "polymer chains through the central citrate unit. "
                "Trifunctional citrate can in principle form up to three "
                "ester/amide bonds, but steric constraints typically limit "
                "bridge formation to 1-2 bonds per citrate."
            ),
            bond_formed="C-O ester (or C-N amide with chitosan -NH2) — second chain",
            bond_broken="O-H or N-H of second nucleophile",
            equation_latex=(
                r"\text{Polymer-OOC-Citrate-COOH} + \text{HO-Polymer'} "
                r"\xrightarrow{\Delta} "
                r"\text{Polymer-OOC-Citrate-COO-Polymer'} + \text{H}_2\text{O}"
            ),
            conditions="80-120 C, continuation of step 1",
            notes=(
                "f_bridge ~ 0.25 because many citrate COOH groups remain "
                "pendant (unreacted). Pendant carboxyls add negative charge "
                "density — may introduce unintended weak cation-exchange "
                "character in chromatography beads."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="reversible_acid",
    critical_notes=[
        "Heat cure at 80 C+ is required — DO NOT attempt at room temperature. "
        "CAUTION: If agarose gel is present, temperatures > 40 C will melt the "
        "agarose helix network. Apply citric acid crosslinking ONLY to pre-formed, "
        "cooled beads, or in a two-step process after agarose gelation.",
        "Ester bonds are hydrolytically unstable at pH > 10 or pH < 3 — limits "
        "use in strong regeneration protocols (NaOH cleaning-in-place).",
        "Citric acid is FDA GRAS and non-toxic — suitable for food, cosmetic, "
        "and biocompatible single-use applications.",
        "Residual pendant carboxyl groups from partially reacted citrate introduce "
        "weak cation-exchange character — characterize before chromatography use.",
    ],
    buffer_system="Aqueous (minimal water to drive equilibrium toward ester)",
    quench_recommendation=(
        "Cool to RT; wash with PBS pH 7.4 to remove unreacted citric acid."
    ),
    wash_protocol=(
        "5 CV water at 60 C (removes citric acid residuals), "
        "then 5 CV PBS pH 7.4. Confirm no citric acid in eluate by UV 210 nm."
    ),
    storage_recommendation=(
        "Citric acid: store at RT, desiccated, indefinitely. "
        "Crosslinked beads: 4 C in PBS; ester stability above pH 3, below pH 10."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 9. Genipin secondary (reagent_profiles.REAGENT_PROFILES["genipin_secondary"])
# ───────────────────────────────────────────────────────────────────────────
# Same mechanism as genipin — create a reference copy.
MECHANISM_REGISTRY["genipin_secondary"] = MechanismDescriptor(
    reagent_key="genipin_secondary",
    display_name="Genipin (secondary crosslinking)",
    overall_equation_latex=MECHANISM_REGISTRY["genipin"].overall_equation_latex,
    mechanism_type="nucleophilic_addition",
    steps=MECHANISM_REGISTRY["genipin"].steps,
    byproducts=MECHANISM_REGISTRY["genipin"].byproducts,
    reversibility="irreversible",
    critical_notes=[
        "Identical mechanism to primary genipin crosslinking (see 'genipin' entry).",
        "This profile is used for the Module 2 secondary functionalization step "
        "(post-activation crosslinking of residual amines on the bead surface).",
        "Blue coloration from pigment formation accumulates with each genipin "
        "treatment — plan optical monitoring accordingly.",
    ],
    svg_filename="genipin_mechanism.svg",
    buffer_system="PBS pH 7.4",
    quench_recommendation=MECHANISM_REGISTRY["genipin"].quench_recommendation,
    wash_protocol=MECHANISM_REGISTRY["genipin"].wash_protocol,
    storage_recommendation=MECHANISM_REGISTRY["genipin"].storage_recommendation,
)

# ───────────────────────────────────────────────────────────────────────────
# 10. Glutaraldehyde secondary (reagent_profiles["glutaraldehyde_secondary"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["glutaraldehyde_secondary"] = MechanismDescriptor(
    reagent_key="glutaraldehyde_secondary",
    display_name="Glutaraldehyde (secondary crosslinking)",
    overall_equation_latex=MECHANISM_REGISTRY["glutaraldehyde"].overall_equation_latex,
    mechanism_type="schiff_base",
    steps=MECHANISM_REGISTRY["glutaraldehyde"].steps,
    byproducts=MECHANISM_REGISTRY["glutaraldehyde"].byproducts,
    reversibility="reversible_acid",
    critical_notes=MECHANISM_REGISTRY["glutaraldehyde"].critical_notes,
    buffer_system="PBS pH 7.0",
    quench_recommendation=MECHANISM_REGISTRY["glutaraldehyde"].quench_recommendation,
    wash_protocol=MECHANISM_REGISTRY["glutaraldehyde"].wash_protocol,
    storage_recommendation=MECHANISM_REGISTRY["glutaraldehyde"].storage_recommendation,
)

# ───────────────────────────────────────────────────────────────────────────
# 11. ECH activation (reagent_profiles["ech_activation"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["ech_activation"] = MechanismDescriptor(
    reagent_key="ech_activation",
    display_name="Epichlorohydrin (OH activation — introduces epoxide)",
    overall_equation_latex=(
        r"\text{Agarose-OH} + \text{ClCH}_2\text{-CH(O)-CH}_2 "
        r"\xrightarrow{\text{NaOH, pH 12}} "
        r"\text{Agarose-O-CH}_2\text{-CH(O)-CH}_2 + \text{HCl(neutralized)}"
    ),
    mechanism_type="epoxide_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Alkaline deprotonation of agarose hydroxyl to alkoxide (-O-), "
                "followed by SN2 attack on the epoxide ring of ECH at the "
                "less-hindered methylene carbon (C-1), forming an ether-linked "
                "chlorohydrin intermediate on the agarose chain."
            ),
            bond_formed="C-O ether (agarose-O-CH2-CHOH-CH2Cl)",
            bond_broken="C-O (epoxide ring at C1)",
            equation_latex=(
                r"\text{Agarose-O}^- + \text{ClCH}_2\text{CH(O)CH}_2 "
                r"\rightarrow \text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{Cl}"
            ),
            conditions="0.5-2 M NaOH, 25-50 C, 2-4 h",
            notes=(
                "Competing hydrolysis of ECH epoxide by water/hydroxide "
                "reduces yield. Keep temperature at 25-37 C and not too "
                "alkaline (pH 11-12) to balance rate vs hydrolysis."
            ),
        ),
        ReactionStep(
            step_number=2,
            description=(
                "Base-catalysed elimination of HCl from the chlorohydrin "
                "intermediate regenerates an epoxide ring now tethered to "
                "agarose via a stable ether bond. This epoxide group is then "
                "reactive toward nucleophiles (amines, thiols) at pH 8-10 "
                "for ligand coupling in subsequent Module 2 steps."
            ),
            bond_formed="C-O (epoxide ring re-formed, now agarose-tethered)",
            bond_broken="C-Cl and O-H (intramolecular E2 elimination of HCl)",
            equation_latex=(
                r"\text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{Cl} "
                r"\xrightarrow{\text{NaOH}} "
                r"\text{Agarose-O-CH}_2\text{-CH(O)-CH}_2 + \text{NaCl}"
            ),
            conditions="0.5 M NaOH, RT, spontaneous during step 1",
            notes=(
                "The product is an epoxide-activated agarose bead. "
                "Epoxide density (umol/mL) is measured by titration with "
                "sodium thiosulfate or Ellman's reagent after thiol coupling. "
                "Hydrolysis of the tethered epoxide (t1/2 ~ 12-48 h at pH 7.4) "
                "sets the working window for subsequent ligand coupling."
            ),
        ),
    ],
    byproducts=["NaCl", "H2O"],
    reversibility="irreversible",
    critical_notes=[
        "TOXIC AND CARCINOGENIC — IARC Group 1. Handle only in certified fume "
        "hood with appropriate PPE. See 'epichlorohydrin' mechanism entry.",
        "Activated epoxide beads must be used promptly — epoxide hydrolyzes at "
        "pH 7.4 with t1/2 ~ 12-48 h. Store at pH 4-5 to slow hydrolysis, "
        "then adjust pH before ligand coupling.",
        "Strongly alkaline activation conditions (pH > 12) may partially dissolve "
        "agarose if applied before adequate gelation. Confirm gel integrity first.",
        "Epoxide density should be validated by colorimetric assay before "
        "proceeding to ligand coupling to ensure reproducible ligand loading.",
    ],
    buffer_system="0.5-2 M NaOH (activation), sodium acetate pH 4.0 (storage)",
    quench_recommendation=(
        "Quench remaining epoxides with 1 M ethanolamine pH 8.5 for 2 h "
        "before use in chromatography."
    ),
    wash_protocol=(
        "Neutralize to pH 7.0, then 10 CV water, then 5 CV 0.5 M NaCl, "
        "then store in acetate buffer pH 4.0 at 4 C until use."
    ),
    storage_recommendation=(
        "Activated beads: sodium acetate pH 4.0 at 4 C; use within 24-48 h. "
        "ECH reagent: sealed bottles at 4-8 C under nitrogen."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 12. DVS activation (reagent_profiles["dvs_activation"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["dvs_activation"] = MechanismDescriptor(
    reagent_key="dvs_activation",
    display_name="Divinyl Sulfone (OH activation — introduces vinyl sulfone)",
    overall_equation_latex=(
        r"\text{Agarose-OH} + \text{CH}_2\text{=CH-SO}_2\text{-CH=CH}_2 "
        r"\xrightarrow{\text{NaOH, pH 11}} "
        r"\text{Agarose-O-CH}_2\text{CH}_2\text{-SO}_2\text{-CH=CH}_2"
    ),
    mechanism_type="michael_addition",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Oxa-Michael addition: agarose alkoxide (-O-) adds to one "
                "vinyl group of DVS at pH 11-12, forming a stable ether-sulfone "
                "tether with a free pendant vinyl sulfone group. Only one vinyl "
                "group reacts — the other remains available for ligand coupling."
            ),
            bond_formed="C-O ether + C-C (Michael adduct at beta-carbon of vinyl)",
            bond_broken="C=C (one vinyl group of DVS)",
            equation_latex=(
                r"\text{Agarose-O}^- + \text{CH}_2\text{=CH-SO}_2\text{-CH=CH}_2 "
                r"\rightarrow "
                r"\text{Agarose-O-CH}_2\text{CH}_2\text{-SO}_2\text{-CH=CH}_2"
            ),
            conditions="0.5 M NaOH + 0.5 M Na2CO3, pH 11-12, 25 C, 1-2 h",
            notes=(
                "Monoactivation is favoured by using dilute DVS (1-10% v/v) "
                "in excess polymer (agarose). Bis-activation (both vinyl groups "
                "reacting with the same chain = waste) is minimized by this "
                "stoichiometry."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="irreversible",
    critical_notes=[
        "TOXIC AND LACHRYMATORY — handle in certified fume hood. See 'dvs' "
        "mechanism entry for full hazard information.",
        "Vinyl sulfone groups hydrolyze in aqueous buffer: t1/2 ~ 2-4 h at "
        "pH 7.4. Use activated beads promptly or store at pH 4-5.",
        "Vinyl sulfone reacts with -SH (thiol) faster than -NH2 (amine) at "
        "neutral pH — preferred for thiol-containing ligands at pH 6-7, "
        "and for amine coupling at pH 9-10.",
        "DVS-activated beads are more hydrolytically stable than ECH-epoxide "
        "beads — longer working window (2-4 h vs 12-48 h).",
    ],
    buffer_system="0.5 M NaOH / 0.5 M Na2CO3 (activation), sodium carbonate pH 9-10 (coupling)",
    quench_recommendation=(
        "Unreacted VS groups hydrolyze spontaneously. For immediate quench: "
        "10 mM 2-mercaptoethanol for 30 min at pH 7.0."
    ),
    wash_protocol=(
        "Neutralize to pH 7.0, then 10 CV water, then 5 CV 0.5 M NaCl, "
        "then adjust to coupling pH (9-10) or store at pH 4-5."
    ),
    storage_recommendation=(
        "Activated beads: 50 mM sodium acetate pH 4.5 at 4 C; use within 4 h. "
        "DVS reagent: 4-8 C under nitrogen, amber glass."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 13. DEAE coupling (reagent_profiles["deae_coupling"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["deae_coupling"] = MechanismDescriptor(
    reagent_key="deae_coupling",
    display_name="DEAE coupling (epoxide + amine nucleophile)",
    overall_equation_latex=(
        r"\text{Agarose-epoxide} + \text{HN(C}_2\text{H}_5\text{)}_2"
        r"\text{(CH}_2\text{)}_2\text{NH}_2 "
        r"\xrightarrow{\text{pH 10-11}} "
        r"\text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{-NH-(CH}_2\text{)}_2"
        r"\text{-N(C}_2\text{H}_5\text{)}_2"
    ),
    mechanism_type="epoxide_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "The primary amine of 2-(diethylamino)ethylamine (DEAE reagent) "
                "acts as the nucleophile attacking the activated epoxide ring "
                "on the ECH-activated agarose bead via SN2 mechanism at "
                "pH 10-11. The ring opens to form a stable secondary amine "
                "linkage (agarose-O-CH2-CHOH-CH2-NH-R), installing the DEAE "
                "diethylaminoethyl group on the surface."
            ),
            bond_formed="C-N (secondary amine, agarose-CH2-NH-DEAE)",
            bond_broken="C-O (epoxide ring)",
            equation_latex=(
                r"\text{Agarose-CH(O)CH}_2 + \text{H}_2\text{N-R}_{DEAE} "
                r"\rightarrow "
                r"\text{Agarose-CH(OH)-CH}_2\text{-NH-R}_{DEAE}"
            ),
            conditions="pH 10-11, 25 C, 4-16 h",
            notes=(
                "Secondary amine linkage is permanently covalent and stable "
                "over pH 2-14. The tertiary amine of DEAE (pKa ~11.5) provides "
                "weak anion exchange capacity at pH < 9."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="irreversible",
    critical_notes=[
        "Requires prior ECH or BDGE activation of agarose -OH groups. "
        "Ensure epoxide density is adequate before coupling.",
        "DEAE is a weak anion exchanger — fully protonated below pH 9, "
        "partially deprotonated above pH 9. Binding capacity is pH-dependent.",
        "DEAE reagent (2-diethylaminoethylamine) is a volatile liquid with "
        "amine odour — handle in fume hood.",
        "Quench remaining epoxides with ethanolamine after DEAE coupling "
        "to prevent non-specific amine binding during chromatography.",
    ],
    buffer_system="0.1 M NaOH / sodium carbonate pH 10-11",
    quench_recommendation="1 M ethanolamine pH 8.5 for 2 h to quench remaining epoxides.",
    wash_protocol="10 CV 0.5 M NaCl, then 10 CV 20 mM Tris-HCl pH 8.0.",
    storage_recommendation="Coupled beads: PBS pH 7.4 + 0.02% NaN3 at 4 C.",
)

# ───────────────────────────────────────────────────────────────────────────
# 14. IDA coupling (reagent_profiles["ida_coupling"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["ida_coupling"] = MechanismDescriptor(
    reagent_key="ida_coupling",
    display_name="IDA coupling (epoxide + iminodiacetic acid)",
    overall_equation_latex=(
        r"\text{Agarose-epoxide} + \text{NH(CH}_2\text{COOH)}_2 "
        r"\xrightarrow{\text{pH 10-11}} "
        r"\text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{-N(CH}_2\text{COOH)}_2"
    ),
    mechanism_type="epoxide_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "The secondary amine of iminodiacetic acid (IDA, "
                "NH(CH2COOH)2) attacks the activated epoxide ring on the "
                "ECH-activated agarose bead at pH 10-11, opening the ring "
                "to form a tertiary amine linkage. The two acetic acid arms "
                "of IDA remain free and provide bidentate (plus the amine) "
                "tridentate metal chelation geometry for IMAC."
            ),
            bond_formed="C-N (tertiary amine, agarose-CH2-N(IDA))",
            bond_broken="C-O (epoxide ring)",
            equation_latex=(
                r"\text{Agarose-CH(O)CH}_2 + \text{HN(CH}_2\text{COOH)}_2 "
                r"\rightarrow "
                r"\text{Agarose-CH(OH)-CH}_2\text{-N(CH}_2\text{COOH)}_2"
            ),
            conditions="pH 10-11, 25 C, 6-16 h",
            notes=(
                "IDA provides a tridentate chelation geometry (N + 2xO-). "
                "Metal loading (Ni2+, Co2+, Cu2+, Zn2+) is done as a "
                "separate step after IDA coupling and washing."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="irreversible",
    critical_notes=[
        "IDA coupling requires prior ECH activation to introduce epoxides. "
        "Optimal epoxide density: 20-40 umol/mL settled gel.",
        "IDA has lower metal affinity than NTA (tridentate vs tetradentate). "
        "May show more non-specific binding and higher metal leaching than NTA.",
        "Quench remaining epoxides with ethanolamine BEFORE metal charging to "
        "prevent non-specific coordination.",
        "IDA must be charged with metal ion (separate nickel_charging step) "
        "before IMAC application. Store uncharged at 4 C.",
    ],
    buffer_system="0.1 M NaHCO3/Na2CO3 pH 10-11",
    quench_recommendation="1 M ethanolamine pH 8.5 for 2 h after IDA coupling.",
    wash_protocol=(
        "10 CV water, then 10 CV 0.5 M NaCl, then 10 CV 20 mM sodium "
        "phosphate pH 7.4 before metal charging."
    ),
    storage_recommendation="IDA-beads: PBS pH 7.4 + 0.02% NaN3 at 4 C (uncharged); 12 months.",
)

# ───────────────────────────────────────────────────────────────────────────
# 15. Phenyl coupling (reagent_profiles["phenyl_coupling"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["phenyl_coupling"] = MechanismDescriptor(
    reagent_key="phenyl_coupling",
    display_name="Phenyl coupling (aniline-epoxide, HIC ligand)",
    overall_equation_latex=(
        r"\text{Agarose-epoxide} + \text{C}_6\text{H}_5\text{NH}_2 "
        r"\xrightarrow{\text{pH 9-10}} "
        r"\text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{-NH-C}_6\text{H}_5"
    ),
    mechanism_type="epoxide_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "The aromatic amine of aniline (phenylamine, a weaker "
                "nucleophile than aliphatic amines due to resonance delocalisation) "
                "attacks the activated epoxide ring at pH 9-10. The ring opens "
                "to install a phenyl group tethered through a secondary amine "
                "linkage — the benzene ring provides the hydrophobic surface "
                "for HIC (hydrophobic interaction chromatography)."
            ),
            bond_formed="C-N (secondary amine, agarose-CH2-NH-phenyl)",
            bond_broken="C-O (epoxide ring)",
            equation_latex=(
                r"\text{Agarose-CH(O)CH}_2 + \text{H}_2\text{N-Ph} "
                r"\rightarrow "
                r"\text{Agarose-CH(OH)-CH}_2\text{-NH-Ph}"
            ),
            conditions="pH 9-9.5, 25 C, 4-16 h",
            notes=(
                "Aniline is a weaker nucleophile — coupling efficiency is lower "
                "than aliphatic amines. Use excess aniline (10-50 mM) to drive "
                "conversion. k_forward ~ 3e-5 m3/(mol*s) (estimated)."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="irreversible",
    critical_notes=[
        "Aniline is TOXIC (carcinogenic, GHS Cat 2). Handle in fume hood with "
        "gloves and eye protection. Thoroughly wash beads before use.",
        "Phenyl HIC beads bind proteins at high salt, elute at low salt. "
        "Binding strength: Butyl < Phenyl < Octyl.",
        "Quench unreacted epoxides with ethanolamine after phenyl coupling.",
        "Aniline may partially hydrolyze off the support at strongly alkaline "
        "pH (>12) over extended periods — avoid harsh alkaline CIP protocols.",
    ],
    buffer_system="0.1 M sodium carbonate pH 9-10",
    quench_recommendation="1 M ethanolamine pH 8.5 for 2 h.",
    wash_protocol="10 CV 0.5 M NaCl, then 10 CV PBS pH 7.4.",
    storage_recommendation="Phenyl-beads: PBS pH 7.4 + 20% ethanol at 4 C.",
)

# ───────────────────────────────────────────────────────────────────────────
# 16. SP coupling (reagent_profiles["sp_coupling"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["sp_coupling"] = MechanismDescriptor(
    reagent_key="sp_coupling",
    display_name="Sulfopropyl coupling (sultone ring-opening, strong CEX)",
    overall_equation_latex=(
        r"\text{Agarose-epoxide} + \underbrace{\text{CH}_2\text{-CH}_2"
        r"\text{-CH}_2\text{-SO}_3^{-}}_{\text{1,3-propane sultone}} "
        r"\xrightarrow{\text{pH 10-11}} "
        r"\text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{-CH}_2\text{-CH}_2"
        r"\text{-CH}_2\text{-SO}_3^-"
    ),
    mechanism_type="sultone_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Under alkaline conditions (pH 10-11), the alkoxide nucleophile "
                "on the ECH-activated agarose or the agarose -O- opens the "
                "sultone ring (1,3-propane sultone, a cyclic sulfonate ester) "
                "via SN2 attack, producing an agarose-tethered propylsulfonate "
                "group (-CH2CH2CH2-SO3-). The sulfonate group is permanently "
                "negatively charged at all pH — providing strong cation exchange "
                "capacity independent of pH."
            ),
            bond_formed="C-O ether (from agarose alkoxide) + C-S-O3- propylsulfonate",
            bond_broken="C-O (sultone ring ester bond)",
            equation_latex=(
                r"\text{Agarose-O}^- + \text{sultone} \rightarrow "
                r"\text{Agarose-O-CH}_2\text{CH}_2\text{CH}_2\text{-SO}_3^-"
            ),
            conditions="pH 10-11, 25 C, 4-16 h",
            notes=(
                "1,3-Propane sultone is TOXIC — handle in fume hood. "
                "Propylsulfonate is a permanent charge at all pH (pKa ~ -1) — "
                "this distinguishes SP (strong CEX) from CM (weak CEX)."
            ),
        ),
    ],
    byproducts=[],
    reversibility="irreversible",
    critical_notes=[
        "1,3-Propane sultone is classified as a probable carcinogen (IARC 2B). "
        "Handle in certified fume hood with nitrile gloves.",
        "Sulfopropyl (SP) is a STRONG cation exchanger: permanently charged at "
        "all pH. Binding capacity is pH-independent between pH 2-12.",
        "Strong CEX requires high salt (0.5-1 M NaCl) for protein elution — "
        "ensure downstream compatibility.",
        "The sp_coupling profile uses epoxide-activated surface. An alternative "
        "route reacts the sultone directly with agarose -OH under alkaline "
        "conditions (direct etherification) without prior ECH activation.",
    ],
    buffer_system="0.5 M NaOH / Na2CO3 pH 10-11",
    quench_recommendation=(
        "Unreacted sultone is hydrolyzed by water at pH > 8. Wash with "
        "5 CV 50 mM NaOH then 10 CV water to ensure complete hydrolysis."
    ),
    wash_protocol="10 CV water, then 10 CV 0.5 M NaCl, then 10 CV 20 mM acetate pH 5.0.",
    storage_recommendation="SP-beads: 20 mM acetate pH 5.0 + 0.02% NaN3 at 4 C.",
)

# ───────────────────────────────────────────────────────────────────────────
# 17. Protein A coupling (reagent_profiles["protein_a_coupling"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["protein_a_coupling"] = MechanismDescriptor(
    reagent_key="protein_a_coupling",
    display_name="Protein A coupling (random amine coupling to epoxide)",
    overall_equation_latex=(
        r"\text{Agarose-epoxide} + \text{Protein A (Lys-NH}_2\text{)} "
        r"\xrightarrow{\text{pH 8-9, 4 C}} "
        r"\text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{-NH-Lys-Protein A}"
    ),
    mechanism_type="epoxide_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Lysine epsilon-amino groups (-NH2, pKa ~10.5) on Protein A "
                "act as nucleophiles attacking the epoxide-activated agarose "
                "surface. At pH 8-9, a fraction of lysine amines are "
                "deprotonated and reactive. The ring-opening SN2 reaction "
                "forms a stable secondary amine linkage between the protein "
                "and the agarose matrix. Random orientation — Protein A is "
                "coupled through any accessible lysine residue."
            ),
            bond_formed="C-N (secondary amine, Lys-CH2-NH-agarose)",
            bond_broken="C-O (epoxide ring)",
            equation_latex=(
                r"\text{Agarose-CH(O)CH}_2 + \text{Protein A-Lys-NH}_2 "
                r"\rightarrow "
                r"\text{Agarose-CH(OH)-CH}_2\text{-NH-Lys-Protein A}"
            ),
            conditions="pH 8.5-9.5, 4 C, 16-24 h",
            notes=(
                "Coupling at 4 C preserves Protein A folding. Random coupling "
                "orientation reduces IgG binding capacity vs oriented coupling. "
                "activity_retention ~ 0.60 (random) vs 0.80 (C-terminal Cys). "
                "Use DADPA spacer arm for improved capacity."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="irreversible",
    critical_notes=[
        "Random lysine coupling leads to ~40% reduced IgG binding capacity "
        "compared to site-directed Protein A-Cys oriented coupling. "
        "Use protein_a_cys_coupling for maximum activity.",
        "Protein A must be in coupling buffer (pH 8.5-9.5) WITHOUT primary "
        "amine competitors (no Tris, glycine, or BSA).",
        "Couple at 4 C to minimize protein denaturation. Avoid PBS if it "
        "contains amine additives.",
        "Characterize IgG binding capacity before column use. Validate "
        "Protein A leaching per regulatory guidelines (ICH Q6B).",
        "Protein A is derived from Staphylococcus aureus — handle as "
        "biohazardous material.",
    ],
    buffer_system="0.1 M NaHCO3 or PBS pH 8.5-9.0 (no amine additives)",
    quench_recommendation=(
        "1 M ethanolamine pH 8.5 for 2 h at 4 C; wash 5 CV PBS pH 7.4."
    ),
    wash_protocol=(
        "5 CV 1 M NaCl in PBS pH 7.4 (remove non-covalently bound protein), "
        "then 5 CV PBS pH 7.4, then 5 CV 0.1 M glycine-HCl pH 2.7 (strip "
        "any non-covalently adsorbed Protein A), then 5 CV PBS pH 7.4."
    ),
    storage_recommendation=(
        "Protein A beads: PBS pH 7.4 + 0.02% NaN3 at 4 C. "
        "Avoid freeze-thaw. Shelf life: 12-24 months with regular validation."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 18. Ethanolamine quench (reagent_profiles["ethanolamine_quench"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["ethanolamine_quench"] = MechanismDescriptor(
    reagent_key="ethanolamine_quench",
    display_name="Ethanolamine quench (blocks remaining epoxides)",
    overall_equation_latex=(
        r"\text{Agarose-epoxide} + \text{HOCH}_2\text{CH}_2\text{NH}_2 "
        r"\rightarrow "
        r"\text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{-NH-CH}_2\text{CH}_2\text{OH}"
    ),
    mechanism_type="epoxide_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Ethanolamine primary amine (pKa ~9.5) reacts with remaining "
                "activated epoxide groups on the agarose bead surface via SN2 "
                "ring-opening, capping them as stable secondary amine-hydroxyl "
                "adducts. This prevents residual epoxides from reacting "
                "non-specifically with proteins during chromatography."
            ),
            bond_formed="C-N (secondary amine quench product)",
            bond_broken="C-O (epoxide ring)",
            equation_latex=(
                r"\text{Agarose-CH(O)CH}_2 + \text{H}_2\text{N-CH}_2\text{CH}_2\text{OH} "
                r"\rightarrow "
                r"\text{Agarose-CH(OH)-CH}_2\text{-NH-CH}_2\text{CH}_2\text{OH}"
            ),
            conditions="1 M ethanolamine pH 8.0-8.5, 25 C, 2-4 h",
            notes=(
                "Standard quench in Sepharose CNBr and ECH-activated resin "
                "protocols. Ethanolamine is non-toxic and water-soluble. "
                "The resulting -NHCH2CH2OH cap is hydrophilic and minimizes "
                "non-specific protein adsorption."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="irreversible",
    critical_notes=[
        "Quenching is a mandatory step after ligand coupling on epoxide-activated "
        "beads. Un-quenched epoxides will react non-specifically with target "
        "proteins during chromatography, reducing purity.",
        "Use 1 M ethanolamine (excess) at pH 8.0-8.5 to ensure complete quench. "
        "Incubate for minimum 2 h.",
        "Wash thoroughly after ethanolamine quench to remove all traces.",
        "Ethanolamine is mildly irritating — standard PPE sufficient.",
    ],
    buffer_system="1 M ethanolamine adjusted to pH 8.0-8.5",
    quench_recommendation="This step IS the quench — apply after all desired ligand coupling is complete.",
    wash_protocol="5 CV water, then 10 CV PBS pH 7.4.",
    storage_recommendation="Quenched beads: PBS pH 7.4 + 0.02% NaN3 at 4 C.",
)

# ───────────────────────────────────────────────────────────────────────────
# 19. Mercaptoethanol quench (reagent_profiles["mercaptoethanol_quench"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["mercaptoethanol_quench"] = MechanismDescriptor(
    reagent_key="mercaptoethanol_quench",
    display_name="2-Mercaptoethanol quench (thiol-VS Michael addition, blocks VS groups)",
    overall_equation_latex=(
        r"\text{Agarose-CH=CH-SO}_2 + \text{HOCH}_2\text{CH}_2\text{SH} "
        r"\rightarrow "
        r"\text{Agarose-CH}_2\text{CH}_2\text{-SO}_2\text{-CH}_2\text{CH}_2"
        r"\text{-S-CH}_2\text{CH}_2\text{OH}"
    ),
    mechanism_type="michael_addition",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Thiol-Michael addition (aza-Michael with S nucleophile): "
                "the thiol group (-SH) of 2-mercaptoethanol acts as a soft "
                "nucleophile and adds to the beta-carbon of the pendant vinyl "
                "sulfone group on DVS-activated agarose. The reaction is fast "
                "at neutral pH (pH 6-7) due to the high nucleophilicity of "
                "thiolate anion (-S-). Forms a stable thioether bond."
            ),
            bond_formed="C-S thioether (Michael adduct)",
            bond_broken="C=C (vinyl sulfone beta-carbon)",
            equation_latex=(
                r"\text{Agarose-SO}_2\text{-CH=CH}_2 + "
                r"\text{HOCH}_2\text{CH}_2\text{S}^- "
                r"\rightarrow "
                r"\text{Agarose-SO}_2\text{-CH}_2\text{CH}_2\text{-S-CH}_2"
                r"\text{CH}_2\text{OH}"
            ),
            conditions="pH 6.0-7.0, 25 C, 1-2 h",
            notes=(
                "Thiol-VS is faster than amine-VS at pH 6-7 (rate ~ 5e-3 "
                "m3/(mol*s)). Use at pH 6-7 to preferentially quench VS over "
                "amine groups. 2-Mercaptoethanol: malodorous and toxic — use "
                "in fume hood."
            ),
        ),
    ],
    byproducts=[],
    reversibility="irreversible",
    critical_notes=[
        "2-Mercaptoethanol is TOXIC and malodorous. Handle in certified fume "
        "hood with nitrile gloves. Avoid skin contact.",
        "Thioether bonds are stable at pH 1-14 and resistant to "
        "denaturing agents (urea, guanidine). Permanent quench.",
        "Apply at pH 6-7 where thiolate reactivity is high but amine "
        "reactivity (pKa-dependent) is low — selectively quenches VS without "
        "blocking coupled amine ligands.",
        "This quench is used after DVS-activation or after VS-group-containing "
        "reagent coupling.",
    ],
    buffer_system="50 mM sodium phosphate pH 6.5",
    quench_recommendation="This step IS the quench for VS groups.",
    wash_protocol="10 CV PBS pH 7.4 to remove all mercaptoethanol.",
    storage_recommendation=(
        "Store 2-mercaptoethanol at 4 C, sealed, in fume hood-accessible "
        "location. Quenched beads: PBS pH 7.4 + 0.02% NaN3 at 4 C."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 20. NaBH4 quench (reagent_profiles["nabh4_quench"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["nabh4_quench"] = MechanismDescriptor(
    reagent_key="nabh4_quench",
    display_name="NaBH4 quench (reductive amination — stabilizes Schiff base)",
    overall_equation_latex=(
        r"\text{Polymer-N=CH-R} + \text{NaBH}_4 "
        r"\xrightarrow{\text{PBS pH 7-8}} "
        r"\text{Polymer-NH-CH}_2\text{-R} + \text{NaBO}_2\text{/H}_2"
    ),
    mechanism_type="reductive_amination",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Hydride (H-) transfer from sodium borohydride (NaBH4) to the "
                "electrophilic carbon of the Schiff base imine (C=N) produces "
                "the corresponding amine via a two-electron reduction. "
                "The unstable, hydrolytically reversible imine (C=N, Schiff base) "
                "formed by glutaraldehyde-amine condensation is converted to a "
                "stable, permanent secondary amine (C-NH) bond."
            ),
            bond_formed="C-H and N-H (secondary amine from imine reduction)",
            bond_broken="C=N (imine / Schiff base reduced to C-NH)",
            equation_latex=(
                r"\text{R-CH=N-R'} + \text{H}^- \text{(from NaBH}_4\text{)} "
                r"\rightarrow \text{R-CH}_2\text{-NH-R'}"
            ),
            conditions="10 mM NaBH4 in PBS pH 7-8, 25 C, 30 min",
            notes=(
                "NaBH4 also reduces free aldehyde groups to alcohols. "
                "Use minimum effective NaBH4 concentration (10-20 mM) to "
                "avoid unwanted reduction of other functional groups. "
                "H2 gas is evolved — ensure adequate ventilation."
            ),
        ),
    ],
    byproducts=["NaBO2", "H2"],
    reversibility="irreversible",
    critical_notes=[
        "FLAMMABLE / CORROSIVE — NaBH4 reacts with water to evolve H2 gas. "
        "Add slowly to aqueous solution in well-ventilated area. Avoid flames.",
        "After NaBH4 reduction, the secondary amine C-NH bond is stable at "
        "pH 1-14 and cannot be reversed by mild acid hydrolysis — unlike the "
        "precursor Schiff base.",
        "Apply this reduction step immediately after glutaraldehyde crosslinking "
        "or after reductive amination ligand coupling while Schiff bases are "
        "still intact.",
        "NaBH4 concentration must not exceed 50 mM — excess borohydride may "
        "reduce other functional groups (esters, NHS esters if present).",
    ],
    buffer_system="PBS pH 7.4 (freshly prepared NaBH4 solution)",
    quench_recommendation=(
        "NaBH4 is quenched by water and acidic pH. Wash 5 CV PBS pH 7.4 "
        "immediately after reaction to remove borate byproducts."
    ),
    wash_protocol="5 CV PBS pH 7.4, then 5 CV 0.5 M NaCl/PBS, then 5 CV PBS pH 7.4.",
    storage_recommendation=(
        "Store NaBH4 as anhydrous powder at RT in desiccated container. "
        "Prepare aqueous solution fresh immediately before use."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 21. Nickel charging (reagent_profiles["nickel_charging"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["nickel_charging"] = MechanismDescriptor(
    reagent_key="nickel_charging",
    display_name="Nickel(II) charging of IDA/NTA chelator",
    overall_equation_latex=(
        r"\text{Agarose-IDA/NTA} + \text{Ni}^{2+} "
        r"\rightleftharpoons "
        r"\text{Agarose-IDA/NTA}\cdot\text{Ni}^{2+} \quad "
        r"(K_a \approx 10^{11} \text{ M}^{-1}\text{ for NTA-Ni})"
    ),
    mechanism_type="metal_chelation",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "Ni2+ ions in solution coordinate to the tridentate (IDA) or "
                "tetradentate (NTA) chelator groups immobilized on the agarose "
                "surface. For IDA: coordination via N (amine) + 2x carboxylate "
                "O. For NTA: coordination via N + 3x carboxylate O. The "
                "remaining coordination sites (2 for IDA, 1 for NTA) are "
                "available for histidine imidazole coordination from the "
                "His-tagged target protein."
            ),
            bond_formed="Ni-N and Ni-O coordination bonds (dative, chelation)",
            bond_broken="Ni-H2O (displacement of water ligands from Ni2+ hexaaquo complex)",
            equation_latex=(
                r"[\text{Ni(H}_2\text{O)}_6]^{2+} + \text{IDA}_{(s)} "
                r"\rightleftharpoons "
                r"\text{IDA-Ni}^{2+}_{(s)} + 3\,\text{H}_2\text{O}"
            ),
            conditions="50 mM NiSO4 in PBS pH 7.0-8.0, 25 C, 30 min",
            notes=(
                "log Ka(NTA-Ni) ~ 11.5; log Ka(IDA-Ni) ~ 8.5 at 25 C "
                "(Martell & Smith Critical Stability Constants). "
                "After loading, wash excess Ni2+ with 5 CV equilibration buffer. "
                "Ni2+ leaches ~0.01-0.1 ppm per cycle — monitor if required."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="equilibrium",
    critical_notes=[
        "Ni2+ is mildly cytotoxic — ensure adequate washing before "
        "biological use. Monitor Ni leaching by ICP-MS if regulated.",
        "Ni2+ loading requires no organic solvent and is carried out at "
        "physiological pH — safe for protein-compatible media.",
        "EDTA (0.1-0.5 M) strips Ni2+ from the column for regeneration "
        "or metal switching (see edta_stripping profile).",
        "Imidazole (0.25-0.5 M) in wash buffer removes weakly-bound "
        "non-His-tagged proteins (competitive elution via imidazole-Ni "
        "coordination).",
        "Store Ni-loaded beads in PBS + 0.02% NaN3 at 4 C. Avoid EDTA-"
        "containing buffers during storage.",
    ],
    buffer_system="50 mM NiSO4 in PBS pH 7.4",
    quench_recommendation=(
        "Wash 5 CV PBS pH 7.4 to remove excess unbound Ni2+."
    ),
    wash_protocol=(
        "5 CV equilibration buffer (20 mM sodium phosphate, 0.5 M NaCl, "
        "pH 7.4), then 5 CV equilibration buffer + 20 mM imidazole "
        "(removes very weak binders)."
    ),
    storage_recommendation=(
        "Ni-IDA/NTA beads: 20 mM sodium phosphate, 0.5 M NaCl, pH 7.4 "
        "+ 0.02% NaN3 at 4 C. Stable 12 months if Ni not stripped."
    ),
)

# ───────────────────────────────────────────────────────────────────────────
# 22. BDGE activation (reagent_profiles["bdge_activation"])
# ───────────────────────────────────────────────────────────────────────────
MECHANISM_REGISTRY["bdge_activation"] = MechanismDescriptor(
    reagent_key="bdge_activation",
    display_name="BDGE activation (bis-epoxide, long-arm 18 A spacer)",
    overall_equation_latex=(
        r"\text{Agarose-OH} + \text{(CH}_2\text{O)CH-CH}_2\text{-O-"
        r"(CH}_2\text{)}_4\text{-O-CH}_2\text{-CH(O)CH}_2 "
        r"\xrightarrow{\text{NaOH, pH 11}} "
        r"\text{Agarose-O-CH}_2\text{CHOH-CH}_2\text{-O-(CH}_2\text{)}_4"
        r"\text{-O-CH}_2\text{-CH(O)-CH}_2"
    ),
    mechanism_type="epoxide_ring_opening",
    steps=[
        ReactionStep(
            step_number=1,
            description=(
                "BDGE (1,4-butanediol diglycidyl ether) is a bis-epoxide "
                "spacer with a 4-carbon butyl chain between two glycidyl "
                "(epoxide) groups. The first epoxide reacts with an agarose "
                "alkoxide (-O-) under alkaline conditions, anchoring the "
                "spacer to the agarose backbone via a stable ether bond and "
                "leaving the second, distal epoxide free for ligand coupling."
            ),
            bond_formed="C-O ether (agarose-O-CH2-CHOH-CH2-O-spacer)",
            bond_broken="C-O (proximal epoxide ring opening)",
            equation_latex=(
                r"\text{Agarose-O}^- + \text{Ep-spacer-Ep} "
                r"\rightarrow "
                r"\text{Agarose-O-CH}_2\text{CHOH-spacer-Ep}"
            ),
            conditions="0.5 M NaOH, pH 11-12, 25 C, 2-4 h",
            notes=(
                "Monoactivation is essential — if both epoxides react with the "
                "same agarose bead, a bis-ether crosslink is formed instead of "
                "a spacer arm. Use dilute BDGE (0.5-2% v/v) with excess "
                "agarose to maximize mono-adduct formation."
            ),
        ),
    ],
    byproducts=["H2O"],
    reversibility="irreversible",
    critical_notes=[
        "BDGE provides an 18 A spacer arm vs ~5 A for ECH — longer spacer "
        "reduces steric hindrance and improves protein ligand coupling efficiency.",
        "Risk of bis-activation (both epoxides reacting with the same chain) "
        "increases at high BDGE concentration — keep BDGE < 5% v/v.",
        "BDGE is an irritant — handle in fume hood with gloves.",
        "Product epoxide density should be validated before ligand coupling "
        "(titration with thiosulfate or Ellman's reagent).",
        "BDGE-activated beads show faster hydrolysis than ECH-epoxide beads "
        "(secondary vs primary epoxide). Use within 24-48 h of activation.",
    ],
    buffer_system="0.5 M NaOH (activation), sodium acetate pH 4.5 (storage after activation)",
    quench_recommendation="Quench remaining epoxides with 1 M ethanolamine pH 8.5 for 2 h.",
    wash_protocol=(
        "Neutralize with HCl to pH 7.0, then 10 CV water, then 5 CV 0.5 M NaCl, "
        "then store in acetate buffer pH 4.5 at 4 C until ligand coupling."
    ),
    storage_recommendation=(
        "BDGE: sealed at RT under nitrogen. "
        "Activated beads: acetate pH 4.5 at 4 C; use within 24 h."
    ),
)


# ═══════════════════════════════════════════════════════════════════════════
#  FALLBACK AUTO-GENERATOR
# ═══════════════════════════════════════════════════════════════════════════


def get_mechanism(reagent_key: str) -> MechanismDescriptor:
    """Return the MechanismDescriptor for *reagent_key*.

    If the key is in MECHANISM_REGISTRY, the manually authored entry is
    returned directly.  Otherwise a minimal fallback descriptor is
    auto-generated from the reagent's existing profile data in
    reagent_library.CROSSLINKERS or
    module2_functionalization.reagent_profiles.REAGENT_PROFILES.

    Parameters
    ----------
    reagent_key:
        Reagent identifier string matching a key in CROSSLINKERS or
        REAGENT_PROFILES.

    Returns
    -------
    MechanismDescriptor
        Either the manually authored entry or an auto-generated minimal one.
    """
    # Fast path — manually authored entry exists.
    if reagent_key in MECHANISM_REGISTRY:
        return MECHANISM_REGISTRY[reagent_key]

    # ── Attempt to build fallback from CROSSLINKERS ──────────────────────
    try:
        from emulsim.reagent_library import CROSSLINKERS

        if reagent_key in CROSSLINKERS:
            profile = CROSSLINKERS[reagent_key]
            fallback_step = ReactionStep(
                step_number=1,
                description=profile.notes[:200] if profile.notes else "No mechanism notes available.",
                bond_formed="(not yet authored)",
                bond_broken="(not yet authored)",
                equation_latex="",
                conditions=(
                    f"T = {profile.T_crosslink_default - 273.15:.0f} C, "
                    f"t = {profile.t_crosslink_default / 3600:.1f} h"
                ),
                notes="Auto-generated fallback from CrosslinkerProfile.notes.",
            )
            return MechanismDescriptor(
                reagent_key=reagent_key,
                display_name=profile.name,
                overall_equation_latex="",
                mechanism_type=profile.mechanism,
                steps=[fallback_step],
                byproducts=[],
                reversibility=(
                    "reversible_ionic" if profile.reversible else "irreversible"
                ),
                critical_notes=[
                    "Mechanism details not yet authored for this reagent.",
                    f"Hazard info: see notes in CrosslinkerProfile for '{reagent_key}'.",
                ],
            )
    except ImportError:
        pass

    # ── Attempt to build fallback from REAGENT_PROFILES ──────────────────
    try:
        from emulsim.module2_functionalization.reagent_profiles import REAGENT_PROFILES

        if reagent_key in REAGENT_PROFILES:
            # Use a distinct variable name from the CROSSLINKERS branch above
            # so mypy does not narrow the type across the two try/except blocks.
            m2_profile = REAGENT_PROFILES[reagent_key]
            note_text = m2_profile.notes[:200] if m2_profile.notes else "No mechanism notes available."
            fallback_step = ReactionStep(
                step_number=1,
                description=note_text,
                bond_formed="(not yet authored)",
                bond_broken="(not yet authored)",
                equation_latex="",
                conditions=(
                    f"pH {m2_profile.ph_optimum:.1f}, "
                    f"T = {m2_profile.temperature_default - 273.15:.0f} C, "
                    f"t = {m2_profile.time_default / 3600:.1f} h"
                ),
                notes="Auto-generated fallback from ReagentProfile.notes.",
            )
            hazard_note = (
                f"Hazard class: {m2_profile.hazard_class}"
                if m2_profile.hazard_class
                else "No hazard classification provided."
            )
            # Derive buffer_system from profile pH (Codex fix 3)
            _ph = getattr(m2_profile, "ph_optimum", 7.4)
            if _ph >= 10.0:
                _buffer = f"0.1 M NaOH/Na2CO3 pH {_ph:.0f}"
            elif _ph <= 5.5:
                _buffer = f"MES buffer pH {_ph:.1f}"
            else:
                _buffer = f"PBS pH {_ph:.1f}"
            return MechanismDescriptor(
                reagent_key=reagent_key,
                display_name=m2_profile.name,
                overall_equation_latex="",
                mechanism_type=m2_profile.chemistry_class if m2_profile.chemistry_class else "unknown",
                steps=[fallback_step],
                byproducts=[],
                reversibility="irreversible",
                critical_notes=[
                    "Mechanism details not yet authored for this reagent.",
                    hazard_note,
                ],
                buffer_system=_buffer,
            )
    except ImportError:
        pass

    # ── Ultimate fallback — unknown reagent ──────────────────────────────
    return MechanismDescriptor(
        reagent_key=reagent_key,
        display_name=reagent_key,
        overall_equation_latex="",
        mechanism_type="unknown",
        steps=[
            ReactionStep(
                step_number=1,
                description="Unknown reagent — no profile or mechanism data found.",
                bond_formed="(unknown)",
                bond_broken="(unknown)",
                equation_latex="",
                conditions="(unknown)",
                notes=f"Reagent key '{reagent_key}' not found in CROSSLINKERS or REAGENT_PROFILES.",
            )
        ],
        byproducts=[],
        reversibility="irreversible",
        critical_notes=[
            f"No data found for reagent key '{reagent_key}'. "
            "Check that the key matches an entry in reagent_library.CROSSLINKERS "
            "or module2_functionalization.reagent_profiles.REAGENT_PROFILES."
        ],
    )
