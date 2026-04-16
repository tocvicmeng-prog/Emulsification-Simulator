"""Protocol generator engine for EmulSim wet-lab protocols.

Generates tailored ProtocolDocument objects from reagent profile data
and user-supplied reaction parameters.

Architecture reference: WN-2 (protocol_generator module spec).
Depends on:
  - protocols.mechanism_data  (get_mechanism, MechanismDescriptor)
  - protocols.protocol_document (ProtocolDocument, ProtocolStep, ReagentRequirement)
  - reagent_library.CROSSLINKERS
  - module2_functionalization.reagent_profiles.REAGENT_PROFILES

DO NOT import streamlit here — this is pure data/logic.
"""

from __future__ import annotations

from datetime import datetime

from emulsim.protocols.mechanism_data import get_mechanism
from emulsim.protocols.protocol_document import (
    ProtocolDocument,
    ProtocolStep,
    ReagentRequirement,
)

# ─── Known molecular weights (g/mol) ──────────────────────────────────────────
_KNOWN_MW: dict[str, float] = {
    "genipin": 226.23,
    "glutaraldehyde": 100.12,
    "edc_nhs": 191.70,          # EDC MW (NHS is co-reagent)
    "pegda_uv": 700.0,
    "tpp": 367.86,
    "epichlorohydrin": 92.52,
    "dvs": 118.15,
    "citric_acid": 192.12,
    "genipin_secondary": 226.23,
    "glutaraldehyde_secondary": 100.12,
    "ech_activation": 92.52,
    "dvs_activation": 118.15,
    "bdge_activation": 202.25,
    "ethanolamine_quench": 61.08,
    "mercaptoethanol_quench": 78.13,
    "nabh4_quench": 37.83,
    "acetic_anhydride_quench": 102.09,
    "nickel_charging": 262.85,  # NiSO4·6H2O
    "cobalt_charging": 237.93,  # CoCl2·6H2O
    "copper_charging": 249.69,  # CuSO4·5H2O
    "zinc_charging": 287.56,    # ZnSO4·7H2O
    "edta_stripping": 372.24,   # EDTA disodium·2H2O
    "tcep_reduction": 286.65,
    "dtt_reduction": 154.25,
    "wash_buffer": 0.0,         # no reagent mass
}

# ─── Hazard → safety warning text ─────────────────────────────────────────────
_HAZARD_WARNINGS: dict[str, str] = {
    "toxic": "WARNING: Toxic reagent. Handle with gloves and in a fume hood.",
    "toxic_carcinogen": (
        "WARNING: CARCINOGENIC. Use in certified fume hood with respiratory "
        "protection. Consult institutional safety officer."
    ),
    "toxic_malodorous": "WARNING: Toxic and malodorous. Handle in fume hood.",
    "flammable_corrosive": (
        "WARNING: Flammable and corrosive. No open flames. Wear face shield."
    ),
    "corrosive_flammable": (
        "WARNING: Flammable and corrosive. No open flames. Wear face shield."
    ),
    "irritant": "CAUTION: Irritant. Wear gloves and safety glasses.",
    "low_hazard": "Standard PPE: lab coat, gloves, safety glasses.",
    "": "Standard PPE: lab coat, gloves, safety glasses.",
}


class ProtocolGenerator:
    """Generates tailored ProtocolDocuments from reagent data + user parameters."""

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def generate(
        self,
        reagent_key: str,
        temperature_K: float,
        time_s: float,
        concentration_mM: float,
        pH: float,
        bead_volume_mL: float = 1.0,
        source: str = "reagent_profiles",  # or "crosslinkers"
    ) -> ProtocolDocument:
        """Generate a complete wet-lab ProtocolDocument.

        Parameters
        ----------
        reagent_key:
            Key into CROSSLINKERS or REAGENT_PROFILES.
        temperature_K:
            Reaction temperature [K].
        time_s:
            Reaction time [s].
        concentration_mM:
            Reagent working concentration [mM].
        pH:
            Reaction pH.
        bead_volume_mL:
            Settled bead volume [mL]. Defaults to 1.0.
        source:
            "crosslinkers" or "reagent_profiles".

        Returns
        -------
        ProtocolDocument
        """
        # 1. Load profile and mechanism
        profile = self._load_profile(reagent_key, source)
        mechanism = self._load_mechanism(reagent_key)

        # 2. Derived quantities
        total_volume_mL = bead_volume_mL * 10.0
        temp_C = temperature_K - 273.15
        time_str = ProtocolDocument._format_time(time_s)
        mw = self._get_mw(reagent_key, profile, source)
        mass_mg = self._calc_mass_mg(concentration_mM, mw, total_volume_mL)

        # 3. Determine buffer
        buffer = mechanism.buffer_system or f"PBS pH {pH:.1f}"

        # 4. Build output sections
        reagent_table = self._build_reagent_table(
            reagent_key, profile, source, mechanism,
            mass_mg, mw, total_volume_mL, pH,
        )
        procedure_steps = self._build_procedure_steps(
            reagent_key, profile, source, mechanism,
            bead_volume_mL, total_volume_mL, temp_C, time_str, pH, buffer,
        )
        qc_targets = self._build_qc_targets(mechanism, profile, source)
        safety_warnings = self._build_safety_warnings(
            profile, source, pH, reagent_key=reagent_key,
        )

        return ProtocolDocument(
            title=mechanism.display_name,
            reagent_key=reagent_key,
            generation_timestamp=datetime.now().isoformat(timespec="seconds"),
            user_temperature_K=temperature_K,
            user_time_s=time_s,
            user_concentration_mM=concentration_mM,
            user_pH=pH,
            user_bead_volume_mL=bead_volume_mL,
            reagent_table=reagent_table,
            procedure_steps=procedure_steps,
            qc_targets=qc_targets,
            safety_warnings=safety_warnings,
            mechanism_summary=mechanism.overall_equation_latex,
        )

    # ------------------------------------------------------------------
    # Profile / mechanism loading
    # ------------------------------------------------------------------

    def _load_profile(self, reagent_key: str, source: str) -> object:
        if source == "crosslinkers":
            from emulsim.reagent_library import CROSSLINKERS
            if reagent_key not in CROSSLINKERS:
                raise ValueError(
                    f"Reagent key {reagent_key!r} not found in CROSSLINKERS. "
                    f"Available keys: {sorted(CROSSLINKERS)}"
                )
            return CROSSLINKERS[reagent_key]
        else:
            from emulsim.module2_functionalization.reagent_profiles import REAGENT_PROFILES
            if reagent_key not in REAGENT_PROFILES:
                raise ValueError(
                    f"Reagent key {reagent_key!r} not found in REAGENT_PROFILES. "
                    f"Available keys: {sorted(REAGENT_PROFILES)}"
                )
            return REAGENT_PROFILES[reagent_key]

    def _load_mechanism(self, reagent_key: str):
        return get_mechanism(reagent_key)

    # ------------------------------------------------------------------
    # MW and mass calculations
    # ------------------------------------------------------------------

    def _get_mw(self, reagent_key: str, profile: object, source: str) -> float:
        """Return MW in g/mol. Returns 0.0 for reagents with no mass (e.g. wash_buffer)."""
        # Known override dict first — 0.0 is a valid sentinel meaning "no reagent mass"
        if reagent_key in _KNOWN_MW:
            return _KNOWN_MW[reagent_key]
        # For reagent_profiles, try ligand_mw field
        if source == "reagent_profiles":
            lmw = getattr(profile, "ligand_mw", 0.0)
            if lmw and lmw > 0.0:
                return lmw
        # Default fallback — return 200.0 g/mol with no crash
        return 200.0

    @staticmethod
    def _calc_mass_mg(concentration_mM: float, mw: float, total_volume_mL: float) -> float:
        """mass [mg] = conc [mM] * MW [g/mol] * volume [mL] / 1000"""
        return concentration_mM * mw * total_volume_mL / 1000.0

    # ------------------------------------------------------------------
    # Reagent table
    # ------------------------------------------------------------------

    def _build_reagent_table(
        self,
        reagent_key: str,
        profile: object,
        source: str,
        mechanism,
        mass_mg: float,
        mw: float,
        total_volume_mL: float,
        pH: float,
    ) -> list[ReagentRequirement]:
        rows: list[ReagentRequirement] = []

        name = getattr(profile, "name", reagent_key)
        cas = getattr(profile, "cas", "N/A")
        mtype = mechanism.mechanism_type

        # Primary reagent (skip if zero MW = no reagent mass, e.g. wash_buffer)
        if mw > 0.0:
            rows.append(ReagentRequirement(
                name=name,
                cas=cas,
                amount=f"{mass_mg:.1f} mg",
                grade=">=95%",
                role=self._infer_role(reagent_key, profile, source),
            ))

        # For EDC/NHS: add NHS as co-reagent (equimolar to EDC)
        if reagent_key == "edc_nhs":
            nhs_mw = 115.09
            nhs_mass_mg = (mass_mg / _KNOWN_MW.get("edc_nhs", 191.70)) * nhs_mw
            rows.append(ReagentRequirement(
                name="N-Hydroxysuccinimide (NHS)",
                cas="6066-82-6",
                amount=f"{nhs_mass_mg:.1f} mg",
                grade=">=98%",
                role="NHS ester stabiliser (co-reagent)",
            ))

        # Buffer
        buffer_name = mechanism.buffer_system or f"PBS pH {pH:.1f}"
        rows.append(ReagentRequirement(
            name=buffer_name,
            cas="N/A",
            amount=f"{total_volume_mL:.1f} mL",
            grade="AR grade",
            role="reaction buffer",
        ))

        # UV photoinitiator for radical-chain / UV crosslinking
        if mtype == "radical_chain":
            rows.append(ReagentRequirement(
                name="Irgacure 2959 (photoinitiator)",
                cas="106797-53-9",
                amount=f"{0.05 * total_volume_mL / 100 * 1000:.1f} mg (0.05% w/v)",
                grade=">=98%",
                role="photoinitiator",
            ))

        # NaOH for strongly alkaline reactions
        if pH > 10.0:
            rows.append(ReagentRequirement(
                name="NaOH (10 M stock)",
                cas="1310-73-2",
                amount="As needed to adjust pH",
                grade="AR grade",
                role="pH adjustment",
            ))

        # Ethanol co-solvent for sparingly soluble reagents
        if reagent_key in {"dvs", "dvs_activation", "epichlorohydrin", "ech_activation"}:
            rows.append(ReagentRequirement(
                name="Ethanol (absolute)",
                cas="64-17-5",
                amount=f"{total_volume_mL * 0.1:.1f} mL (10% v/v)",
                grade="HPLC grade",
                role="co-solvent",
            ))

        return rows

    @staticmethod
    def _infer_role(reagent_key: str, profile: object, source: str) -> str:
        if source == "crosslinkers":
            mech = getattr(profile, "mechanism", "")
            if "uv" in mech or "uv" in reagent_key:
                return "UV crosslinker"
            if "ionic" in mech:
                return "ionic crosslinker"
            return "chemical crosslinker"
        fmode = getattr(profile, "functional_mode", "")
        if fmode:
            return fmode.replace("_", " ")
        rtype = getattr(profile, "reaction_type", "")
        return rtype if rtype else "reagent"

    # ------------------------------------------------------------------
    # Procedure steps
    # ------------------------------------------------------------------

    def _build_procedure_steps(
        self,
        reagent_key: str,
        profile: object,
        source: str,
        mechanism,
        bead_volume_mL: float,
        total_volume_mL: float,
        temp_C: float,
        time_str: str,
        pH: float,
        buffer: str,
    ) -> list[ProtocolStep]:
        steps: list[ProtocolStep] = []
        mtype = mechanism.mechanism_type
        temp_str = f"{temp_C:.1f} deg C"
        step_n = 1

        # ── Step 1: Prepare solution ───────────────────────────────────
        name = getattr(profile, "name", reagent_key)
        if mtype == "radical_chain":
            dissolve_note = "Dissolve PEGDA in buffer, then add photoinitiator. Protect from ambient light."
            dissolve_action = (
                f"Dissolve {name} in {buffer} to working concentration. "
                "Add Irgacure 2959 to 0.05% w/v. Degas under N2 for 5 min."
            )
        elif mtype == "ionic_interaction":
            dissolve_action = (
                f"Dissolve {name} in ultrapure water at the target concentration. "
                "Filter through 0.22 µm membrane. Check pH."
            )
            dissolve_note = "TPP solution must be fresh; gelation is instantaneous on contact."
        elif reagent_key in {"ech_activation", "dvs_activation", "epichlorohydrin", "dvs"}:
            dissolve_action = (
                f"Prepare {total_volume_mL:.1f} mL of {name} solution in 0.5 M NaOH. "
                "Mix in fume hood."
            )
            dissolve_note = "Work in fume hood. Prepare fresh immediately before use."
        else:
            dissolve_action = (
                f"Dissolve {name} in {buffer} to the working concentration "
                f"({total_volume_mL:.1f} mL total). Mix until clear."
            )
            dissolve_note = mechanism.critical_notes[0] if mechanism.critical_notes else ""
        steps.append(ProtocolStep(
            step_number=step_n,
            action=dissolve_action,
            duration="5–10 min",
            temperature="RT",
            notes=dissolve_note,
        ))
        step_n += 1

        # ── Step 2: Wash beads ─────────────────────────────────────────
        steps.append(ProtocolStep(
            step_number=step_n,
            action=(
                f"Wash {bead_volume_mL:.1f} mL settled beads 2x with "
                f"{total_volume_mL:.1f} mL {buffer}."
            ),
            duration="5 min per wash",
            temperature="RT",
            notes="Remove storage buffer completely before reaction.",
        ))
        step_n += 1

        # ── Step 3: Add reagent ────────────────────────────────────────
        if mtype == "ionic_interaction":
            add_action = (
                f"Add {total_volume_mL:.1f} mL {name} solution dropwise to beads. "
                "Gelation is instantaneous."
            )
            add_note = "Maintain gentle agitation during addition."
        else:
            add_action = (
                f"Add {total_volume_mL:.1f} mL reagent solution to beads. "
                "Ensure complete mixing."
            )
            add_note = ""
        steps.append(ProtocolStep(
            step_number=step_n,
            action=add_action,
            duration="Immediate",
            temperature=temp_str,
            notes=add_note,
        ))
        step_n += 1

        # ── Step 4: React (incubate or UV exposure) ────────────────────
        if mtype == "radical_chain":
            react_action = (
                f"Expose bead suspension to UV light (365 nm, >=10 mW/cm^2) "
                f"for {time_str} with gentle swirling."
            )
            react_note = "Keep sample in ice bath during UV exposure to limit thermal damage."
        elif mtype == "ionic_interaction":
            react_action = (
                "Incubate with gentle end-over-end rotation. "
                "Ionic crosslinks form immediately; extended incubation improves homogeneity."
            )
            react_note = "pH stability: keep pH 4–6 for chitosan-TPP."
        elif mtype in {"reductive_amination"}:
            react_action = (
                f"Incubate at {temp_str} for {time_str} with gentle end-over-end rotation. "
                "Add NaBH3CN (5 mM) after 30 min for Schiff base reduction."
            )
            react_note = "Keep pH 7–8 for optimal imine formation."
        elif mtype == "metal_chelation":
            react_action = (
                f"Incubate bead suspension in metal salt solution at {temp_str} for {time_str}."
            )
            react_note = "Use 50 mM metal sulfate/chloride in binding buffer (pH 8)."
        else:
            conditions_str = ""
            if mechanism.steps:
                conditions_str = mechanism.steps[0].conditions
            react_action = (
                f"Incubate at {temp_str} for {time_str} with gentle end-over-end rotation."
            )
            react_note = conditions_str if conditions_str else ""
        steps.append(ProtocolStep(
            step_number=step_n,
            action=react_action,
            duration=time_str,
            temperature=temp_str,
            notes=react_note,
        ))
        step_n += 1

        # ── Step 5: Quench (if applicable) ────────────────────────────
        quench_rec = mechanism.quench_recommendation
        if quench_rec and mtype not in {"ionic_interaction"}:
            steps.append(ProtocolStep(
                step_number=step_n,
                action=f"Quench: {quench_rec}",
                duration="30 min",
                temperature=temp_str,
                notes="Add quench reagent to stop reaction and cap remaining reactive groups.",
            ))
            step_n += 1

        # ── Step 6: Post-reaction wash ─────────────────────────────────
        wash_rec = mechanism.wash_protocol or (
            f"Wash 5x with {total_volume_mL:.1f} mL PBS"
        )
        steps.append(ProtocolStep(
            step_number=step_n,
            action=wash_rec,
            duration="10 min per wash",
            temperature="RT",
            notes="Remove unreacted reagent completely. Check wash effluent for absence of reagent.",
        ))
        step_n += 1

        # ── Step 7: Storage ────────────────────────────────────────────
        storage = (
            mechanism.storage_recommendation
            or "Store in PBS + 0.02% NaN3 at 4°C"
        )
        steps.append(ProtocolStep(
            step_number=step_n,
            action=f"Storage: {storage}",
            duration="Ongoing",
            temperature="4 deg C",
            notes="Check bead integrity before use in chromatographic application.",
        ))

        return steps

    # ------------------------------------------------------------------
    # QC targets
    # ------------------------------------------------------------------

    def _build_qc_targets(self, mechanism, profile: object, source: str) -> list[str]:
        targets: list[str] = []
        mtype = mechanism.mechanism_type

        if mtype in {
            "nucleophilic_addition", "schiff_base", "carbodiimide_activation",
            "radical_chain", "michael_addition", "ester_bond_formation",
            "epoxide_ring_opening", "sultone_ring_opening",
        }:
            fmode = getattr(profile, "functional_mode", "")
            rtype = getattr(profile, "reaction_type", "") if source == "reagent_profiles" else ""

            if fmode == "crosslinker" or (source == "crosslinkers"):
                targets.append("Target crosslink conversion: >50%")
            elif fmode == "activator" or rtype == "activation":
                targets.append(
                    "Target activation density: check with TNBS assay (amine) or thiol assay"
                )
            elif fmode in {"iex_ligand", "affinity_ligand", "hic_ligand", "imac_chelator"}:
                is_macro = getattr(profile, "is_macromolecule", False)
                activity = getattr(profile, "activity_retention", None)
                if is_macro and activity is not None:
                    targets.append(
                        f"Verify by BCA protein assay; expected retention: {activity:.0%}"
                    )
                else:
                    targets.append("Verify ligand density by appropriate assay")
            elif fmode == "quencher" or rtype == "quenching":
                targets.append(
                    "Verify no residual reactive groups by appropriate assay"
                )
            else:
                targets.append("Verify ligand density by appropriate assay")

        elif mtype == "ionic_interaction":
            targets.append("Target crosslink conversion: >50%")

        elif mtype == "metal_chelation":
            targets.append("Verify metal loading by ICP-MS or colorimetric assay")

        elif mtype == "reductive_amination":
            targets.append("Verify ligand density by appropriate assay")

        else:
            targets.append("Verify reaction completion by appropriate assay")

        # Universal check
        targets.append("Visual inspection: no aggregation or debris")
        return targets

    # ------------------------------------------------------------------
    # Safety warnings
    # ------------------------------------------------------------------

    # Crosslinker-specific hazards (CrosslinkerProfile has no hazard_class field)
    _CROSSLINKER_HAZARDS: dict[str, str] = {
        "glutaraldehyde": "toxic",
        "epichlorohydrin": "toxic_carcinogen",
        "dvs": "toxic",
        "formaldehyde": "toxic_carcinogen",
        "citric_acid": "irritant",
    }

    def _build_safety_warnings(
        self, profile: object, source: str, pH: float,
        reagent_key: str = "",
    ) -> list[str]:
        warnings: list[str] = []

        hazard = getattr(profile, "hazard_class", "")
        if hazard is None:
            hazard = ""

        # CrosslinkerProfile has no hazard_class — use per-key lookup
        if not hazard and source == "crosslinkers":
            hazard = self._CROSSLINKER_HAZARDS.get(reagent_key, "")

        warning_text = _HAZARD_WARNINGS.get(hazard, _HAZARD_WARNINGS[""])
        warnings.append(warning_text)

        if pH > 10.0:
            warnings.append(
                "CAUTION: Strongly alkaline solution. Risk of chemical burns."
            )

        return warnings
