"""Module 2 functionalization orchestrator.

Phase B+: Full 5-workflow Module 2 with backend validation.
Architecture: docs/18_m2_expansion_final_plan.md.

Provides:
- FunctionalMicrosphere: complete description of a functionalized microsphere
- FunctionalMediaContract: stable M2->M3 bridge with ligand density mapping
- ModificationOrchestrator: sequential execution with ACS tracking + validation
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from .acs import ACSSiteType, ACSProfile, initialize_acs_from_m1
from .modification_steps import (
    ModificationResult,
    ModificationStep,
    ModificationStepType,
    solve_modification_step,
)
from .reagent_profiles import REAGENT_PROFILES
from .surface_area import AccessibleSurfaceModel

from typing import Optional

from ..datatypes import ModelEvidenceTier, ModelManifest

if TYPE_CHECKING:
    from emulsim.datatypes import M1ExportContract

logger = logging.getLogger(__name__)


# ─── ProteinPretreatmentState (v5.9.2) ───────────────────────────────

@dataclass
class ProteinPretreatmentState:
    """State of protein after disulfide reduction pretreatment."""
    protein_key: str = ""
    free_thiol_fraction: float = 0.0
    activity_after_reduction: float = 1.0
    reductant_used: str = ""
    excess_reductant_removed: bool = False
    time_since_reduction_s: float = 0.0
    warnings: list[str] = field(default_factory=list)


# ─── FunctionalMicrosphere ────────────────────────────────────────────

@dataclass
class FunctionalMicrosphere:
    """Complete description of a functionalized microsphere.

    Holds the M1 contract, surface model, current ACS state, and the
    full modification history with provenance.

    Attributes:
        m1_contract: Source Module 1 export contract.
        surface_model: Computed accessible surface area model.
        acs_profiles: Current ACS inventory (mutated by modification steps).
        modification_history: Ordered list of completed modification results.
        G_DN_updated: Updated double-network shear modulus [Pa] after
            secondary crosslinking.
        E_star_updated: Updated effective Young's modulus [Pa].
    """
    m1_contract: M1ExportContract
    surface_model: AccessibleSurfaceModel
    acs_profiles: dict[ACSSiteType, ACSProfile]
    modification_history: list[ModificationResult] = field(default_factory=list)
    G_DN_updated: float = 0.0
    E_star_updated: float = 0.0
    # v5.9.2: Protein pretreatment state (set by PROTEIN_PRETREATMENT steps)
    pretreatment_state: ProteinPretreatmentState = field(default_factory=ProteinPretreatmentState)
    # v5.9.4: Residual reagent concentrations after washing
    residual_concentrations: dict = field(default_factory=dict)
    # v6.1 (Node 4): composite manifest summarising all step manifests; the
    # weakest tier across the modification history wins, mirroring how
    # RunReport.compute_min_tier() rolls up M1 levels.
    model_manifest: Optional[ModelManifest] = None

    def validate(self) -> list[str]:
        """Validate ACS conservation across all profiles.

        Returns:
            List of violation messages (empty = all valid).
        """
        errors: list[str] = []
        for profile in self.acs_profiles.values():
            errors.extend(profile.validate())
        return errors


# ─── FunctionalMediaContract (M2→M3 bridge, audit F7) ───────────────

@dataclass
class FunctionalMediaContract:
    """Stable interface between Module 2 and Module 3 for process simulation.

    Maps functional ligand density to estimated chromatographic capacity.
    M3 should consume this contract rather than relying on default isotherm
    parameters when ligand data is available.
    """
    # ── Pass-through from M1 ──
    bead_d50: float = 0.0              # [m]
    porosity: float = 0.0              # [-]
    pore_size_mean: float = 0.0        # [m]

    # ── From M2 ──
    ligand_type: str = "none"          # "iex_anion", "iex_cation", "affinity",
                                       # "imac", "hic", "none"
    installed_ligand: str = ""         # e.g., "DEAE", "Protein A", "Phenyl"
    functional_ligand_density: float = 0.0  # [mol/m^2]
    total_coupled_density: float = 0.0      # [mol/m^2]
    charge_density: float = 0.0        # [mol/m^2] for IEX
    active_protein_density: float = 0.0  # [mol/m^2] for affinity

    # ── Mechanical (pass-through) ──
    G_DN_updated: float = 0.0          # [Pa]
    E_star_updated: float = 0.0        # [Pa]

    # ── Derived M3 parameter estimates ──
    estimated_q_max: float = 0.0       # [mol/m^3 bed] mapped from ligand density
    q_max_confidence: str = "not_mapped"  # "mapped_estimated" | "not_mapped"
    q_max_mapping_notes: str = ""

    # ── Area basis (audit F13) ──
    ligand_density_area_basis: str = ""  # "reagent_accessible", "ligand_accessible", "external"
    q_max_area_basis_note: str = ""      # Human-readable note on q_max derivation

    # ── Binding model hint for M3 (audit F15) ──
    binding_model_hint: str = ""         # Passed from reagent profile for M3 routing

    # ── v5.9.0 FMC v2 fields ──
    reagent_accessible_area_per_bed_volume: float = 0.0  # [m2/m3 bed]
    ligand_accessible_area_per_bed_volume: float = 0.0   # [m2/m3 bed]
    capacity_area_basis: str = ""        # "reagent_accessible" or "ligand_accessible"
    activity_retention_uncertainty: float = 0.0  # +/- on activity_retention
    q_max_lower: float = 0.0            # [mol/m3 bed] lower bound
    q_max_upper: float = 0.0            # [mol/m3 bed] upper bound
    m3_support_level: str = "not_mapped"  # "mapped_quantitative", "mapped_estimated",
                                          # "not_mapped", "requires_user_calibration"
    final_ligand_profile_key: str = ""   # Key of the last coupling reagent profile
    process_state_requirements: str = "" # "salt_concentration", "imidazole", ""
    residual_reagent_warnings: list[str] = field(default_factory=list)

    # ── Trust ──
    confidence_tier: str = "semi_quantitative"
    warnings: list[str] = field(default_factory=list)
    # v6.1 (Node 4): structured evidence record. confidence_tier above is kept
    # as the legacy string field that downstream UI/tests already read; the
    # manifest carries the typed enum, valid_domain, and diagnostics used by
    # RunReport and the trust-aware optimizer (consensus plan Sprint 5).
    model_manifest: Optional[ModelManifest] = None

    def validate_units(self) -> list[str]:
        """Node 10 (F11): boundary-level unit/range sanity checks for M2->M3.

        Same intent as M1ExportContract.validate_units: catch silent unit
        mismatches at the contract boundary. Returns a list of violations
        (empty = pass).
        """
        violations: list[str] = []

        # Geometry inherited from M1.
        if self.bead_d50 != 0.0 and not (1e-7 <= self.bead_d50 <= 1e-2):
            violations.append(
                f"bead_d50={self.bead_d50:g} m outside [1e-7, 1e-2]; wrong unit?"
            )
        if not (0.0 <= self.porosity <= 1.0):
            violations.append(
                f"porosity={self.porosity:g} outside [0, 1]."
            )
        if self.pore_size_mean != 0.0 and not (1e-9 <= self.pore_size_mean <= 1e-5):
            violations.append(
                f"pore_size_mean={self.pore_size_mean:g} m outside [1e-9, 1e-5]."
            )

        # Ligand densities are mol/m^2. Realistic values 1e-7 to 1e-3 mol/m^2.
        for name, val in (
            ("functional_ligand_density", self.functional_ligand_density),
            ("total_coupled_density", self.total_coupled_density),
            ("charge_density", self.charge_density),
            ("active_protein_density", self.active_protein_density),
        ):
            if val < 0.0:
                violations.append(f"{name}={val:g} mol/m^2 must be non-negative.")
            elif val > 1.0:
                violations.append(
                    f"{name}={val:g} mol/m^2 exceeds 1; mol/m^2 vs mol/m^3 confusion?"
                )

        # estimated_q_max is mol/m^3 of bed. Realistic capacities are
        # 0-2000 mol/m^3 for IEX, 0-100 for affinity. Block obviously wrong.
        if self.estimated_q_max < 0.0:
            violations.append(
                f"estimated_q_max={self.estimated_q_max:g} mol/m^3 must be non-negative."
            )
        if self.estimated_q_max > 1e5:
            violations.append(
                f"estimated_q_max={self.estimated_q_max:g} mol/m^3 exceeds 1e5; "
                "wrong unit?"
            )
        if self.q_max_lower > self.q_max_upper and self.q_max_upper > 0.0:
            violations.append(
                f"q_max_lower={self.q_max_lower:g} > q_max_upper={self.q_max_upper:g}; "
                "uncertainty bounds inverted."
            )

        # Mechanical inherited from M1 (Pa).
        if self.G_DN_updated != 0.0 and not (1.0 <= self.G_DN_updated <= 1e9):
            violations.append(
                f"G_DN_updated={self.G_DN_updated:g} Pa outside [1, 1e9]."
            )
        if self.E_star_updated != 0.0 and not (1.0 <= self.E_star_updated <= 1e10):
            violations.append(
                f"E_star_updated={self.E_star_updated:g} Pa outside [1, 1e10]."
            )

        if not (0.0 <= self.activity_retention_uncertainty <= 1.0):
            violations.append(
                f"activity_retention_uncertainty={self.activity_retention_uncertainty:g} "
                "outside [0, 1]."
            )

        return violations


def build_functional_media_contract(
    microsphere: FunctionalMicrosphere,
) -> FunctionalMediaContract:
    """Build M2→M3 bridge contract from a functionalized microsphere.

    Scans modification history to identify the installed ligand type and
    computes functional ligand density for M3 capacity estimation.

    Args:
        microsphere: Completed FunctionalMicrosphere from orchestrator.

    Returns:
        FunctionalMediaContract with ligand mapping.
    """
    contract = microsphere.m1_contract
    surface = microsphere.surface_model
    warnings: list[str] = []

    # Find the last coupling step to determine ligand type
    ligand_type = "none"
    installed_ligand = ""
    _last_coupling_rp = None
    functional_density = 0.0
    coupled_density = 0.0
    confidence = "not_mapped"
    q_max_est = 0.0
    q_max_notes = ""

    for mr in microsphere.modification_history:
        step = mr.step
        if step.step_type in (ModificationStepType.LIGAND_COUPLING,
                              ModificationStepType.PROTEIN_COUPLING):
            rp = REAGENT_PROFILES.get(step.reagent_key)
            if rp is not None:
                installed_ligand = getattr(rp, 'installed_ligand', step.reagent_key)
                fm = getattr(rp, 'functional_mode', '')

                # Map functional_mode to ligand_type
                # Audit F14: use charge_type for IEX instead of string heuristic
                _ct = getattr(rp, 'charge_type', '')
                _mode_map = {
                    "iex_ligand": "iex_anion" if _ct == "anion" else (
                        "iex_cation" if _ct == "cation" else "iex_anion"),
                    "affinity_ligand": "affinity",
                    "imac_chelator": "imac",
                    "hic_ligand": "hic",
                    "gst_affinity": "gst_affinity",
                    "biotin_affinity": "biotin_affinity",
                    "heparin_affinity": "heparin_affinity",
                }
                ligand_type = _mode_map.get(fm, "none")
                _last_coupling_rp = rp  # carry for density area selection

    # Compute densities from ACS profiles (audit F13: track area basis)
    _is_macro = getattr(_last_coupling_rp, 'is_macromolecule', False) if _last_coupling_rp is not None else False
    _area_basis = "reagent_accessible"
    for _st, profile in microsphere.acs_profiles.items():
        if profile.ligand_coupled_sites > 0:
            if _is_macro and surface.ligand_accessible_area > 0:
                area = surface.ligand_accessible_area
                _area_basis = "ligand_accessible"
            else:
                area = surface.reagent_accessible_area if surface.reagent_accessible_area > 0 else surface.ligand_accessible_area
                _area_basis = "reagent_accessible" if surface.reagent_accessible_area > 0 else "ligand_accessible"
            if area > 0:
                coupled_density = profile.ligand_coupled_sites / area
                functional_density = profile.ligand_functional_sites / area

    # ── v5.9.0: Compute accessible area per bed volume ──
    d_p = contract.bead_d50
    eps_b = 0.38  # default bed porosity
    _binding_hint = getattr(_last_coupling_rp, 'binding_model_hint', '') if _last_coupling_rp else ''
    _m3_support = getattr(_last_coupling_rp, 'm3_support_level', 'not_mapped') if _last_coupling_rp else 'not_mapped'
    _final_key = getattr(_last_coupling_rp, 'name', '') if _last_coupling_rp else ''
    _act_unc = getattr(_last_coupling_rp, 'activity_retention_uncertainty', 0.0) if _last_coupling_rp else 0.0

    # Accessible area per bed volume (v5.9.0 WN-0b)
    _reagent_a_v = 0.0
    _ligand_a_v = 0.0
    if d_p > 0:
        import math as _math_fmc
        V_particle = (4.0 / 3.0) * _math_fmc.pi * (d_p / 2.0) ** 3
        if V_particle > 0:
            particles_per_bed_vol = (1.0 - eps_b) / V_particle
            _reagent_a_v = surface.reagent_accessible_area * particles_per_bed_vol
            _ligand_a_v = surface.ligand_accessible_area * particles_per_bed_vol

    # Choose appropriate a_v for q_max based on molecule size
    if _is_macro:
        a_v_for_qmax = _ligand_a_v
        _cap_basis = "ligand_accessible"
    else:
        a_v_for_qmax = _reagent_a_v if _reagent_a_v > 0 else _ligand_a_v
        _cap_basis = "reagent_accessible" if _reagent_a_v > 0 else "ligand_accessible"

    _q_max_area_note = f"Area basis: {_area_basis}. Bed a_v({_cap_basis}): {a_v_for_qmax:.0f} m2/m3."

    if d_p > 0 and functional_density > 0 and a_v_for_qmax > 0:

        if ligand_type in ("iex_anion", "iex_cation"):
            q_max_est = functional_density * a_v_for_qmax
            confidence = "mapped_estimated"
            q_max_notes = (
                f"IEX: q_max = ligand_density * a_v = {q_max_est:.2f} mol/m^3. "
                f"a_v = {a_v_for_qmax:.0f} m^2/m^3 ({_cap_basis}). {_q_max_area_note} "
                f"Bed porosity = {eps_b:.2f}."
            )
        elif ligand_type == "affinity":
            binding_stoich = 2.0  # Protein A/G: ~2 IgG per ligand
            q_max_est = functional_density * a_v_for_qmax * binding_stoich
            confidence = "mapped_estimated"
            q_max_notes = (
                f"Affinity (Fc): q_max = density * a_v * stoich({binding_stoich:.0f}). "
                f"{_q_max_area_note} Ranking only."
            )
        elif ligand_type == "imac":
            binding_stoich = 1.0
            _mi = getattr(_last_coupling_rp, 'metal_ion', 'Ni2+') if _last_coupling_rp else 'Ni2+'
            _mlf = getattr(_last_coupling_rp, 'metal_loaded_fraction', 1.0) if _last_coupling_rp else 1.0
            q_max_est = functional_density * a_v_for_qmax * binding_stoich * _mlf
            confidence = "mapped_estimated"
            q_max_notes = (
                f"IMAC: q_max = density * a_v * stoich * metal_frac. "
                f"Assumes fully {_mi}-loaded (frac={_mlf:.0%}), "
                f"no leaching, no competing chelators. {_q_max_area_note}"
            )
        elif ligand_type == "gst_affinity":
            binding_stoich = 1.0  # 1 GST per glutathione
            q_max_est = functional_density * a_v_for_qmax * binding_stoich
            confidence = "mapped_estimated"
            q_max_notes = (
                f"GST affinity: 1:1 GST:glutathione. {_q_max_area_note}"
            )
        elif ligand_type == "biotin_affinity":
            binding_stoich = 2.5  # Audit F7: capped from theoretical 4
            q_max_est = functional_density * a_v_for_qmax * binding_stoich
            confidence = "mapped_estimated"
            q_max_notes = (
                f"Biotin affinity: stoich={binding_stoich} (capped from theoretical 4 "
                f"due to steric occlusion). Near-irreversible binding. {_q_max_area_note}"
            )
        elif ligand_type == "heparin_affinity":
            binding_stoich = 1.0  # approximate, highly target-dependent
            q_max_est = functional_density * a_v_for_qmax * binding_stoich
            confidence = "mapped_estimated"
            q_max_notes = (
                f"Heparin affinity: stoich~1 (approximate, target-dependent). "
                f"Mixed affinity + cation exchange. {_q_max_area_note}"
            )
        elif ligand_type == "hic":
            q_max_notes = (
                f"HIC: q_max not mappable from ligand density alone. "
                f"Requires salt-dependent adsorption isotherm. {_q_max_area_note}"
            )
        else:
            q_max_notes = f"Ligand type '{ligand_type}' — q_max mapping not implemented."
            if ligand_type != "none":
                warnings.append(q_max_notes)
    else:
        if ligand_type != "none" and functional_density > 0:
            q_max_notes = "q_max not computed: bead diameter is zero."
        elif ligand_type != "none":
            q_max_notes = "q_max not computed: no functional ligand density."

    # Determine confidence tier
    _ranking_types = {"affinity", "biotin_affinity", "heparin_affinity"}
    _conf_tier = "ranking_only" if ligand_type in _ranking_types else "semi_quantitative"

    # Compute q_max uncertainty bounds from activity_retention_uncertainty
    _q_lower = 0.0
    _q_upper = 0.0
    if q_max_est > 0 and _act_unc > 0:
        # q_max scales linearly with activity_retention
        _act_ret = getattr(_last_coupling_rp, 'activity_retention', 1.0) if _last_coupling_rp else 1.0
        if _act_ret > 0:
            _q_lower = q_max_est * max(_act_ret - _act_unc, 0.0) / _act_ret
            _q_upper = q_max_est * min(_act_ret + _act_unc, 1.0) / _act_ret

    # Process state requirements for M3 routing
    _proc_req = ""
    if _binding_hint == "charge_exchange":
        _proc_req = "salt_concentration"
    elif _binding_hint == "metal_chelation":
        _proc_req = "imidazole"

    # Build FMC manifest by combining the microsphere composite with the
    # FMC-level mapping evidence (e.g., ranking_only ligand types, m3_support).
    fmc_manifest = _build_fmc_manifest(
        microsphere_manifest=microsphere.model_manifest,
        confidence_tier_str=_conf_tier,
        ligand_type=ligand_type,
        m3_support_level=_m3_support,
        q_max_est=q_max_est,
        q_max_notes=q_max_notes,
        warnings=warnings,
    )

    fmc = FunctionalMediaContract(
        bead_d50=contract.bead_d50,
        porosity=contract.porosity,
        pore_size_mean=contract.pore_size_mean,
        ligand_type=ligand_type,
        installed_ligand=installed_ligand,
        functional_ligand_density=functional_density,
        total_coupled_density=coupled_density,
        charge_density=functional_density if ligand_type.startswith("iex") else 0.0,
        active_protein_density=functional_density if ligand_type in ("affinity", "biotin_affinity") else 0.0,
        G_DN_updated=microsphere.G_DN_updated,
        E_star_updated=microsphere.E_star_updated,
        estimated_q_max=q_max_est,
        q_max_confidence=confidence,
        q_max_mapping_notes=q_max_notes,
        ligand_density_area_basis=_area_basis,
        q_max_area_basis_note=_q_max_area_note,
        binding_model_hint=_binding_hint,
        # v5.9.0 FMC v2 fields
        reagent_accessible_area_per_bed_volume=_reagent_a_v,
        ligand_accessible_area_per_bed_volume=_ligand_a_v,
        capacity_area_basis=_cap_basis,
        activity_retention_uncertainty=_act_unc,
        q_max_lower=_q_lower,
        q_max_upper=_q_upper,
        m3_support_level=_m3_support,
        final_ligand_profile_key=_final_key,
        process_state_requirements=_proc_req,
        confidence_tier=_conf_tier,
        warnings=warnings,
        model_manifest=fmc_manifest,
    )

    # Node 10 (F11): boundary unit/range check at the M2->M3 contract.
    _unit_violations = fmc.validate_units()
    if _unit_violations:
        logger.warning(
            "FunctionalMediaContract failed %d unit/range check(s):\n  %s",
            len(_unit_violations),
            "\n  ".join(_unit_violations),
        )
        fmc.warnings.extend(
            f"FMC unit check: {v}" for v in _unit_violations
        )

    return fmc


def _build_fmc_manifest(
    microsphere_manifest: Optional[ModelManifest],
    confidence_tier_str: str,
    ligand_type: str,
    m3_support_level: str,
    q_max_est: float,
    q_max_notes: str,
    warnings: list[str],
) -> ModelManifest:
    """Build the FMC manifest combining microsphere evidence with FMC mapping.

    The FMC tier is the weaker of the upstream microsphere tier and the FMC's
    own q_max-mapping confidence. This is the M2->M3 evidence handoff:
    M3 should never claim better evidence than the FMC's worst input.

    Tier mapping for the legacy `confidence_tier` string:
      ranking_only          -> QUALITATIVE_TREND
      semi_quantitative     -> SEMI_QUANTITATIVE
      mapped_estimated      -> SEMI_QUANTITATIVE (q_max from accessibility model,
                              not from a measured isotherm)
      requires_user_calibration -> SEMI_QUANTITATIVE with diagnostic flag
    """
    _STR_TO_TIER = {
        "ranking_only": ModelEvidenceTier.QUALITATIVE_TREND,
        "semi_quantitative": ModelEvidenceTier.SEMI_QUANTITATIVE,
        "mapped_estimated": ModelEvidenceTier.SEMI_QUANTITATIVE,
        "requires_user_calibration": ModelEvidenceTier.SEMI_QUANTITATIVE,
        "not_mapped": ModelEvidenceTier.UNSUPPORTED,
        "validated": ModelEvidenceTier.VALIDATED_QUANTITATIVE,
    }
    fmc_own_tier = _STR_TO_TIER.get(confidence_tier_str, ModelEvidenceTier.SEMI_QUANTITATIVE)

    # If no upstream manifest, we can only attest to FMC-level evidence.
    if microsphere_manifest is None:
        composite_tier = fmc_own_tier
    else:
        # Weakest wins (largest index in the tier list = weakest)
        _ORDER = list(ModelEvidenceTier)
        composite_idx = max(
            _ORDER.index(microsphere_manifest.evidence_tier),
            _ORDER.index(fmc_own_tier),
        )
        composite_tier = _ORDER[composite_idx]

    return ModelManifest(
        model_name=f"M2.FMC.{ligand_type or 'none'}",
        evidence_tier=composite_tier,
        diagnostics={
            "fmc_confidence_tier": confidence_tier_str,
            "m3_support_level": m3_support_level,
            "estimated_q_max": float(q_max_est),
            "ligand_type": ligand_type,
            "n_warnings": len(warnings),
        },
        assumptions=[q_max_notes] if q_max_notes else [],
    )


# ─── ModificationOrchestrator ─────────────────────────────────────────

class ModificationOrchestrator:
    """Execute sequential modification steps on a microsphere.

    Usage:
        orchestrator = ModificationOrchestrator()
        result = orchestrator.run(contract, steps)
        assert result.validate() == []
    """

    def run(
        self,
        contract: M1ExportContract,
        steps: list[ModificationStep],
    ) -> FunctionalMicrosphere:
        """Execute all modification steps sequentially, tracking ACS.

        Algorithm:
            1. Build AccessibleSurfaceModel from M1 contract.
            2. Initialize ACS inventory from contract + surface model.
            3. For each step: look up reagent, solve, record result.
            4. Accumulate G_DN updates from crosslinking steps.
            5. Return FunctionalMicrosphere with full provenance.

        Args:
            contract: M1ExportContract from Module 1.
            steps: Ordered list of modification steps to execute.

        Returns:
            FunctionalMicrosphere with updated ACS and modification history.

        Raises:
            KeyError: If a step references an unknown reagent_key.
            ValueError: If a step targets a missing ACS type.
        """
        # --- Build surface model ---
        surface_model = AccessibleSurfaceModel.from_m1_export(contract)

        # --- Initialize ACS from M1 ---
        acs_profiles = initialize_acs_from_m1(contract, surface_model)

        logger.info(
            "ModificationOrchestrator: %d steps, %d ACS types initialized.",
            len(steps), len(acs_profiles),
        )

        # --- Validate workflow ordering (audit F8: backend, not UI-only) ---
        _validate_workflow_ordering(steps, acs_profiles)

        # --- Execute steps sequentially ---
        history: list[ModificationResult] = []
        total_delta_G = 0.0

        for i, step in enumerate(steps):
            # Look up reagent profile
            if step.reagent_key not in REAGENT_PROFILES:
                raise KeyError(
                    f"Step {i}: unknown reagent_key '{step.reagent_key}'. "
                    f"Available: {list(REAGENT_PROFILES.keys())}"
                )
            reagent_profile = REAGENT_PROFILES[step.reagent_key]

            logger.info(
                "Step %d/%d: %s with %s on %s",
                i + 1, len(steps), step.step_type.value,
                step.reagent_key, step.target_acs.value,
            )

            result = solve_modification_step(
                step=step,
                acs_state=acs_profiles,
                surface_model=surface_model,
                reagent_profile=reagent_profile,
            )
            history.append(result)
            total_delta_G += result.delta_G_DN

            logger.info(
                "Step %d: conversion=%.4f, delta_G=%.1f Pa",
                i + 1, result.conversion, result.delta_G_DN,
            )

        # --- Compute updated mechanical properties ---
        G_DN_base = contract.G_DN
        G_DN_updated = G_DN_base + total_delta_G
        # E* ~ 3*G for incompressible rubber (Poisson's ratio ~ 0.5)
        E_star_updated = 3.0 * G_DN_updated

        logger.info(
            "Orchestrator complete: G_DN %.1f -> %.1f Pa (delta=%.1f Pa)",
            G_DN_base, G_DN_updated, total_delta_G,
        )

        # Composite manifest for the microsphere — weakest step tier wins.
        composite_manifest = _build_microsphere_manifest(history)

        return FunctionalMicrosphere(
            m1_contract=contract,
            surface_model=surface_model,
            acs_profiles=acs_profiles,
            modification_history=history,
            G_DN_updated=G_DN_updated,
            E_star_updated=E_star_updated,
            model_manifest=composite_manifest,
        )


def _build_microsphere_manifest(
    history: list[ModificationResult],
) -> ModelManifest:
    """Aggregate per-step manifests into a single FunctionalMicrosphere manifest.

    Tier rule: the weakest tier across the history wins (UNSUPPORTED >
    QUALITATIVE_TREND > SEMI_QUANTITATIVE > CALIBRATED_LOCAL >
    VALIDATED_QUANTITATIVE). Diagnostics aggregate per-step counts so a
    consumer can tell whether any step was a fallback.

    Empty history returns an UNSUPPORTED manifest — the orchestrator was
    invoked with no chemistry, so no functional surface evidence exists.
    """
    if not history:
        return ModelManifest(
            model_name="M2.composite",
            evidence_tier=ModelEvidenceTier.UNSUPPORTED,
            diagnostics={"n_steps": 0},
            assumptions=["No modification steps executed."],
        )

    _ORDER = list(ModelEvidenceTier)  # validated ... unsupported, in order
    worst_idx = 0
    step_summaries: list[dict] = []
    assumptions: list[str] = []
    for i, mr in enumerate(history):
        m = mr.model_manifest
        if m is None:
            # Defensive: a step without a manifest is treated as the worst tier
            # so the composite cannot be silently upgraded by missing data.
            worst_idx = max(worst_idx, _ORDER.index(ModelEvidenceTier.UNSUPPORTED))
            step_summaries.append({"step": i + 1, "missing_manifest": True})
            continue
        worst_idx = max(worst_idx, _ORDER.index(m.evidence_tier))
        step_summaries.append({
            "step": i + 1,
            "name": m.model_name,
            "tier": m.evidence_tier.value,
            "conversion": m.diagnostics.get("conversion"),
        })
        assumptions.extend(m.assumptions)

    return ModelManifest(
        model_name="M2.composite",
        evidence_tier=_ORDER[worst_idx],
        diagnostics={
            "n_steps": len(history),
            "steps": step_summaries,
            "weakest_tier": _ORDER[worst_idx].value,
        },
        # De-dupe assumptions while preserving order
        assumptions=list(dict.fromkeys(assumptions)),
    )


# ─── Backend Workflow Validation (audit F8) ──────────────────────────

_COUPLING_TYPES = {
    ModificationStepType.LIGAND_COUPLING,
    ModificationStepType.PROTEIN_COUPLING,
}

# Types that require activated sites on target (includes SPACER_ARM)
_REQUIRES_ACTIVATED = _COUPLING_TYPES | {
    ModificationStepType.QUENCHING,
    ModificationStepType.SPACER_ARM,
}

# Rule 4: Allowed reaction_type values per step type (Codex P1-1 fix)
_STEP_ALLOWED_REACTION_TYPES: dict[ModificationStepType, set[str]] = {
    ModificationStepType.SECONDARY_CROSSLINKING: {"crosslinking"},
    ModificationStepType.ACTIVATION: {"activation"},
    ModificationStepType.LIGAND_COUPLING: {"coupling"},
    ModificationStepType.PROTEIN_COUPLING: {"protein_coupling"},
    ModificationStepType.QUENCHING: {"blocking"},
    ModificationStepType.SPACER_ARM: {"spacer_arm", "heterobifunctional"},
    ModificationStepType.METAL_CHARGING: {"metal_charging", "metal_stripping"},
    ModificationStepType.PROTEIN_PRETREATMENT: {"protein_pretreatment"},
    ModificationStepType.WASHING: {"washing"},
}


def _validate_workflow_ordering(
    steps: list[ModificationStep],
    acs_profiles: dict[ACSSiteType, ACSProfile],
) -> None:
    """Validate step ordering before execution (backend authority).

    Rules:
        1. Coupling/quenching requires activated sites on target ACS type.
        2. No steps after quenching on the same target ACS type.
        3. Reagent-target ACS type must match reagent profile's target_acs.
        4. Reagent reaction_type must be compatible with step type.

    Raises:
        ValueError: On blocker-level violations.

    Logs warnings for advisory issues.
    """
    quenched_targets: set[ACSSiteType] = set()

    for i, step in enumerate(steps):
        idx = i + 1

        # Rule 2: No steps after quenching on same target.
        # NOTE: Cross-target workflows (e.g., Quench(EPOXIDE) then Crosslink(AMINE))
        # are intentionally allowed — quenching one group does not block chemistry
        # on a chemically distinct group (validated against wetlab practice).
        if step.target_acs in quenched_targets:
            raise ValueError(
                f"Step {idx}: {step.step_type.value} targets "
                f"{step.target_acs.value} which was already quenched. "
                f"No further chemistry is possible on blocked sites."
            )

        # Rule 1: Coupling/quenching/spacer_arm requires activated sites
        if step.step_type in _REQUIRES_ACTIVATED:
            target_profile = acs_profiles.get(step.target_acs)
            if target_profile is None:
                # Target ACS type doesn't exist yet — may be created by prior step
                prior_produces = any(
                    s.product_acs == step.target_acs
                    for s in steps[:i]
                    if s.step_type in (ModificationStepType.ACTIVATION,
                                       ModificationStepType.SPACER_ARM)
                )
                if not prior_produces:
                    raise ValueError(
                        f"Step {idx}: {step.step_type.value} targets "
                        f"{step.target_acs.value} but no prior activation step "
                        f"produces this site type and it doesn't exist in M1 ACS."
                    )
            elif target_profile.activated_sites <= 0 and target_profile.remaining_activated <= 0:
                # Check if a prior step in this batch will activate it
                prior_activates = any(
                    s.product_acs == step.target_acs
                    for s in steps[:i]
                    if s.step_type in (ModificationStepType.ACTIVATION,
                                       ModificationStepType.SPACER_ARM)
                )
                if not prior_activates:
                    logger.warning(
                        "Step %d: %s targets %s with 0 activated sites. "
                        "Coupling/quenching will have zero conversion.",
                        idx, step.step_type.value, step.target_acs.value,
                    )

        # Rule 3: Reagent profile compatibility
        if step.reagent_key in REAGENT_PROFILES:
            rp = REAGENT_PROFILES[step.reagent_key]
            if rp.target_acs != step.target_acs:
                raise ValueError(
                    f"Step {idx}: reagent '{step.reagent_key}' targets "
                    f"{rp.target_acs.value} but step targets {step.target_acs.value}. "
                    f"Reagent-target mismatch."
                )

        # Rule 4: Reagent reaction_type must match step type (Codex P1-1 fix)
        if step.reagent_key in REAGENT_PROFILES:
            rp_r4 = REAGENT_PROFILES[step.reagent_key]
            allowed = _STEP_ALLOWED_REACTION_TYPES.get(step.step_type, set())
            if allowed and rp_r4.reaction_type not in allowed:
                raise ValueError(
                    f"Step {idx}: reagent '{step.reagent_key}' has reaction_type "
                    f"'{rp_r4.reaction_type}' which is incompatible with step type "
                    f"'{step.step_type.value}'. Allowed reaction_types: {sorted(allowed)}."
                )

        # Track quenching
        if step.step_type == ModificationStepType.QUENCHING:
            quenched_targets.add(step.target_acs)
