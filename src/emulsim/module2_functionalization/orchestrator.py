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
from .reagent_profiles import REAGENT_PROFILES, ReagentProfile
from .surface_area import AccessibleSurfaceModel

if TYPE_CHECKING:
    from emulsim.datatypes import M1ExportContract

logger = logging.getLogger(__name__)


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

    # ── Trust ──
    confidence_tier: str = "semi_quantitative"
    warnings: list[str] = field(default_factory=list)


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
                _mode_map = {
                    "iex_ligand": "iex_anion" if "anion" in installed_ligand.lower() or "deae" in installed_ligand.lower()
                                  else "iex_cation",
                    "affinity_ligand": "affinity",
                    "imac_chelator": "imac",
                    "hic_ligand": "hic",
                }
                ligand_type = _mode_map.get(fm, "none")

    # Compute densities from ACS profiles
    for _st, profile in microsphere.acs_profiles.items():
        if profile.ligand_coupled_sites > 0:
            area = surface.ligand_accessible_area if surface.ligand_accessible_area > 0 else surface.reagent_accessible_area
            if area > 0:
                coupled_density = profile.ligand_coupled_sites / area
                functional_density = profile.ligand_functional_sites / area

    # Estimate q_max for IEX and affinity
    if ligand_type in ("iex_anion", "iex_cation") and functional_density > 0:
        # Ionic capacity = functional_density * specific_surface_area * (1 - bed_porosity)
        # Approximate specific surface area: 6*(1-eps)/(d_p)
        d_p = contract.bead_d50
        eps_b = 0.38  # default bed porosity
        if d_p > 0:
            a_v = 6.0 * (1.0 - eps_b) / d_p  # [m^2/m^3 bed]
            q_max_est = functional_density * a_v  # [mol/m^3 bed]
            confidence = "mapped_estimated"
            q_max_notes = (
                f"IEX: q_max = ligand_density * a_v. "
                f"a_v = 6(1-eps)/d_p = {a_v:.0f} m^2/m^3. "
                f"Assumes bed porosity = {eps_b:.2f}."
            )
    elif ligand_type == "affinity" and functional_density > 0:
        d_p = contract.bead_d50
        eps_b = 0.38
        if d_p > 0:
            a_v = 6.0 * (1.0 - eps_b) / d_p
            binding_stoich = 2.0  # 1 Protein A binds 2 IgG Fc
            q_max_est = functional_density * a_v * binding_stoich
            confidence = "mapped_estimated"
            q_max_notes = (
                f"Affinity: q_max = functional_density * a_v * stoich. "
                f"Binding stoichiometry = {binding_stoich:.0f} (Protein A:IgG). "
                f"Ranking only — calibrate with frontal analysis."
            )
    else:
        q_max_notes = f"Ligand type '{ligand_type}' — q_max mapping not implemented."
        if ligand_type != "none":
            warnings.append(q_max_notes)

    return FunctionalMediaContract(
        bead_d50=contract.bead_d50,
        porosity=contract.porosity,
        pore_size_mean=contract.pore_size_mean,
        ligand_type=ligand_type,
        installed_ligand=installed_ligand,
        functional_ligand_density=functional_density,
        total_coupled_density=coupled_density,
        charge_density=functional_density if ligand_type.startswith("iex") else 0.0,
        active_protein_density=functional_density if ligand_type == "affinity" else 0.0,
        G_DN_updated=microsphere.G_DN_updated,
        E_star_updated=microsphere.E_star_updated,
        estimated_q_max=q_max_est,
        q_max_confidence=confidence,
        q_max_mapping_notes=q_max_notes,
        confidence_tier="ranking_only" if ligand_type == "affinity" else "semi_quantitative",
        warnings=warnings,
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

        return FunctionalMicrosphere(
            m1_contract=contract,
            surface_model=surface_model,
            acs_profiles=acs_profiles,
            modification_history=history,
            G_DN_updated=G_DN_updated,
            E_star_updated=E_star_updated,
        )


# ─── Backend Workflow Validation (audit F8) ──────────────────────────

_COUPLING_TYPES = {
    ModificationStepType.LIGAND_COUPLING,
    ModificationStepType.PROTEIN_COUPLING,
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

    Raises:
        ValueError: On blocker-level violations.

    Logs warnings for advisory issues.
    """
    quenched_targets: set[ACSSiteType] = set()

    for i, step in enumerate(steps):
        idx = i + 1

        # Rule 2: No steps after quenching on same target
        if step.target_acs in quenched_targets:
            raise ValueError(
                f"Step {idx}: {step.step_type.value} targets "
                f"{step.target_acs.value} which was already quenched. "
                f"No further chemistry is possible on blocked sites."
            )

        # Rule 1: Coupling/quenching requires activated sites
        if step.step_type in _COUPLING_TYPES or step.step_type == ModificationStepType.QUENCHING:
            target_profile = acs_profiles.get(step.target_acs)
            if target_profile is None:
                # Target ACS type doesn't exist yet — it may be created by
                # a prior activation step. Check if any prior step produces it.
                prior_produces = any(
                    s.product_acs == step.target_acs
                    for s in steps[:i]
                    if s.step_type == ModificationStepType.ACTIVATION
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
                    if s.step_type == ModificationStepType.ACTIVATION
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

        # Track quenching
        if step.step_type == ModificationStepType.QUENCHING:
            quenched_targets.add(step.target_acs)
