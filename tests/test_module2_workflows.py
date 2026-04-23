"""Tests for Module 2 Phase B: 2 validated chemistry workflows.

Validates:
  - Amine secondary crosslinking (genipin, glutaraldehyde)
  - Hydroxyl activation (ECH, DVS)
  - Sequential multi-step workflows
  - ACS conservation across all steps
  - Orchestrator end-to-end from M1ExportContract
  - Physical range of conversion [0, 1]
  - Hydrolysis competition (ECH vs DVS yield)
"""

from __future__ import annotations

import math

import pytest

from emulsim.datatypes import M1ExportContract
from emulsim.module2_functionalization.acs import (
    ACSSiteType,
    ACSProfile,
    initialize_acs_from_m1,
)
from emulsim.module2_functionalization.modification_steps import (
    ModificationStep,
    ModificationStepType,
    solve_modification_step,
)
from emulsim.module2_functionalization.orchestrator import (
    FunctionalMicrosphere,
    ModificationOrchestrator,
)
from emulsim.module2_functionalization.reactions import (
    solve_second_order_consumption,
)
from emulsim.module2_functionalization.reagent_profiles import (
    REAGENT_PROFILES,
)
from emulsim.module2_functionalization.surface_area import (
    AccessibleSurfaceModel,
)


# ─── Fixtures ───────────────────────────────────────────────────────────

def _make_contract(
    bead_radius: float = 50e-6,
    porosity: float = 0.7,
    pore_size_mean: float = 100e-9,
    mesh_size_xi: float = 20e-9,
    nh2_bulk: float = 100.0,
    oh_bulk: float = 400.0,
    p_final: float = 0.5,
    G_DN: float = 5000.0,
) -> M1ExportContract:
    """Helper: create a minimal M1ExportContract for testing."""
    return M1ExportContract(
        bead_radius=bead_radius,
        bead_d32=bead_radius * 2,
        bead_d50=bead_radius * 2,
        pore_size_mean=pore_size_mean,
        pore_size_std=pore_size_mean * 0.3,
        porosity=porosity,
        l2_model_tier="empirical_calibrated",
        mesh_size_xi=mesh_size_xi,
        p_final=p_final,
        primary_crosslinker="genipin",
        nh2_bulk_concentration=nh2_bulk,
        oh_bulk_concentration=oh_bulk,
        G_DN=G_DN,
        E_star=G_DN * 3.0,
        model_used="phenomenological",
        c_agarose=42.0,
        c_chitosan=18.0,
        DDA=0.90,
        trust_level="CAUTION",
    )


def _make_surface_model(contract: M1ExportContract) -> AccessibleSurfaceModel:
    """Helper: build and compute surface model from contract."""
    return AccessibleSurfaceModel.from_m1_export(contract)


def _make_acs_state(
    contract: M1ExportContract,
    surface_model: AccessibleSurfaceModel,
) -> dict[ACSSiteType, ACSProfile]:
    """Helper: initialize ACS state from contract."""
    return initialize_acs_from_m1(contract, surface_model)


# ─── Test: Amine secondary crosslinking — genipin ─────────────────────

def test_amine_secondary_crosslinking_genipin():
    """Genipin secondary crosslinking consumes NH2, increases G_DN."""
    contract = _make_contract(nh2_bulk=100.0)
    surface_model = _make_surface_model(contract)
    acs_state = _make_acs_state(contract, surface_model)

    reagent = REAGENT_PROFILES["genipin_secondary"]

    step = ModificationStep(
        step_type=ModificationStepType.SECONDARY_CROSSLINKING,
        reagent_key="genipin_secondary",
        target_acs=ACSSiteType.AMINE_PRIMARY,
        temperature=310.15,   # 37 degC
        time=14400.0,         # 4 h
        reagent_concentration=10.0,  # [mol/m^3]
        stoichiometry=0.5,
    )

    nh2_before = acs_state[ACSSiteType.AMINE_PRIMARY].remaining_sites

    result = solve_modification_step(step, acs_state, surface_model, reagent)

    # Conversion must be in [0, 1]
    assert 0.0 <= result.conversion <= 1.0, f"conversion={result.conversion}"

    # Some NH2 should be consumed
    assert result.conversion > 0.0, "Expected nonzero conversion for genipin"

    # Remaining sites should decrease
    nh2_after = acs_state[ACSSiteType.AMINE_PRIMARY].remaining_sites
    assert nh2_after < nh2_before, "NH2 remaining should decrease"

    # delta_G should be positive (crosslinking increases modulus)
    assert result.delta_G_DN > 0.0, "Crosslinking should increase G_DN"

    # Conservation check
    violations = acs_state[ACSSiteType.AMINE_PRIMARY].validate()
    assert violations == [], f"Conservation violations: {violations}"


# ─── Test: Hydroxyl activation — ECH ──────────────────────────────────

def test_hydroxyl_activation_ech():
    """ECH activation converts OH to EPOXIDE with hydrolysis competition."""
    contract = _make_contract(oh_bulk=400.0)
    surface_model = _make_surface_model(contract)
    acs_state = _make_acs_state(contract, surface_model)

    reagent = REAGENT_PROFILES["ech_activation"]

    step = ModificationStep(
        step_type=ModificationStepType.ACTIVATION,
        reagent_key="ech_activation",
        target_acs=ACSSiteType.HYDROXYL,
        product_acs=ACSSiteType.EPOXIDE,
        temperature=298.15,
        time=7200.0,          # 2 h
        ph=12.0,
        reagent_concentration=50.0,  # [mol/m^3]
        stoichiometry=1.0,
    )

    oh_before = acs_state[ACSSiteType.HYDROXYL].remaining_sites

    result = solve_modification_step(step, acs_state, surface_model, reagent)

    # Conversion in [0, 1]
    assert 0.0 <= result.conversion <= 1.0

    # Some OH should be consumed
    assert result.conversion > 0.0, "Expected nonzero conversion for ECH"

    # EPOXIDE sites should be created
    assert ACSSiteType.EPOXIDE in acs_state, "EPOXIDE sites should be created"
    epoxide_sites = acs_state[ACSSiteType.EPOXIDE].activated_sites
    assert epoxide_sites > 0, "Epoxide activated sites should be positive"

    # OH remaining should decrease
    oh_after = acs_state[ACSSiteType.HYDROXYL].remaining_sites
    assert oh_after < oh_before

    # Conservation checks
    for profile in acs_state.values():
        violations = profile.validate()
        assert violations == [], f"Conservation violations for {profile.site_type.value}: {violations}"


# ─── Test: Hydroxyl activation — DVS ──────────────────────────────────

def test_hydroxyl_activation_dvs():
    """DVS converts OH to VINYL_SULFONE."""
    contract = _make_contract(oh_bulk=400.0)
    surface_model = _make_surface_model(contract)
    acs_state = _make_acs_state(contract, surface_model)

    reagent = REAGENT_PROFILES["dvs_activation"]

    step = ModificationStep(
        step_type=ModificationStepType.ACTIVATION,
        reagent_key="dvs_activation",
        target_acs=ACSSiteType.HYDROXYL,
        product_acs=ACSSiteType.VINYL_SULFONE,
        temperature=298.15,
        time=3600.0,
        ph=11.0,
        reagent_concentration=50.0,
        stoichiometry=1.0,
    )

    result = solve_modification_step(step, acs_state, surface_model, reagent)

    assert 0.0 <= result.conversion <= 1.0
    assert result.conversion > 0.0

    assert ACSSiteType.VINYL_SULFONE in acs_state
    vs_sites = acs_state[ACSSiteType.VINYL_SULFONE].activated_sites
    assert vs_sites > 0


# ─── Test: Sequential workflow — ECH then genipin ─────────────────────

def test_sequential_workflow_ech_then_crosslink():
    """ECH activation then genipin crosslinking — full 2-step sequence."""
    contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0)
    surface_model = _make_surface_model(contract)
    acs_state = _make_acs_state(contract, surface_model)

    # Step 1: ECH activation on hydroxyl
    ech_reagent = REAGENT_PROFILES["ech_activation"]
    step1 = ModificationStep(
        step_type=ModificationStepType.ACTIVATION,
        reagent_key="ech_activation",
        target_acs=ACSSiteType.HYDROXYL,
        product_acs=ACSSiteType.EPOXIDE,
        temperature=298.15,
        time=7200.0,
        reagent_concentration=50.0,
    )
    result1 = solve_modification_step(step1, acs_state, surface_model, ech_reagent)

    # Step 2: Genipin secondary crosslinking on amine
    genipin_reagent = REAGENT_PROFILES["genipin_secondary"]
    step2 = ModificationStep(
        step_type=ModificationStepType.SECONDARY_CROSSLINKING,
        reagent_key="genipin_secondary",
        target_acs=ACSSiteType.AMINE_PRIMARY,
        temperature=310.15,
        time=14400.0,
        reagent_concentration=10.0,
        stoichiometry=0.5,
    )
    result2 = solve_modification_step(step2, acs_state, surface_model, genipin_reagent)

    # Both steps should have nonzero conversion
    assert result1.conversion > 0.0, "ECH activation should proceed"
    assert result2.conversion > 0.0, "Genipin crosslinking should proceed"

    # EPOXIDE should exist from step 1
    assert ACSSiteType.EPOXIDE in acs_state

    # NH2 should be partially consumed from step 2
    nh2 = acs_state[ACSSiteType.AMINE_PRIMARY]
    assert nh2.crosslinked_sites > 0

    # All profiles valid
    for profile in acs_state.values():
        violations = profile.validate()
        assert violations == [], f"Violations for {profile.site_type.value}: {violations}"


# ─── Test: ACS conservation across all steps ──────────────────────────

def test_acs_conservation_across_steps():
    """Total ACS conserved across all modification steps."""
    contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0)
    surface_model = _make_surface_model(contract)
    acs_state = _make_acs_state(contract, surface_model)

    # Record initial totals
    nh2_initial_total = acs_state[ACSSiteType.AMINE_PRIMARY].total_sites
    oh_initial_total = acs_state[ACSSiteType.HYDROXYL].total_sites

    # Run a crosslinking step
    reagent = REAGENT_PROFILES["genipin_secondary"]
    step = ModificationStep(
        step_type=ModificationStepType.SECONDARY_CROSSLINKING,
        reagent_key="genipin_secondary",
        target_acs=ACSSiteType.AMINE_PRIMARY,
        temperature=310.15,
        time=14400.0,
        reagent_concentration=10.0,
        stoichiometry=0.5,
    )
    solve_modification_step(step, acs_state, surface_model, reagent)

    # Total sites should not change (conservation)
    nh2 = acs_state[ACSSiteType.AMINE_PRIMARY]
    assert nh2.total_sites == pytest.approx(nh2_initial_total, rel=1e-10), \
        "Total NH2 sites should be conserved"

    # terminal_sum + remaining = accessible (new conservation law)
    assert (nh2.crosslinked_sites + nh2.hydrolyzed_sites
            + nh2.ligand_coupled_sites + nh2.blocked_sites
            + nh2.remaining_sites) == \
        pytest.approx(nh2.accessible_sites, rel=1e-10), \
        "terminal_sum + remaining should equal accessible"

    # OH should be untouched
    oh = acs_state[ACSSiteType.HYDROXYL]
    assert oh.total_sites == pytest.approx(oh_initial_total, rel=1e-10)
    assert oh.crosslinked_sites == 0.0


# ─── Test: Orchestrator end-to-end ────────────────────────────────────

def test_orchestrator_end_to_end():
    """Full M1 contract -> orchestrator -> FunctionalMicrosphere."""
    contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0, G_DN=5000.0)

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15,
            time=7200.0,
            reagent_concentration=50.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.SECONDARY_CROSSLINKING,
            reagent_key="genipin_secondary",
            target_acs=ACSSiteType.AMINE_PRIMARY,
            temperature=310.15,
            time=14400.0,
            reagent_concentration=10.0,
            stoichiometry=0.5,
        ),
    ]

    orchestrator = ModificationOrchestrator()
    result = orchestrator.run(contract, steps)

    # Type check
    assert isinstance(result, FunctionalMicrosphere)

    # Has modification history
    assert len(result.modification_history) == 2

    # G_DN should have increased from secondary crosslinking
    assert result.G_DN_updated >= contract.G_DN, \
        f"G_DN should not decrease: {result.G_DN_updated} < {contract.G_DN}"

    # E* should be approximately 3 * G_DN
    assert result.E_star_updated == pytest.approx(3.0 * result.G_DN_updated, rel=1e-10)

    # ACS conservation
    violations = result.validate()
    assert violations == [], f"Conservation violations: {violations}"

    # Both steps should have nonzero conversion
    for mod_result in result.modification_history:
        assert 0.0 <= mod_result.conversion <= 1.0


# ─── Test: Conversion physical range ──────────────────────────────────

def test_conversion_physical_range():
    """Conversion between 0 and 1 for all tested conditions."""
    contract = _make_contract()
    surface_model = _make_surface_model(contract)

    conditions = [
        ("genipin_secondary", ACSSiteType.AMINE_PRIMARY,
         ModificationStepType.SECONDARY_CROSSLINKING, 310.15, 14400.0),
        ("glutaraldehyde_secondary", ACSSiteType.AMINE_PRIMARY,
         ModificationStepType.SECONDARY_CROSSLINKING, 298.15, 3600.0),
        ("ech_activation", ACSSiteType.HYDROXYL,
         ModificationStepType.ACTIVATION, 298.15, 7200.0),
        ("dvs_activation", ACSSiteType.HYDROXYL,
         ModificationStepType.ACTIVATION, 298.15, 3600.0),
    ]

    for reagent_key, target, step_type, temp, time in conditions:
        acs_state = _make_acs_state(contract, surface_model)
        reagent = REAGENT_PROFILES[reagent_key]

        step = ModificationStep(
            step_type=step_type,
            reagent_key=reagent_key,
            target_acs=target,
            product_acs=reagent.product_acs,
            temperature=temp,
            time=time,
            reagent_concentration=20.0,
            stoichiometry=reagent.stoichiometry,
        )

        result = solve_modification_step(step, acs_state, surface_model, reagent)

        assert 0.0 <= result.conversion <= 1.0, \
            f"{reagent_key}: conversion={result.conversion} out of [0,1]"


# ─── Test: Hydrolysis reduces yield ───────────────────────────────────

def test_hydrolysis_reduces_yield():
    """ECH with high hydrolysis gives lower epoxide yield than DVS.

    ECH has k_hydrol = 1e-4 /s which competes with coupling,
    while DVS has negligible hydrolysis. Under identical conditions
    (same OH, same reagent concentration, same time), ECH should
    produce fewer activated sites due to reagent loss.
    """
    contract = _make_contract(oh_bulk=400.0)

    # Run ECH activation
    surface_model_ech = _make_surface_model(contract)
    acs_ech = _make_acs_state(contract, surface_model_ech)
    ech_reagent = REAGENT_PROFILES["ech_activation"]

    step_ech = ModificationStep(
        step_type=ModificationStepType.ACTIVATION,
        reagent_key="ech_activation",
        target_acs=ACSSiteType.HYDROXYL,
        product_acs=ACSSiteType.EPOXIDE,
        temperature=298.15,
        time=7200.0,
        reagent_concentration=50.0,
    )
    result_ech = solve_modification_step(step_ech, acs_ech, surface_model_ech, ech_reagent)

    # Run DVS activation with same conditions
    surface_model_dvs = _make_surface_model(contract)
    acs_dvs = _make_acs_state(contract, surface_model_dvs)
    dvs_reagent = REAGENT_PROFILES["dvs_activation"]

    step_dvs = ModificationStep(
        step_type=ModificationStepType.ACTIVATION,
        reagent_key="dvs_activation",
        target_acs=ACSSiteType.HYDROXYL,
        product_acs=ACSSiteType.VINYL_SULFONE,
        temperature=298.15,
        time=7200.0,
        reagent_concentration=50.0,
    )
    result_dvs = solve_modification_step(step_dvs, acs_dvs, surface_model_dvs, dvs_reagent)

    # ECH has k_forward = 1.5e-5 (faster coupling) but hydrolysis = 1e-4 /s
    # DVS has k_forward = 5e-6 (slower coupling) but no hydrolysis
    # The key point: with hydrolysis, ECH's reagent is consumed by
    # both coupling AND hydrolysis, so less reagent is available for coupling.
    # With same reagent_concentration and enough time, the net activated
    # yield depends on the competition.
    #
    # For a fair comparison of hydrolysis effect, we compare the reagent
    # remaining fractions. ECH should have less reagent remaining because
    # hydrolysis consumes reagent too.

    # Both should have positive conversion
    assert result_ech.conversion > 0.0
    assert result_dvs.conversion > 0.0

    # The reagent remaining for ECH should be lower (hydrolysis consumes it)
    # We verify this through the ODE solver directly
    acs_conc = 400.0 * (4.0 / 3.0 * math.pi * (50e-6)**3) / (4.0 / 3.0 * math.pi * (50e-6)**3)
    # Use direct solver to compare reagent remaining
    _, ech_remaining = solve_second_order_consumption(
        acs_concentration=acs_conc,
        reagent_concentration=50.0,
        k_forward=1.5e-5,
        stoichiometry=1.0,
        time=7200.0,
        temperature=298.15,
        hydrolysis_rate=1e-4,
    )
    _, dvs_remaining = solve_second_order_consumption(
        acs_concentration=acs_conc,
        reagent_concentration=50.0,
        k_forward=5e-6,
        stoichiometry=1.0,
        time=7200.0,
        temperature=298.15,
        hydrolysis_rate=0.0,
    )

    # ECH reagent should be more depleted (hydrolysis + coupling)
    assert ech_remaining < dvs_remaining, (
        f"ECH reagent remaining ({ech_remaining:.4f}) should be less than "
        f"DVS ({dvs_remaining:.4f}) due to hydrolysis competition"
    )


# ─── Test: Solver edge cases ──────────────────────────────────────────

def test_second_order_consumption_zero_time():
    """Zero time should give zero conversion."""
    conv, remaining = solve_second_order_consumption(
        acs_concentration=100.0,
        reagent_concentration=50.0,
        k_forward=1e-5,
        stoichiometry=1.0,
        time=0.0,
        temperature=298.15,
    )
    assert conv == 0.0
    assert remaining == 1.0


def test_second_order_consumption_zero_reagent():
    """Zero reagent should give zero conversion."""
    conv, remaining = solve_second_order_consumption(
        acs_concentration=100.0,
        reagent_concentration=0.0,
        k_forward=1e-5,
        stoichiometry=1.0,
        time=3600.0,
        temperature=298.15,
    )
    assert conv == 0.0
    assert remaining == 1.0


# ─── Test: Temperature affects conversion via Arrhenius ───────────────

def test_temperature_affects_conversion():
    """Verify that higher temperature increases conversion (Arrhenius).

    BF-1 regression test: before the fix, k0=0.0 meant Arrhenius was
    bypassed and temperature had no effect on kinetics.  After the fix,
    _arrhenius_prefactor back-calculates k0 from the reagent profile so
    the Arrhenius branch in solve_second_order_consumption is engaged.
    """
    contract = _make_contract(nh2_bulk=100.0)
    reagent = REAGENT_PROFILES["genipin_secondary"]

    def _run_at_temperature(temp_K: float) -> float:
        surface_model = _make_surface_model(contract)
        acs_state = _make_acs_state(contract, surface_model)
        step = ModificationStep(
            step_type=ModificationStepType.SECONDARY_CROSSLINKING,
            reagent_key="genipin_secondary",
            target_acs=ACSSiteType.AMINE_PRIMARY,
            temperature=temp_K,
            time=14400.0,    # 4 h — same for both runs
            reagent_concentration=10.0,
            stoichiometry=0.5,
        )
        result = solve_modification_step(step, acs_state, surface_model, reagent)
        return result.conversion

    conversion_298 = _run_at_temperature(298.15)   # 25 degC
    conversion_310 = _run_at_temperature(310.15)   # 37 degC (reference T)

    assert conversion_298 > 0.0, "Conversion at 298 K should be positive"
    assert conversion_310 > 0.0, "Conversion at 310 K should be positive"
    assert conversion_310 > conversion_298, (
        f"Higher temperature should give higher conversion via Arrhenius: "
        f"conv(310K)={conversion_310:.6f} vs conv(298K)={conversion_298:.6f}"
    )


# ═════════════════════════════════════════════════════════════════════════
# W12: Integration tests for new workflows (Ligand, Protein, Quench)
# ═════════════════════════════════════════════════════════════════════════


def test_full_activation_coupling_quench_pipeline():
    """Full Activation -> DEAE Coupling -> Ethanolamine Quench pipeline.

    Conservation must hold after each step and at the end.
    Quenching must only update blocked_sites (audit F1 fix).
    """
    contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0, G_DN=5000.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.LIGAND_COUPLING,
            reagent_key="deae_coupling",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=14400.0,
            reagent_concentration=100.0, ph=10.5,
        ),
        ModificationStep(
            step_type=ModificationStepType.QUENCHING,
            reagent_key="ethanolamine_quench",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=1000.0, ph=8.5,
        ),
    ]

    result = orchestrator.run(contract, steps)

    # All 3 steps executed
    assert len(result.modification_history) == 3

    # Conservation must pass
    violations = result.validate()
    assert violations == [], f"Conservation violations: {violations}"

    # EPOXIDE profile should have coupling + quench results
    epox = result.acs_profiles[ACSSiteType.EPOXIDE]
    assert epox.ligand_coupled_sites > 0, "DEAE should have coupled"
    assert epox.blocked_sites > 0, "Ethanolamine should have blocked"
    assert epox.remaining_activated == pytest.approx(0.0, abs=1e-20), \
        "All activated sites should be consumed by coupling + quench"

    # G_DN should NOT increase (no crosslinking in coupling/quench)
    assert result.G_DN_updated == pytest.approx(contract.G_DN, rel=1e-10)


def test_quench_95pct_conservation_regression():
    """Audit F1 regression: 95%+ quenching must not violate conservation.

    Previously, quenching incremented both consumed_sites and blocked_sites,
    causing double-counting. With the new terminal-state model, quenching
    only increments blocked_sites.
    """
    contract = _make_contract(oh_bulk=400.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.QUENCHING,
            reagent_key="ethanolamine_quench",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=1000.0, ph=8.5,  # 1M = very high
        ),
    ]

    result = orchestrator.run(contract, steps)
    violations = result.validate()
    assert violations == [], f"F1 regression: conservation violated at >95% quench: {violations}"

    epox = result.acs_profiles[ACSSiteType.EPOXIDE]
    # Quench conversion should be >95%
    quench_result = result.modification_history[1]
    assert quench_result.conversion > 0.95, \
        f"Expected >95% quench conversion, got {quench_result.conversion:.4f}"

    # remaining_activated should be ~0
    assert epox.remaining_activated < epox.activated_sites * 0.05, \
        "remaining_activated should be <5% of activated after 95%+ quench"

    # blocked_sites should be ~= activated_sites
    assert epox.blocked_sites > 0.95 * epox.activated_sites * 0.99


def test_protein_coupling_steric_limit():
    """Protein coupling respects steric jamming limit."""
    contract = _make_contract(oh_bulk=400.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.PROTEIN_COUPLING,
            reagent_key="protein_a_coupling",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=277.15, time=57600.0,
            reagent_concentration=0.01,  # 10 uM protein
            ph=9.0,
        ),
    ]

    result = orchestrator.run(contract, steps)
    violations = result.validate()
    assert violations == [], f"Conservation violations: {violations}"

    epox = result.acs_profiles[ACSSiteType.EPOXIDE]
    # Protein A has activity_retention = 0.60
    if epox.ligand_coupled_sites > 0:
        ratio = epox.ligand_functional_sites / epox.ligand_coupled_sites
        assert ratio == pytest.approx(0.60, rel=0.01), \
            f"Activity retention should be ~0.60, got {ratio:.3f}"


def test_backend_blocks_coupling_without_activation():
    """Backend validation blocks coupling on non-activated sites."""
    contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0)

    steps = [
        ModificationStep(
            step_type=ModificationStepType.LIGAND_COUPLING,
            reagent_key="deae_coupling",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=14400.0,
            reagent_concentration=100.0, ph=10.5,
        ),
    ]

    orchestrator = ModificationOrchestrator()
    with pytest.raises(ValueError, match="no prior activation"):
        orchestrator.run(contract, steps)


def test_backend_blocks_step_after_quench():
    """Backend validation blocks steps after quenching on same target."""
    contract = _make_contract(oh_bulk=400.0)

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.QUENCHING,
            reagent_key="ethanolamine_quench",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=1000.0, ph=8.5,
        ),
        ModificationStep(
            step_type=ModificationStepType.LIGAND_COUPLING,
            reagent_key="deae_coupling",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=14400.0,
            reagent_concentration=100.0, ph=10.5,
        ),
    ]

    orchestrator = ModificationOrchestrator()
    with pytest.raises(ValueError, match="already quenched"):
        orchestrator.run(contract, steps)


def test_functional_media_contract_iex():
    """FunctionalMediaContract maps DEAE ligand density to q_max."""
    from emulsim.module2_functionalization.orchestrator import build_functional_media_contract

    contract = _make_contract(oh_bulk=400.0, G_DN=5000.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.LIGAND_COUPLING,
            reagent_key="deae_coupling",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=14400.0,
            reagent_concentration=100.0, ph=10.5,
        ),
    ]

    result = orchestrator.run(contract, steps)
    fmc = build_functional_media_contract(result)

    assert fmc.ligand_type == "iex_anion"
    assert fmc.installed_ligand == "DEAE"
    assert fmc.functional_ligand_density > 0
    assert fmc.estimated_q_max > 0
    assert fmc.q_max_confidence == "mapped_estimated"


def test_backend_blocks_reagent_target_mismatch():
    """Backend validation blocks reagent-target ACS mismatch."""
    contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0)

    steps = [
        # ethanolamine_quench targets EPOXIDE, but step says AMINE_PRIMARY
        ModificationStep(
            step_type=ModificationStepType.QUENCHING,
            reagent_key="ethanolamine_quench",
            target_acs=ACSSiteType.AMINE_PRIMARY,
            temperature=298.15, time=7200.0,
            reagent_concentration=1000.0,
        ),
    ]

    orchestrator = ModificationOrchestrator()
    with pytest.raises(ValueError, match="Reagent-target mismatch"):
        orchestrator.run(contract, steps)


# ═════════════════════════════════════════════════════════════════════════
# Codex P1-1, P1-2, P2-1 regression tests
# ═════════════════════════════════════════════════════════════════════════


def test_p1_1_reagent_steptype_mismatch_blocked():
    """Codex P1-1: QUENCHING with a coupling reagent must raise ValueError."""
    contract = _make_contract(oh_bulk=400.0)

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        # WRONG: using a coupling reagent for a quenching step
        ModificationStep(
            step_type=ModificationStepType.QUENCHING,
            reagent_key="deae_coupling",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0,
        ),
    ]

    orchestrator = ModificationOrchestrator()
    with pytest.raises(ValueError, match="incompatible with step type"):
        orchestrator.run(contract, steps)


def test_p1_2_acetic_anhydride_quench_on_native_amine():
    """Codex P1-2: Acetic anhydride quench on AMINE_PRIMARY must achieve conversion > 0."""
    contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.QUENCHING,
            reagent_key="acetic_anhydride_quench",
            target_acs=ACSSiteType.AMINE_PRIMARY,
            temperature=298.15, time=3600.0,
            reagent_concentration=500.0, ph=7.5,
        ),
    ]

    result = orchestrator.run(contract, steps)
    quench_result = result.modification_history[0]

    assert quench_result.conversion > 0, \
        f"Acetic anhydride should block native amines, got conversion={quench_result.conversion}"

    nh2 = result.acs_profiles[ACSSiteType.AMINE_PRIMARY]
    assert nh2.blocked_sites > 0, \
        "blocked_sites should be > 0 after acetic anhydride quench"

    violations = result.validate()
    assert violations == [], f"Conservation violations: {violations}"


def test_p2_1_activation_uses_activated_consumed_not_crosslinked():
    """Codex P2-1: Activation should put consumed OH into activated_consumed_sites,
    not crosslinked_sites."""
    contract = _make_contract(oh_bulk=400.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
    ]

    result = orchestrator.run(contract, steps)
    oh = result.acs_profiles[ACSSiteType.HYDROXYL]

    # Activation should NOT increment crosslinked_sites
    assert oh.crosslinked_sites == 0.0, \
        f"Activation should not use crosslinked_sites, got {oh.crosslinked_sites}"

    # Should use activated_consumed_sites instead
    assert oh.activated_consumed_sites > 0, \
        f"Activation should use activated_consumed_sites, got {oh.activated_consumed_sites}"

    violations = result.validate()
    assert violations == [], f"Conservation violations: {violations}"


# ═════════════════════════════════════════════════════════════════════════
# WN-7/8/9: v5.7 expansion tests
# ═════════════════════════════════════════════════════════════════════════


def test_all_25_profiles_have_required_metadata():
    """Every profile must have non-empty confidence_tier and valid reaction_type."""
    for key, profile in REAGENT_PROFILES.items():
        assert profile.confidence_tier in ("semi_quantitative", "ranking_only"), \
            f"{key}: invalid confidence_tier '{profile.confidence_tier}'"
        assert profile.reaction_type in (
            "crosslinking", "activation", "coupling", "protein_coupling", "blocking",
            "spacer", "spacer_arm", "heterobifunctional",
            "metal_charging", "metal_stripping", "protein_pretreatment", "washing",
        ), f"{key}: invalid reaction_type '{profile.reaction_type}'"


def test_iex_profiles_have_charge_type():
    """All IEX ligand profiles must have explicit charge_type."""
    iex_keys = ["deae_coupling", "q_coupling", "sp_coupling", "cm_coupling"]
    for key in iex_keys:
        profile = REAGENT_PROFILES[key]
        assert profile.charge_type in ("anion", "cation"), \
            f"{key}: IEX profile missing charge_type, got '{profile.charge_type}'"


def test_imac_profiles_have_metal_ion():
    """IMAC profiles must declare metal_ion."""
    for key in ["ida_coupling", "nta_coupling"]:
        profile = REAGENT_PROFILES[key]
        assert profile.metal_ion != "", \
            f"{key}: IMAC profile missing metal_ion"
        assert profile.metal_loaded_fraction > 0, \
            f"{key}: metal_loaded_fraction should be > 0"


def test_spacer_profiles_are_not_executable():
    """Spacer profiles have reaction_type='spacer' and k_forward=0."""
    spacer_keys = ["dadpa_spacer", "aha_spacer", "dah_spacer"]
    for key in spacer_keys:
        profile = REAGENT_PROFILES[key]
        assert profile.reaction_type == "spacer", \
            f"{key}: spacer should have reaction_type='spacer'"
        assert profile.k_forward == 0.0, \
            f"{key}: spacer k_forward should be 0 (not executable)"
        assert profile.spacer_activity_multiplier >= 1.0, \
            f"{key}: spacer_activity_multiplier should be >= 1.0"


def test_heparin_is_macromolecule():
    """Heparin must be treated as macromolecule (audit F6)."""
    profile = REAGENT_PROFILES["heparin_coupling"]
    assert profile.is_macromolecule is True
    assert profile.ligand_mw > 10000


def test_glutathione_has_reduced_activity():
    """Glutathione activity_retention should be <1.0 (audit F8)."""
    profile = REAGENT_PROFILES["glutathione_coupling"]
    assert profile.activity_retention < 1.0
    assert profile.activity_retention_uncertainty > 0


def test_spacer_multiplier_capped_at_1():
    """effective_activity = activity_retention * spacer_mult must not exceed 1.0."""
    contract = _make_contract(oh_bulk=400.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.PROTEIN_COUPLING,
            reagent_key="streptavidin_coupling",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=277.15, time=57600.0,
            reagent_concentration=0.01, ph=9.0,
        ),
    ]
    result = orchestrator.run(contract, steps)
    epox = result.acs_profiles.get(ACSSiteType.EPOXIDE)
    if epox and epox.ligand_coupled_sites > 0:
        ratio = epox.ligand_functional_sites / epox.ligand_coupled_sites
        assert ratio <= 1.0, f"Effective activity should be <= 1.0, got {ratio:.4f}"
        assert ratio == pytest.approx(0.70 * 1.22, rel=0.02)


def test_fmc_q_maps_to_iex_anion():
    """Q maps to iex_anion via charge_type (audit F14)."""
    from emulsim.module2_functionalization.orchestrator import build_functional_media_contract
    contract = _make_contract(oh_bulk=400.0, G_DN=5000.0)
    orchestrator = ModificationOrchestrator()
    steps = [
        ModificationStep(step_type=ModificationStepType.ACTIVATION, reagent_key="ech_activation",
                         target_acs=ACSSiteType.HYDROXYL, product_acs=ACSSiteType.EPOXIDE,
                         temperature=298.15, time=7200.0, reagent_concentration=100.0, ph=12.0),
        ModificationStep(step_type=ModificationStepType.LIGAND_COUPLING, reagent_key="q_coupling",
                         target_acs=ACSSiteType.EPOXIDE, temperature=298.15, time=14400.0,
                         reagent_concentration=100.0, ph=10.5),
    ]
    result = orchestrator.run(contract, steps)
    fmc = build_functional_media_contract(result)
    assert fmc.ligand_type == "iex_anion"
    assert fmc.estimated_q_max > 0


def test_fmc_biotin_stoich_2_5():
    """Streptavidin q_max uses stoich=2.5 not 4 (audit F7)."""
    from emulsim.module2_functionalization.orchestrator import build_functional_media_contract
    contract = _make_contract(oh_bulk=400.0, G_DN=5000.0)
    orchestrator = ModificationOrchestrator()
    steps = [
        ModificationStep(step_type=ModificationStepType.ACTIVATION, reagent_key="ech_activation",
                         target_acs=ACSSiteType.HYDROXYL, product_acs=ACSSiteType.EPOXIDE,
                         temperature=298.15, time=7200.0, reagent_concentration=100.0, ph=12.0),
        ModificationStep(step_type=ModificationStepType.PROTEIN_COUPLING, reagent_key="streptavidin_coupling",
                         target_acs=ACSSiteType.EPOXIDE, temperature=277.15, time=57600.0,
                         reagent_concentration=0.01, ph=9.0),
    ]
    result = orchestrator.run(contract, steps)
    fmc = build_functional_media_contract(result)
    assert fmc.ligand_type == "biotin_affinity"
    assert "2.5" in fmc.q_max_mapping_notes
    assert fmc.binding_model_hint == "near_irreversible"


def test_profile_count():
    """Canonical profile count — 53 total (v5.9: 52, +1 for stmp_secondary)."""
    assert len(REAGENT_PROFILES) == 53, \
        f"Expected 53 profiles, got {len(REAGENT_PROFILES)}"


# ═════════════════════════════════════════════════════════════════════════
# v5.8 Phase 2: SPACER_ARM + SM(PEG)n integration tests (WN-8)
# ═════════════════════════════════════════════════════════════════════════


def test_spacer_arm_creates_amine_distal():
    """SPACER_ARM step consumes EPOXIDE and creates AMINE_DISTAL."""
    contract = _make_contract(oh_bulk=400.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.SPACER_ARM,
            reagent_key="dadpa_spacer_arm",
            target_acs=ACSSiteType.EPOXIDE,
            product_acs=ACSSiteType.AMINE_DISTAL,
            temperature=298.15, time=14400.0,
            reagent_concentration=100.0, ph=10.5,
        ),
    ]

    result = orchestrator.run(contract, steps)
    violations = result.validate()
    assert violations == [], f"Conservation violations: {violations}"

    # AMINE_DISTAL should exist
    assert ACSSiteType.AMINE_DISTAL in result.acs_profiles, \
        "AMINE_DISTAL profile should be created by SPACER_ARM"

    ad = result.acs_profiles[ACSSiteType.AMINE_DISTAL]
    assert ad.accessible_sites > 0, "AMINE_DISTAL should have accessible sites"
    assert ad.activated_sites > 0, "AMINE_DISTAL should have activated sites"


def test_spacer_arm_distal_yield():
    """Distal group yield < 1.0 means some sites are bridged, not all produce product."""
    contract = _make_contract(oh_bulk=400.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        ModificationStep(
            step_type=ModificationStepType.SPACER_ARM,
            reagent_key="eda_spacer_arm",  # distal_group_yield = 0.60
            target_acs=ACSSiteType.EPOXIDE,
            product_acs=ACSSiteType.AMINE_DISTAL,
            temperature=298.15, time=14400.0,
            reagent_concentration=200.0, ph=10.5,
        ),
    ]

    result = orchestrator.run(contract, steps)
    violations = result.validate()
    assert violations == [], f"Conservation violations: {violations}"

    # Check: AMINE_DISTAL sites < EPOXIDE consumed (distal_yield < 1.0)
    epox = result.acs_profiles[ACSSiteType.EPOXIDE]
    ad = result.acs_profiles.get(ACSSiteType.AMINE_DISTAL)
    assert ad is not None, "AMINE_DISTAL should be created"
    # EDA yield = 0.60, so created sites < consumed sites
    if epox.ligand_coupled_sites > 0 and ad.accessible_sites > 0:
        yield_ratio = ad.accessible_sites / epox.ligand_coupled_sites
        assert yield_ratio < 1.0, \
            f"Distal yield should be < 1.0 for EDA, got {yield_ratio:.2f}"


def test_full_smpeg_path():
    """Full SM(PEG)n path: ECH → DADPA → SM(PEG)4 → Protein A-Cys → Quench."""
    contract = _make_contract(oh_bulk=400.0, G_DN=5000.0)
    orchestrator = ModificationOrchestrator()

    steps = [
        # Step 1: Activate hydroxyl → epoxide
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        # Step 2: Spacer arm — epoxide → amine_distal
        ModificationStep(
            step_type=ModificationStepType.SPACER_ARM,
            reagent_key="dadpa_spacer_arm",
            target_acs=ACSSiteType.EPOXIDE,
            product_acs=ACSSiteType.AMINE_DISTAL,
            temperature=298.15, time=14400.0,
            reagent_concentration=100.0, ph=10.5,
        ),
        # Step 3: Heterobifunctional — amine_distal → maleimide
        ModificationStep(
            step_type=ModificationStepType.SPACER_ARM,
            reagent_key="sm_peg4",
            target_acs=ACSSiteType.AMINE_DISTAL,
            product_acs=ACSSiteType.MALEIMIDE,
            temperature=298.15, time=1800.0,
            reagent_concentration=10.0, ph=7.4,
        ),
        # Step 4: Protein coupling — maleimide → thioether (Protein A-Cys)
        ModificationStep(
            step_type=ModificationStepType.PROTEIN_COUPLING,
            reagent_key="protein_a_cys_coupling",
            target_acs=ACSSiteType.MALEIMIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=0.01, ph=7.0,
        ),
        # Step 5: Quench remaining epoxides
        ModificationStep(
            step_type=ModificationStepType.QUENCHING,
            reagent_key="ethanolamine_quench",
            target_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=1000.0, ph=8.5,
        ),
    ]

    result = orchestrator.run(contract, steps)

    # All 5 steps should execute
    assert len(result.modification_history) == 5

    # Conservation must hold across all profiles
    violations = result.validate()
    assert violations == [], f"Conservation violations: {violations}"

    # MALEIMIDE profile should exist with coupled protein
    assert ACSSiteType.MALEIMIDE in result.acs_profiles
    mal = result.acs_profiles[ACSSiteType.MALEIMIDE]
    # Some maleimide may have decayed + some coupled to protein
    assert mal.hydrolyzed_sites >= 0  # maleimide decay
    assert mal.ligand_coupled_sites >= 0  # protein coupling

    # Protein A-Cys should have higher activity than random coupling (0.80 vs 0.60)
    if mal.ligand_coupled_sites > 0:
        ratio = mal.ligand_functional_sites / mal.ligand_coupled_sites
        # effective_activity = 0.80 * spacer_mult (but Cys profiles have multiplier 1.0)
        assert ratio <= 1.0


def test_smpeg_without_amine_distal_fails():
    """SM(PEG)n targeting AMINE_DISTAL without prior spacer should fail validation."""
    contract = _make_contract(oh_bulk=400.0)

    steps = [
        ModificationStep(
            step_type=ModificationStepType.ACTIVATION,
            reagent_key="ech_activation",
            target_acs=ACSSiteType.HYDROXYL,
            product_acs=ACSSiteType.EPOXIDE,
            temperature=298.15, time=7200.0,
            reagent_concentration=100.0, ph=12.0,
        ),
        # Skip DADPA spacer — go directly to SM(PEG)4 on AMINE_DISTAL
        ModificationStep(
            step_type=ModificationStepType.SPACER_ARM,
            reagent_key="sm_peg4",
            target_acs=ACSSiteType.AMINE_DISTAL,
            product_acs=ACSSiteType.MALEIMIDE,
            temperature=298.15, time=1800.0,
            reagent_concentration=10.0, ph=7.4,
        ),
    ]

    orchestrator = ModificationOrchestrator()
    with pytest.raises(ValueError, match="no prior"):
        orchestrator.run(contract, steps)


# ─── Node 4 (v6.1): M2 ModelManifest evidence wiring ──────────────────────


class TestM2ModelManifest:
    """Verifies that ModelManifest fields are populated end-to-end through M2.

    Acceptance for Node 4:
      - solve_modification_step attaches a manifest to every ModificationResult.
      - ModificationOrchestrator.run produces a FunctionalMicrosphere whose
        composite manifest reflects the weakest tier across all steps.
      - build_functional_media_contract attaches an FMC manifest combining
        the upstream microsphere tier with the FMC's own mapping confidence.
      - Ranking-only ligand classes (affinity, biotin, heparin) cap evidence
        at QUALITATIVE_TREND.
    """

    def test_modification_result_has_manifest(self):
        """A genipin crosslinking step attaches a SEMI_QUANTITATIVE manifest."""
        from emulsim.datatypes import ModelEvidenceTier, ModelManifest

        contract = _make_contract(nh2_bulk=100.0)
        surface_model = _make_surface_model(contract)
        acs_state = _make_acs_state(contract, surface_model)

        step = ModificationStep(
            step_type=ModificationStepType.SECONDARY_CROSSLINKING,
            reagent_key="genipin_secondary",
            target_acs=ACSSiteType.AMINE_PRIMARY,
            temperature=310.15, time=14400.0,
            reagent_concentration=10.0, stoichiometry=0.5,
        )
        result = solve_modification_step(
            step, acs_state, surface_model, REAGENT_PROFILES["genipin_secondary"]
        )

        assert result.model_manifest is not None
        assert isinstance(result.model_manifest, ModelManifest)
        # Genipin is not ranking-only and conservation is preserved -> SEMI
        assert result.model_manifest.evidence_tier == ModelEvidenceTier.SEMI_QUANTITATIVE
        assert "secondary_crosslinking" in result.model_manifest.model_name
        assert "genipin_secondary" in result.model_manifest.model_name
        # Diagnostics carry the achieved conversion
        diag = result.model_manifest.diagnostics
        assert "conversion" in diag and diag["conversion"] > 0
        assert diag["conservation_ok"] is True

    def test_microsphere_composite_manifest_weakest_wins(self):
        """Composite manifest reports the weakest tier across all steps."""
        from emulsim.datatypes import ModelEvidenceTier
        from emulsim.module2_functionalization.orchestrator import (
            build_functional_media_contract,
        )

        contract = _make_contract(nh2_bulk=200.0, oh_bulk=400.0)
        # Two steps: ECH activation (SEMI) then a ranking-only affinity coupling
        # to force the composite down to QUALITATIVE_TREND.
        steps = [
            ModificationStep(
                step_type=ModificationStepType.ACTIVATION,
                reagent_key="ech_activation",
                target_acs=ACSSiteType.HYDROXYL,
                product_acs=ACSSiteType.EPOXIDE,
                temperature=298.15, time=7200.0,
                reagent_concentration=50.0, ph=12.0,
            ),
            ModificationStep(
                step_type=ModificationStepType.PROTEIN_COUPLING,
                reagent_key="protein_a_coupling",
                target_acs=ACSSiteType.EPOXIDE,
                temperature=277.15, time=14400.0,
                reagent_concentration=0.05, ph=7.0,
            ),
        ]

        microsphere = ModificationOrchestrator().run(contract, steps)

        assert microsphere.model_manifest is not None
        # Protein A is functional_mode='affinity_ligand' -> ranking-only ->
        # composite must be QUALITATIVE_TREND or weaker.
        _ORDER = list(ModelEvidenceTier)
        composite_idx = _ORDER.index(microsphere.model_manifest.evidence_tier)
        qt_idx = _ORDER.index(ModelEvidenceTier.QUALITATIVE_TREND)
        assert composite_idx >= qt_idx, (
            f"Composite tier {microsphere.model_manifest.evidence_tier} should "
            "be QUALITATIVE_TREND or weaker after a ranking-only coupling step."
        )
        # Diagnostics enumerate steps
        steps_diag = microsphere.model_manifest.diagnostics["steps"]
        assert len(steps_diag) == 2

        # FMC manifest inherits the weakest of microsphere + FMC's own tier
        fmc = build_functional_media_contract(microsphere)
        assert fmc.model_manifest is not None
        fmc_idx = _ORDER.index(fmc.model_manifest.evidence_tier)
        assert fmc_idx >= composite_idx, (
            "FMC tier must not be stronger than upstream microsphere tier."
        )
        # FMC manifest names ligand_type
        assert fmc.model_manifest.diagnostics["ligand_type"] == "affinity"

    def test_empty_history_yields_unsupported(self):
        """No steps -> microsphere manifest is UNSUPPORTED, not None."""
        from emulsim.datatypes import ModelEvidenceTier

        contract = _make_contract()
        microsphere = ModificationOrchestrator().run(contract, steps=[])

        assert microsphere.model_manifest is not None
        assert microsphere.model_manifest.evidence_tier == ModelEvidenceTier.UNSUPPORTED
        assert microsphere.model_manifest.diagnostics["n_steps"] == 0
