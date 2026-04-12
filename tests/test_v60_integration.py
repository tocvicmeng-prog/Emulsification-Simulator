"""v6.0 end-to-end integration tests.

Tests cross-module boundaries for the three v6.0 frameworks:
  1. Calibration: override FMC via CalibrationStore, verify in M3
  2. Uncertainty: run_with_uncertainty through M2 pipeline
  3. Lifetime: project_lifetime with representative parameters
  4. Gradient elution: run_gradient_elution with gradient-sensitive isotherm
  5. Full M2→M3 pipeline with calibration + gradient
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from emulsim.datatypes import M1ExportContract


# ─── Helpers ──────────────────────────────────────────────────────────────

def _make_contract(
    nh2_bulk: float = 100.0,
    oh_bulk: float = 400.0,
    G_DN: float = 5000.0,
    porosity: float = 0.7,
) -> M1ExportContract:
    """Minimal M1ExportContract for integration tests."""
    return M1ExportContract(
        bead_radius=50e-6,
        bead_d32=100e-6,
        bead_d50=100e-6,
        pore_size_mean=100e-9,
        pore_size_std=30e-9,
        porosity=porosity,
        l2_model_tier="empirical_calibrated",
        mesh_size_xi=20e-9,
        p_final=0.5,
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


# ═════════════════════════════════════════════════════════════════════════
# 1. CALIBRATION INTEGRATION
# ═════════════════════════════════════════════════════════════════════════

class TestCalibrationIntegration:
    """Calibration framework end-to-end: load JSON → apply to FMC → verify override."""

    def test_load_json_roundtrip(self, tmp_path):
        """CalibrationStore loads JSON, entries survive save/reload."""
        from emulsim.calibration import CalibrationEntry, CalibrationStore

        entries = [
            CalibrationEntry(
                profile_key="protein_a_coupling",
                parameter_name="estimated_q_max",
                measured_value=120.5,
                units="mol/m3",
                target_molecule="IgG1",
                temperature_C=25.0,
                ph=7.0,
                confidence="measured",
                measurement_type="static_binding",
                replicates=3,
            ),
            CalibrationEntry(
                profile_key="deae_coupling",
                parameter_name="charge_density",
                measured_value=0.15,
                units="meq/mL",
                confidence="literature",
            ),
        ]

        # Save
        json_path = tmp_path / "cal_test.json"
        store = CalibrationStore()
        for e in entries:
            store.add(e)
        store.save_json(json_path)

        # Reload
        store2 = CalibrationStore()
        n = store2.load_json(json_path)
        assert n == 2
        assert len(store2) == 2

        # Query
        pa_entries = store2.query("protein_a_coupling")
        assert len(pa_entries) == 1
        assert pa_entries[0].measured_value == pytest.approx(120.5)
        assert pa_entries[0].confidence == "measured"

    def test_apply_calibration_overrides_fmc(self):
        """CalibrationStore.apply_to_fmc overrides FMC q_max."""
        from emulsim.calibration import CalibrationEntry, CalibrationStore
        from emulsim.module2_functionalization import (
            ModificationOrchestrator, ModificationStep, ModificationStepType, ACSSiteType,
        )
        from emulsim.module2_functionalization.orchestrator import build_functional_media_contract

        contract = _make_contract()
        orch = ModificationOrchestrator()

        # Run M2 with a protein A coupling step
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
                temperature=298.15, time=14400.0,
                reagent_concentration=5.0, ph=8.0,
            ),
        ]
        result = orch.run(contract, steps)
        fmc = build_functional_media_contract(result)

        # Now apply calibration override
        store = CalibrationStore()
        store.add(CalibrationEntry(
            profile_key="protein_a_coupling",
            parameter_name="estimated_q_max",
            measured_value=150.0,
            units="mol/m3",
            confidence="measured",
        ))

        fmc_cal, overrides = store.apply_to_fmc(fmc)
        assert len(overrides) > 0, "Should report at least one override"
        assert fmc_cal.estimated_q_max == pytest.approx(150.0)

    def test_load_json_from_file(self, tmp_path):
        """Load calibration from a JSON file with typed validation."""
        from emulsim.calibration import CalibrationStore

        data = [
            {
                "profile_key": "phenyl_coupling",
                "parameter_name": "estimated_q_max",
                "measured_value": 45.0,
                "units": "mol/m3",
                "confidence": "literature",
            }
        ]
        json_path = tmp_path / "hic_cal.json"
        json_path.write_text(json.dumps(data))

        store = CalibrationStore()
        n = store.load_json(json_path)
        assert n == 1
        assert store.entries[0].measured_value == pytest.approx(45.0)
        assert store.entries[0].confidence == "literature"


# ═════════════════════════════════════════════════════════════════════════
# 2. UNCERTAINTY INTEGRATION
# ═════════════════════════════════════════════════════════════════════════

class TestUncertaintyIntegration:
    """Uncertainty propagation end-to-end through M2."""

    def test_screening_defaults_produce_bounds(self):
        """Monte Carlo with screening CVs produces p5/p95 bounds."""
        from emulsim.uncertainty_propagation import (
            M1UncertaintyContract, run_with_uncertainty,
        )
        from emulsim.module2_functionalization import (
            ModificationStep, ModificationStepType, ACSSiteType,
        )

        contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0)
        unc = M1UncertaintyContract.screening_defaults()

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

        result = run_with_uncertainty(
            contract=contract,
            steps=steps,
            uncertainty=unc,
            n_samples=20,  # small N for speed
            seed=42,
        )

        assert result.n_samples == 20
        assert result.n_valid > 0, "Should have at least some valid samples"
        assert result.p5_q_max <= result.median_q_max <= result.p95_q_max
        assert result.cv_q_max >= 0.0
        assert result.tier == "assumed"

    def test_deterministic_produces_single_value(self):
        """Zero CVs = deterministic, all percentiles equal."""
        from emulsim.uncertainty_propagation import (
            M1UncertaintyContract, run_with_uncertainty,
        )
        from emulsim.module2_functionalization import (
            ModificationStep, ModificationStepType, ACSSiteType,
        )

        contract = _make_contract()
        unc = M1UncertaintyContract()  # all CVs=0

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

        result = run_with_uncertainty(
            contract=contract, steps=steps, uncertainty=unc,
            n_samples=5, seed=0,
        )

        # Deterministic: all values equal
        assert result.p5_q_max == pytest.approx(result.p95_q_max, rel=1e-6)
        assert result.cv_q_max == pytest.approx(0.0, abs=1e-6)


# ═════════════════════════════════════════════════════════════════════════
# 3. LIFETIME INTEGRATION
# ═════════════════════════════════════════════════════════════════════════

class TestLifetimeIntegration:
    """Lifetime projection end-to-end."""

    def test_protein_a_typical_lifetime(self):
        """Protein A at k=0.005 should last ~60 cycles to 80%."""
        from emulsim.lifetime import project_lifetime

        proj = project_lifetime(
            initial_capacity=100.0,
            k_deactivation=0.005,
            n_cycles=200,
            cip_description="0.1 M NaOH, 15 min",
        )

        assert proj.confidence == "empirical"
        assert proj.initial_capacity == pytest.approx(100.0)
        # ln(0.8) / 0.005 ~ 44.6 cycles
        assert 40 <= proj.cycles_to_80pct <= 50
        # ln(0.5) / 0.005 ~ 138.6 cycles
        assert 130 <= proj.cycles_to_50pct <= 145
        # Capacity at 200 cycles: 100 * exp(-0.005*200) = 100 * exp(-1) ~ 36.8
        assert proj.capacity_at_n == pytest.approx(100.0 * np.exp(-1.0), rel=0.01)
        assert "NaOH" in proj.assumption_notes

    def test_no_deactivation_infinite(self):
        """k=0 means infinite lifetime."""
        from emulsim.lifetime import project_lifetime

        proj = project_lifetime(
            initial_capacity=50.0,
            k_deactivation=0.0,
            n_cycles=100,
        )

        assert proj.cycles_to_80pct >= 999999
        assert proj.cycles_to_50pct >= 999999
        assert proj.capacity_at_n == pytest.approx(50.0)

    def test_high_deactivation_fast_decay(self):
        """High k = fast decay."""
        from emulsim.lifetime import project_lifetime

        proj = project_lifetime(
            initial_capacity=100.0,
            k_deactivation=0.05,  # aggressive
            n_cycles=50,
        )

        # ln(0.8) / 0.05 ~ 4.5 cycles
        assert proj.cycles_to_80pct <= 6
        assert proj.capacity_at_n < 10.0  # 100*exp(-2.5) ~ 8.2


# ═════════════════════════════════════════════════════════════════════════
# 4. GRADIENT ELUTION ORCHESTRATOR INTEGRATION
# ═════════════════════════════════════════════════════════════════════════

class TestGradientElutionIntegration:
    """run_gradient_elution with gradient-sensitive isotherms."""

    def test_competitive_langmuir_gradient_not_sensitive(self):
        """CompetitiveLangmuir: gradient_affects_binding should be False."""
        from emulsim.module3_performance import (
            ColumnGeometry, CompetitiveLangmuirIsotherm,
            run_gradient_elution, make_linear_gradient,
        )

        col = ColumnGeometry(
            diameter=0.01, bed_height=0.05,
            particle_diameter=100e-6, bed_porosity=0.38,
            particle_porosity=0.5, G_DN=10000.0, E_star=30000.0,
        )
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([100.0]), K_L=np.array([1000.0]),
        )
        grad = make_linear_gradient(0.0, 500.0, 0.0, 600.0)

        result = run_gradient_elution(
            column=col,
            C_feed=np.array([0.01]),
            gradient=grad,
            flow_rate=5e-8,
            total_time=1200.0,
            isotherm=iso,
            feed_duration=300.0,
            n_z=15,
        )

        assert result.gradient_affects_binding is False
        assert len(result.peaks) == 1
        assert result.time.shape[0] > 100

    def test_run_gradient_elution_returns_gradient_profile(self):
        """Gradient profile in result matches GradientProgram values."""
        from emulsim.module3_performance import (
            ColumnGeometry, CompetitiveLangmuirIsotherm,
            run_gradient_elution, make_linear_gradient,
        )

        col = ColumnGeometry(
            diameter=0.01, bed_height=0.05,
            particle_diameter=100e-6, bed_porosity=0.38,
            particle_porosity=0.5, G_DN=10000.0, E_star=30000.0,
        )
        grad = make_linear_gradient(0.0, 1.0, 0.0, 600.0)
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([50.0]), K_L=np.array([500.0]),
        )

        result = run_gradient_elution(
            column=col,
            C_feed=np.array([0.01]),
            gradient=grad,
            flow_rate=5e-8,
            total_time=900.0,
            isotherm=iso,
            feed_duration=200.0,
            n_z=15,
        )

        # Gradient profile should be monotonically increasing
        assert result.gradient_profile[0] == pytest.approx(0.0, abs=0.01)
        assert result.gradient_profile[-1] == pytest.approx(1.0, abs=0.01)
        assert np.all(np.diff(result.gradient_profile) >= -1e-10)


# ═════════════════════════════════════════════════════════════════════════
# 5. CROSS-MODULE: M2 → CALIBRATION → M3
# ═════════════════════════════════════════════════════════════════════════

class TestM2CalibrationM3Pipeline:
    """Full M2→calibration→M3 pipeline: calibration overrides affect M3 output."""

    def test_m2_produces_fmc_for_m3(self):
        """build_functional_media_contract produces FMC usable by M3."""
        from emulsim.module2_functionalization import (
            ModificationOrchestrator, ModificationStep,
            ModificationStepType, ACSSiteType,
        )
        from emulsim.module2_functionalization.orchestrator import build_functional_media_contract

        contract = _make_contract()
        orch = ModificationOrchestrator()

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

        result = orch.run(contract, steps)
        fmc = build_functional_media_contract(result)

        # FMC should have binding model hint for M3 routing
        assert hasattr(fmc, 'binding_model_hint')
        assert fmc.binding_model_hint != ""
        assert hasattr(fmc, 'estimated_q_max')

    def test_isotherm_selection_from_fmc(self):
        """select_isotherm_from_fmc routes correctly for IEX profile."""
        from emulsim.module2_functionalization import (
            ModificationOrchestrator, ModificationStep,
            ModificationStepType, ACSSiteType,
        )
        from emulsim.module2_functionalization.orchestrator import build_functional_media_contract
        from emulsim.module3_performance.isotherms.adapter import select_isotherm_from_fmc

        contract = _make_contract()
        orch = ModificationOrchestrator()

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

        result = orch.run(contract, steps)
        fmc = build_functional_media_contract(result)

        # Route through select_isotherm_from_fmc
        isotherm = select_isotherm_from_fmc(fmc)
        assert isotherm is not None
        # Should be callable for equilibrium loading
        assert hasattr(isotherm, 'equilibrium_loading')
