"""Tests for Node 10 (v6.1, F11): targeted unit/range assertions at boundaries.

Acceptance for Node 10:
  1. M1ExportContract.validate_units catches obvious unit confusion (µm vs m,
     mM vs mol/m^3, kPa vs Pa, etc.) but accepts every valid run.
  2. FunctionalMediaContract.validate_units does the same at the M2->M3 boundary.
  3. The orchestrators surface violations as warnings without crashing —
     downstream consumers see them in trust_warnings / warnings lists.
"""

from __future__ import annotations

import pytest

from emulsim.datatypes import M1ExportContract
from emulsim.module2_functionalization.orchestrator import FunctionalMediaContract


def _good_m1_contract(**overrides) -> M1ExportContract:
    base = dict(
        bead_radius=50e-6,
        bead_d32=90e-6,
        bead_d50=100e-6,
        pore_size_mean=80e-9,
        pore_size_std=20e-9,
        porosity=0.7,
        l2_model_tier="empirical_uncalibrated",
        mesh_size_xi=10e-9,
        p_final=0.3,
        primary_crosslinker="genipin",
        nh2_bulk_concentration=80.0,
        oh_bulk_concentration=300.0,
        G_DN=8000.0,
        E_star=24000.0,
        model_used="phenomenological",
        c_agarose=42.0,
        c_chitosan=18.0,
        DDA=0.9,
        trust_level="CAUTION",
    )
    base.update(overrides)
    return M1ExportContract(**base)


class TestM1ContractUnitChecks:
    def test_realistic_contract_passes(self):
        contract = _good_m1_contract()
        assert contract.validate_units() == []

    def test_bead_radius_in_microns_is_caught(self):
        # 50 instead of 50e-6: someone passed "microns" through "metres" slot
        contract = _good_m1_contract(bead_radius=50.0, bead_d50=100.0)
        violations = contract.validate_units()
        assert any("bead_radius" in v for v in violations)

    def test_pore_in_nanometres_value_passed_as_metres(self):
        # 80 instead of 80e-9
        contract = _good_m1_contract(pore_size_mean=80.0)
        violations = contract.validate_units()
        assert any("pore_size_mean" in v for v in violations)

    def test_porosity_above_one_caught(self):
        contract = _good_m1_contract(porosity=1.4)  # bad
        violations = contract.validate_units()
        assert any("porosity" in v for v in violations)

    def test_modulus_in_gpa_caught(self):
        # 1e10 = 10 GPa — outside the 1 Pa..1 GPa range; flagged.
        contract = _good_m1_contract(G_DN=1e10)
        violations = contract.validate_units()
        assert any("G_DN" in v for v in violations)

    def test_amine_concentration_in_kg_per_m3_caught(self):
        # 18 kg/m^3 chitosan -> ~80 mol/m^3 NH2. Passing 18000 mol/m^3
        # almost certainly means kg/m^3 was sent through the mol/m^3 slot.
        contract = _good_m1_contract(nh2_bulk_concentration=1.5e4)
        violations = contract.validate_units()
        assert any("nh2_bulk_concentration" in v for v in violations)


class TestFMCUnitChecks:
    def test_default_fmc_passes(self):
        fmc = FunctionalMediaContract()
        assert fmc.validate_units() == []

    def test_ligand_density_in_mol_per_m3_caught(self):
        # ligand densities are mol/m^2; values > 1 strongly suggest mol/m^3.
        fmc = FunctionalMediaContract(functional_ligand_density=42.0)
        violations = fmc.validate_units()
        assert any("functional_ligand_density" in v for v in violations)

    def test_inverted_qmax_bounds_caught(self):
        fmc = FunctionalMediaContract(
            estimated_q_max=100.0,
            q_max_lower=80.0, q_max_upper=20.0,
        )
        violations = fmc.validate_units()
        assert any("q_max_lower" in v for v in violations)

    def test_negative_qmax_caught(self):
        fmc = FunctionalMediaContract(estimated_q_max=-1.0)
        violations = fmc.validate_units()
        assert any("estimated_q_max" in v for v in violations)


class TestOrchestratorBoundaryWarnings:
    """Whole-pipeline smoke: violations surface as warnings, not crashes."""

    def test_smoke_run_validates_clean(self, tmp_path):
        from emulsim.config import load_config
        from emulsim.pipeline.orchestrator import (
            PipelineOrchestrator, export_for_module2,
        )
        from pathlib import Path
        cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
        if not cfg.exists():
            pytest.skip("fast_smoke.toml missing")
        params = load_config(cfg)
        orch = PipelineOrchestrator(output_dir=tmp_path)
        result = orch.run_single(params)
        # The fast-smoke run must produce a contract that passes its own checks.
        contract = export_for_module2(result, _stub_trust(), props=None)
        assert contract.validate_units() == [], (
            "Fast-smoke run produced an M1ExportContract that fails unit checks; "
            "either the smoke baseline drifted or the validator is too strict."
        )


def _stub_trust():
    """Minimal TrustAssessment-shaped stub for export_for_module2."""
    class _T:
        level = "CAUTION"
        warnings: list = []
    return _T()
