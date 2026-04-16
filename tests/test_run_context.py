"""Tests for Node 7 (v6.1): RunContext-driven calibration injection.

Acceptance criteria for Node 7:
  1. Backward compat — orchestrator.run_single without run_context produces
     the same outputs as v6.0 (verified indirectly by the existing
     test_smoke + test_v60_integration suites passing unchanged).
  2. With run_context.calibration_store populated, calibration entries
     targeting L1 KernelConfig fields are applied before solving and the
     change is reflected in the actual L1 result (different d32 vs no
     calibration).
  3. RunReport.diagnostics records the override list so downstream
     consumers (UI, JSON export, optimizer) can see what was calibrated.
  4. has_calibration_for("L1") gates the injection — runs against a store
     with no L1 entries reproduce the no-context baseline.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from emulsim.calibration.calibration_data import CalibrationEntry
from emulsim.calibration.calibration_store import CalibrationStore
from emulsim.config import load_config
from emulsim.datatypes import RunContext
from emulsim.pipeline.orchestrator import PipelineOrchestrator


@pytest.fixture(scope="module")
def smoke_params():
    """SimulationParameters from configs/fast_smoke.toml — runs in <1 s."""
    repo_root = Path(__file__).resolve().parents[1]
    cfg = repo_root / "configs" / "fast_smoke.toml"
    if not cfg.exists():
        pytest.skip(f"fast_smoke.toml not found at {cfg}")
    return load_config(cfg)


def test_baseline_without_run_context(smoke_params, tmp_path):
    """run_single with run_context=None reproduces v6.0 behaviour."""
    orch = PipelineOrchestrator(output_dir=tmp_path)
    result = orch.run_single(smoke_params)
    assert result.run_report is not None
    # Baseline runs do not record any calibrations
    diag = result.run_report.diagnostics
    assert diag.get("calibration_count", 0) == 0
    assert "calibrations_applied" not in diag


def test_calibration_store_overrides_l1_kernel(smoke_params, tmp_path):
    """An L1 calibration entry mutates the kernel constant before solving."""
    orch = PipelineOrchestrator(output_dir=tmp_path)

    # Baseline: no calibration
    baseline = orch.run_single(smoke_params)
    baseline_d32 = baseline.emulsification.d32

    # Build a CalibrationStore with an L1 entry that doubles breakage_C1.
    # The dispatch routing fix (Node 3) means this constant now actually
    # reaches the legacy solver.
    store = CalibrationStore()
    store.add(CalibrationEntry(
        profile_key="rotor_stator_legacy",
        parameter_name="breakage_C1",
        measured_value=2.0,            # default is 0.986
        units="-",
        confidence="medium",
        source_reference="synthetic_test",
        target_module="L1",
        fit_method="manual",
    ))
    ctx = RunContext(calibration_store=store, run_id="cal_test")

    calibrated = orch.run_single(smoke_params, run_context=ctx)
    cal_d32 = calibrated.emulsification.d32

    # The override list is in RunReport.diagnostics
    diag = calibrated.run_report.diagnostics
    assert diag["calibration_count"] == 1
    assert any("breakage_C1" in s for s in diag["calibrations_applied"])

    # Doubling the breakage prefactor should change d32. Allow either direction
    # (depends on regime); assert only that it is actually different.
    assert cal_d32 != baseline_d32, (
        f"L1 calibration had no effect on d32 (baseline={baseline_d32:.4e}, "
        f"calibrated={cal_d32:.4e}). Either dispatch routing regressed (Node 3) "
        f"or apply_to_model_params is not finding the attribute."
    )


def test_empty_store_matches_baseline(smoke_params, tmp_path):
    """A RunContext with an empty CalibrationStore reproduces baseline d32."""
    orch = PipelineOrchestrator(output_dir=tmp_path)
    baseline = orch.run_single(smoke_params)

    ctx = RunContext(calibration_store=CalibrationStore())
    cal = orch.run_single(smoke_params, run_context=ctx)

    assert cal.emulsification.d32 == baseline.emulsification.d32
    assert cal.run_report.diagnostics.get("calibration_count", 0) == 0


def test_unrelated_module_does_not_trigger_l1_branch(smoke_params, tmp_path):
    """A store with only an M2 entry does not touch L1 kernels or props."""
    orch = PipelineOrchestrator(output_dir=tmp_path)
    baseline = orch.run_single(smoke_params)

    store = CalibrationStore()
    store.add(CalibrationEntry(
        profile_key="some_reagent",
        parameter_name="estimated_q_max",
        measured_value=42.0,
        units="mol/m3",
        confidence="high",
        source_reference="external_assay",
        target_module="M2",
    ))
    ctx = RunContext(calibration_store=store)
    cal = orch.run_single(smoke_params, run_context=ctx)

    # No L1/L2/L3/L4 calibration was applied during the M1 pipeline run.
    assert cal.emulsification.d32 == baseline.emulsification.d32
    assert cal.run_report.diagnostics.get("calibration_count", 0) == 0
