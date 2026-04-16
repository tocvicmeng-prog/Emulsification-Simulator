"""Tests for Node 16 (v7.0, P2a): ProcessDossier scaffold.

Acceptance for Node 16:
  1. ProcessDossier.from_run can build a dossier from any FullResult
     without requiring optional inputs (calibration_store, assays, target).
  2. Calibration entries from a CalibrationStore are snapshotted as
     JSON-friendly dicts (decoupled from the live store).
  3. JSON export round-trips: writing to disk and reading back gives the
     same dict.
  4. Environment capture records emulsim version + python version + at
     least numpy.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from emulsim.calibration.calibration_data import CalibrationEntry
from emulsim.calibration.calibration_store import CalibrationStore
from emulsim.config import load_config
from emulsim.pipeline.orchestrator import PipelineOrchestrator
from emulsim.process_dossier import ProcessDossier, TargetProductProfile


@pytest.fixture(scope="module")
def smoke_result(tmp_path_factory):
    cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
    if not cfg.exists():
        pytest.skip("fast_smoke.toml missing")
    params = load_config(cfg)
    out = tmp_path_factory.mktemp("dossier_smoke")
    return PipelineOrchestrator(output_dir=out).run_single(params)


def test_from_run_minimal(smoke_result):
    """ProcessDossier builds with only a FullResult."""
    dossier = ProcessDossier.from_run(smoke_result)
    assert dossier.run_id != ""
    assert dossier.timestamp_utc != ""
    assert dossier.calibration_entries == []
    assert dossier.assay_records == []
    assert dossier.target_profile is None


def test_calibration_snapshot(smoke_result):
    """Calibration entries are deep-copied as dicts at dossier creation."""
    store = CalibrationStore()
    store.add(CalibrationEntry(
        profile_key="rotor_stator_legacy",
        parameter_name="breakage_C1",
        measured_value=1.2,
        units="-",
        confidence="medium",
        source_reference="synthetic",
        target_module="L1",
    ))
    dossier = ProcessDossier.from_run(smoke_result, calibration_store=store)
    assert len(dossier.calibration_entries) == 1
    e = dossier.calibration_entries[0]
    assert e["parameter_name"] == "breakage_C1"
    assert e["measured_value"] == 1.2

    # Mutating the live store after dossier creation must NOT affect dossier.
    store.add(CalibrationEntry(
        profile_key="x", parameter_name="y", measured_value=99.0,
        units="-", confidence="low", source_reference="z", target_module="L1",
    ))
    assert len(dossier.calibration_entries) == 1


def test_target_profile_round_trip(smoke_result, tmp_path):
    """TargetProductProfile survives JSON export."""
    tpp = TargetProductProfile(
        application="size_exclusion_chromatography",
        target_d50_um=80.0,
        target_pore_nm=100.0,
        target_G_DN_kPa=20.0,
    )
    dossier = ProcessDossier.from_run(smoke_result, target_profile=tpp)
    path = dossier.export_json(tmp_path / "dossier.json")
    assert path.exists()
    with open(path) as f:
        loaded = json.load(f)
    assert loaded["target_profile"]["application"] == "size_exclusion_chromatography"
    assert loaded["target_profile"]["target_d50_um"] == 80.0


def test_environment_captured(smoke_result):
    """Env snapshot includes versions for the must-have stack."""
    dossier = ProcessDossier.from_run(smoke_result)
    env = dossier.environment
    assert "python_version" in env
    assert "emulsim_version" in env
    assert "numpy_version" in env  # required dep
    assert "timestamp_utc" in env


def test_json_export_includes_run_report(smoke_result, tmp_path):
    """Dossier JSON carries the v6.1 RunReport summary inline."""
    dossier = ProcessDossier.from_run(smoke_result)
    path = dossier.export_json(tmp_path / "dossier_with_rr.json")
    with open(path) as f:
        d = json.load(f)
    assert "run_report" in d
    assert d["run_report"]["min_evidence_tier"] != ""
    assert "model_names" in d["run_report"]
    # Result summary carries L1-L4 headline numbers in JSON-friendly units
    assert "L1" in d["result_summary"]
    assert "d32_um" in d["result_summary"]["L1"]
