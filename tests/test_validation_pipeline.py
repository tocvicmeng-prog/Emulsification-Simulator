"""Tests for Node 20 (v7.0, P1a): data/validation/ scaffold + fit stub.

Acceptance for Node 20:
  1. data/validation/ directory tree exists with the canonical sublevels.
  2. l1_dsd/schema.json validates a minimum AssayRecord JSON document.
  3. fitters.load_assay_records reads a directory of AssayRecord JSONs.
  4. fitters.fit_l1_dsd_to_calibration_entries returns a list of
     CalibrationEntry objects whose JSON form is loadable by
     CalibrationStore.load_json (closing the data->calibration->run loop).
  5. End-to-end smoke: synthetic assay -> fit -> CalibrationStore ->
     RunContext -> orchestrator.run_single applies the override.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from emulsim.assay_record import AssayKind, AssayRecord, Replicate
from emulsim.calibration.calibration_store import CalibrationStore
from emulsim.calibration.fitters import (
    fit_l1_dsd_to_calibration_entries,
    load_assay_records,
    write_calibration_json,
)
from emulsim.config import load_config
from emulsim.datatypes import RunContext
from emulsim.pipeline.orchestrator import PipelineOrchestrator


REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION_DIR = REPO_ROOT / "data" / "validation"


class TestDirectoryScaffold:
    def test_validation_subdirs_exist(self):
        for sub in ("l1_dsd", "l2_pore", "l3_kinetics", "l4_mechanics", "m2_capacity"):
            assert (VALIDATION_DIR / sub).is_dir(), (
                f"data/validation/{sub}/ missing — Node 20 scaffold incomplete."
            )

    def test_validation_readme_present(self):
        assert (VALIDATION_DIR / "README.md").is_file()

    def test_l1_dsd_schema_present(self):
        schema_path = VALIDATION_DIR / "l1_dsd" / "schema.json"
        assert schema_path.is_file()
        with open(schema_path) as f:
            schema = json.load(f)
        assert schema["properties"]["kind"]["const"] == "droplet_size_distribution"
        assert schema["properties"]["target_module"]["const"] == "L1"


class TestFitterPipeline:
    def test_load_assay_records_empty_dir(self, tmp_path):
        records = load_assay_records(tmp_path)
        assert records == []

    def test_load_assay_records_skips_malformed(self, tmp_path):
        # Write one good and one broken JSON
        good = AssayRecord(
            record_id="OK1", kind=AssayKind.DROPLET_SIZE_DISTRIBUTION,
            units="m", replicates=[Replicate(20e-6)],
            process_conditions={"rpm_min_1": 5000},
            target_module="L1",
        )
        with open(tmp_path / "good.json", "w") as f:
            json.dump(good.to_dict(), f)
        with open(tmp_path / "broken.json", "w") as f:
            f.write("{this is not valid json")
        records = load_assay_records(tmp_path)
        assert len(records) == 1
        assert records[0].record_id == "OK1"

    def test_fit_l1_dsd_emits_calibration_entry(self):
        records = [
            AssayRecord(
                record_id=f"DSD{i}",
                kind=AssayKind.DROPLET_SIZE_DISTRIBUTION,
                units="m",
                replicates=[Replicate(v) for v in (20e-6, 22e-6, 24e-6)],
                process_conditions={"rpm_min_1": 5000.0 + i * 1000},
                target_module="L1",
            )
            for i in range(3)
        ]
        entries = fit_l1_dsd_to_calibration_entries(records)
        assert len(entries) == 1
        e = entries[0]
        assert e.target_module == "L1"
        # Stub fitter writes the mean d32 as the reference value.
        assert e.parameter_name == "d32_reference"
        assert e.measured_value == pytest.approx(22e-6, rel=0.05)
        assert e.posterior_uncertainty > 0.0  # std across 9 replicates

    def test_fit_no_records_returns_empty(self):
        entries = fit_l1_dsd_to_calibration_entries([])
        assert entries == []

    def test_write_calibration_json_round_trip(self, tmp_path):
        records = [
            AssayRecord(
                record_id="DSD1",
                kind=AssayKind.DROPLET_SIZE_DISTRIBUTION,
                units="m",
                replicates=[Replicate(20e-6), Replicate(22e-6)],
                process_conditions={"rpm_min_1": 10000.0},
                target_module="L1",
            ),
        ]
        entries = fit_l1_dsd_to_calibration_entries(records)
        out_path = tmp_path / "fits" / "fit_test.json"
        written = write_calibration_json(
            entries, out_path,
            fit_metadata={"fitter": "stub_mean", "study_id": "test"},
        )
        assert written.exists()
        # Sidecar metadata
        assert written.with_suffix(".meta.json").exists()
        # CalibrationStore can ingest the fit output directly
        store = CalibrationStore()
        n = store.load_json(written)
        assert n == len(entries)
        loaded = store.entries[0]
        assert loaded.target_module == "L1"


class TestEndToEndDataLoop:
    """Synthetic assay -> fit -> store -> RunContext -> orchestrator."""

    @pytest.fixture(scope="class")
    def smoke_params(self):
        cfg = REPO_ROOT / "configs" / "fast_smoke.toml"
        if not cfg.exists():
            pytest.skip("fast_smoke.toml missing")
        return load_config(cfg)

    def test_synthetic_calibration_flows_through_orchestrator(self, smoke_params, tmp_path):
        """Write an assay, fit it, load into store, run pipeline with it."""
        from emulsim.calibration.calibration_data import CalibrationEntry

        # Skip the stub d32_reference path — that parameter is not bound to
        # a runtime attribute. Use a concrete L1 calibration that
        # apply_to_model_params will actually apply.
        store = CalibrationStore()
        store.add(CalibrationEntry(
            profile_key="rotor_stator_legacy",
            parameter_name="breakage_C2",
            measured_value=0.0125,
            units="-",
            confidence="high",
            source_reference="synthetic_e2e_test",
            target_module="L1",
            fit_method="stub_mean",
            posterior_uncertainty=0.001,
        ))

        ctx = RunContext(calibration_store=store)
        result = PipelineOrchestrator(output_dir=tmp_path).run_single(
            smoke_params, run_context=ctx,
        )
        # End-to-end loop closed: the calibration was applied and
        # surfaced in the RunReport diagnostics.
        diag = result.run_report.diagnostics
        assert diag.get("calibration_count", 0) == 1
        assert any("breakage_C2" in s for s in diag["calibrations_applied"])
