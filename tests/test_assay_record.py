"""Tests for Node 17 (v7.0, P2b): AssayRecord public data model.

Acceptance for Node 17:
  1. AssayRecord round-trips through JSON (to_dict + from_dict identity).
  2. Replicate statistics (mean, std, cv, n_replicates) work over flagged
     and unflagged replicates.
  3. Outlier-flagged replicates are excluded from statistics.
  4. Schema version is captured in the serialised form.
  5. AssayKind enum covers the canonical wet-lab measurements named in
     docs/35 §9 and docs/10 §5 (the calibration campaign).
"""

from __future__ import annotations

import json
import math

import pytest

from emulsim.assay_record import (
    ASSAY_SCHEMA_VERSION,
    AssayKind,
    AssayRecord,
    Replicate,
)


# ─── Schema coverage ──────────────────────────────────────────────────────


class TestKindCoverage:
    def test_all_wet_lab_kinds_present(self):
        names = {k.value for k in AssayKind}
        for required in (
            "droplet_size_distribution",
            "interfacial_tension",
            "dispersed_viscosity",
            "pore_size",
            "porosity",
            "gelation_onset",
            "crosslink_conversion",
            "bulk_modulus",
            "static_binding_capacity",
            "dynamic_binding_capacity",
            "ligand_density",
            "activity_retention",
        ):
            assert required in names, f"AssayKind missing canonical kind {required!r}"


# ─── Replicate statistics ─────────────────────────────────────────────────


class TestReplicateStats:
    def _three_rep(self, values, flags=None):
        flags = flags or [""] * len(values)
        return AssayRecord(
            record_id="t",
            kind=AssayKind.DROPLET_SIZE_DISTRIBUTION,
            units="m",
            replicates=[Replicate(v, flag=f) for v, f in zip(values, flags)],
            process_conditions={"rpm": 5000.0},
        )

    def test_mean_std_cv(self):
        rec = self._three_rep([10.0, 12.0, 14.0])
        assert rec.mean() == pytest.approx(12.0, abs=1e-12)
        assert rec.std() == pytest.approx(2.0, abs=1e-9)  # sample std of 10/12/14
        assert rec.cv() == pytest.approx(2.0 / 12.0, abs=1e-9)
        assert rec.n_replicates() == 3

    def test_outlier_excluded(self):
        rec = self._three_rep([10.0, 12.0, 14.0, 1000.0],
                               flags=["", "", "", "outlier"])
        assert rec.mean() == pytest.approx(12.0)
        assert rec.n_replicates() == 3

    def test_single_replicate_zero_std(self):
        rec = self._three_rep([42.0])
        assert rec.std() == 0.0
        assert rec.n_replicates() == 1

    def test_zero_mean_yields_nan_cv(self):
        rec = self._three_rep([0.0, 0.0, 0.0])
        assert math.isnan(rec.cv())


# ─── JSON round-trip ──────────────────────────────────────────────────────


class TestJSONRoundTrip:
    def test_roundtrip_preserves_fields(self):
        original = AssayRecord(
            record_id="LAB-2026-04-17-001",
            kind=AssayKind.INTERFACIAL_TENSION,
            units="N/m",
            replicates=[
                Replicate(value=4.8e-3, std=2e-4),
                Replicate(value=5.0e-3, std=2e-4),
                Replicate(value=5.2e-3, std=2e-4),
            ],
            process_conditions={"surfactant_conc_kg_m3": 20.0, "T_K": 363.15},
            sample_id="agarose_lot_A",
            instrument="KrussDSA100_S/N12345",
            operator="TC",
            notebook_ref="LN-2026-04-17-p42",
            target_module="L1",
            notes="triplicate pendant drop after 60s equilibration",
        )
        d = original.to_dict()
        json_str = json.dumps(d)  # must be JSON-serialisable
        restored_dict = json.loads(json_str)
        restored = AssayRecord.from_dict(restored_dict)

        assert restored.record_id == original.record_id
        assert restored.kind == original.kind
        assert restored.units == original.units
        assert len(restored.replicates) == len(original.replicates)
        assert restored.replicates[1].value == pytest.approx(5.0e-3)
        assert restored.process_conditions["surfactant_conc_kg_m3"] == 20.0
        assert restored.target_module == "L1"

    def test_schema_version_in_dict(self):
        rec = AssayRecord(
            record_id="x", kind=AssayKind.PORE_SIZE, units="m",
        )
        d = rec.to_dict()
        assert d["schema_version"] == ASSAY_SCHEMA_VERSION
        assert d["kind"] == "pore_size"

    def test_unknown_schema_version_does_not_crash(self):
        """v7.1+ migration hook: tolerate unknown schema_version field."""
        d = {
            "schema_version": "9.99",
            "record_id": "x",
            "kind": "pore_size",
            "units": "m",
            "replicates": [{"value": 80e-9, "std": 0.0, "flag": "", "notes": ""}],
            "process_conditions": {},
            "sample_id": "", "instrument": "", "operator": "",
            "notebook_ref": "", "timestamp_utc": "", "target_module": "",
            "notes": "",
        }
        rec = AssayRecord.from_dict(d)
        assert rec.kind == AssayKind.PORE_SIZE


# ─── Integration with ProcessDossier ──────────────────────────────────────


class TestDossierIntegration:
    def test_dossier_carries_assay_records(self, tmp_path):
        """ProcessDossier (Node 16) should accept AssayRecord lists."""
        from emulsim.config import load_config
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.process_dossier import ProcessDossier
        from pathlib import Path

        cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
        if not cfg.exists():
            pytest.skip("fast_smoke.toml missing")
        params = load_config(cfg)
        result = PipelineOrchestrator(output_dir=tmp_path).run_single(params)

        records = [
            AssayRecord(
                record_id="DSD-001",
                kind=AssayKind.DROPLET_SIZE_DISTRIBUTION,
                units="m",
                replicates=[Replicate(22e-6), Replicate(23e-6), Replicate(21e-6)],
                process_conditions={"rpm": 10000.0},
            ),
        ]
        dossier = ProcessDossier.from_run(result, assay_records=records)
        assert len(dossier.assay_records) == 1
        assert dossier.assay_records[0]["record_id"] == "DSD-001"
        # Round-trip via JSON
        path = dossier.export_json(tmp_path / "dossier_with_assays.json")
        with open(path) as f:
            loaded = json.load(f)
        assert loaded["assay_records"][0]["kind"] == "droplet_size_distribution"
