"""Tests for Node 24 (v7.0.1, audit N4): batch / dossier / ingest CLI subcommands.

Acceptance for Node 24:
  1. ``python -m emulsim --help`` lists batch, dossier, ingest as commands.
  2. Each command's ``--help`` produces non-empty output without crash.
  3. ``batch`` runs the fast_smoke config and prints quantile-resolved
     bead radii + mass fractions.
  4. ``dossier`` runs the fast_smoke config and writes a valid
     ProcessDossier JSON to the requested path.
  5. ``ingest L1`` reads AssayRecord JSONs from a directory and writes a
     CalibrationStore-loadable fit JSON.
  6. ``--quantiles`` parser rejects out-of-range / malformed input with a
     clear error.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SMOKE_CFG = REPO_ROOT / "configs" / "fast_smoke.toml"


def _run_cli(*args, cwd=None):
    """Invoke `python -m emulsim ...` and return CompletedProcess."""
    return subprocess.run(
        [sys.executable, "-m", "emulsim", *args],
        cwd=cwd or REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )


class TestCliRegistration:
    def test_top_level_help_lists_new_commands(self):
        cp = _run_cli("--help")
        assert cp.returncode == 0
        for cmd in ("batch", "dossier", "ingest"):
            assert cmd in cp.stdout, f"CLI --help missing {cmd!r}: {cp.stdout!r}"

    def test_each_subcommand_help_works(self):
        for cmd in ("batch", "dossier", "ingest"):
            cp = _run_cli(cmd, "--help")
            assert cp.returncode == 0, f"`emulsim {cmd} --help` failed: {cp.stderr}"
            assert "usage" in cp.stdout.lower()


class TestBatchCommand:
    def test_batch_runs_fast_smoke(self, tmp_path):
        if not SMOKE_CFG.exists():
            pytest.skip("fast_smoke.toml missing")
        cp = _run_cli(
            "batch", str(SMOKE_CFG),
            "--quantiles", "0.25,0.50,0.75",
            "--output", str(tmp_path),
            "--quiet",
        )
        assert cp.returncode == 0, cp.stderr
        assert "Batch Variability Results" in cp.stdout
        assert "Per-quantile representative bead radii" in cp.stdout
        # Quantile lines must be present
        for q in ("0.250", "0.500", "0.750"):
            assert f"q={q}" in cp.stdout

    def test_batch_rejects_bad_quantiles(self, tmp_path):
        if not SMOKE_CFG.exists():
            pytest.skip("fast_smoke.toml missing")
        cp = _run_cli(
            "batch", str(SMOKE_CFG),
            "--quantiles", "0.25,1.5,0.75",  # 1.5 is out of (0,1)
            "--output", str(tmp_path),
            "--quiet",
        )
        assert cp.returncode != 0
        assert "0,1" in cp.stderr or "must" in cp.stderr.lower()


class TestDossierCommand:
    def test_dossier_writes_valid_json(self, tmp_path):
        if not SMOKE_CFG.exists():
            pytest.skip("fast_smoke.toml missing")
        out = tmp_path / "dossier.json"
        cp = _run_cli(
            "dossier", str(SMOKE_CFG),
            "--output", str(out),
            "--notes", "v7.0.1 CLI smoke",
            "--quiet",
        )
        assert cp.returncode == 0, cp.stderr
        assert out.exists()
        with open(out) as f:
            d = json.load(f)
        # Schema sanity
        assert d["dossier_kind"] == "ProcessDossier"
        assert d["schema_version"] == "1.0"
        assert d["notes"] == "v7.0.1 CLI smoke"
        assert "L1" in d["result_summary"]
        assert d["run_report"]["min_evidence_tier"] != ""
        # CLI summary printed the dossier path
        assert "dossier written to" in cp.stdout


class TestIngestCommand:
    def test_ingest_l1_round_trip(self, tmp_path):
        from emulsim.assay_record import AssayKind, AssayRecord, Replicate

        # Build a synthetic AssayRecord JSON
        assay_dir = tmp_path / "assays"
        assay_dir.mkdir()
        rec = AssayRecord(
            record_id="CLI-INGEST-001",
            kind=AssayKind.DROPLET_SIZE_DISTRIBUTION,
            units="m",
            replicates=[Replicate(20e-6), Replicate(22e-6), Replicate(24e-6)],
            process_conditions={"rpm_min_1": 5000.0},
            target_module="L1",
        )
        with open(assay_dir / "dsd_001.json", "w") as f:
            json.dump(rec.to_dict(), f)

        out = tmp_path / "fit.json"
        cp = _run_cli(
            "ingest", "L1",
            "--assay-dir", str(assay_dir),
            "--output", str(out),
        )
        assert cp.returncode == 0, cp.stderr
        assert out.exists()
        # Output must be loadable by CalibrationStore
        from emulsim.calibration.calibration_store import CalibrationStore
        store = CalibrationStore()
        n = store.load_json(out)
        assert n >= 1
        assert store.entries[0].target_module == "L1"
        # Sidecar metadata file is written
        assert out.with_suffix(".meta.json").exists()

    def test_ingest_missing_dir_clear_error(self, tmp_path):
        cp = _run_cli(
            "ingest", "L1",
            "--assay-dir", str(tmp_path / "does_not_exist"),
        )
        assert cp.returncode != 0
        assert "does not exist" in cp.stderr

    def test_ingest_empty_dir_warns(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        out = tmp_path / "fit.json"
        cp = _run_cli(
            "ingest", "L1",
            "--assay-dir", str(empty),
            "--output", str(out),
        )
        # Empty dir is a warning, not an error
        assert cp.returncode == 0, cp.stderr
        assert "no AssayRecord JSONs" in cp.stdout
        assert not out.exists()


class TestUncertaintyCommand:
    """Audit N5 (Node 25): uncertainty CLI surfaces --n-jobs and --engine."""

    def test_unified_engine_default(self, tmp_path):
        if not SMOKE_CFG.exists():
            pytest.skip("fast_smoke.toml missing")
        cp = _run_cli(
            "uncertainty", str(SMOKE_CFG),
            "--n-samples", "3", "--seed", "42",
        )
        assert cp.returncode == 0, cp.stderr
        # Unified output schema markers
        assert "Unified UQ" in cp.stdout
        assert "kinds=[material_property]" in cp.stdout
        assert "engine=unified" in cp.stdout

    def test_legacy_engine_flag_routes_through_unified(self, tmp_path):
        """Node 30: --engine legacy runs the merged engine with no
        posterior injection (scripts get the same MaterialProperties
        perturbations they got in v7.0.x, but in the unified output
        schema). The legacy UncertaintyResult header is gone."""
        if not SMOKE_CFG.exists():
            pytest.skip("fast_smoke.toml missing")
        cp = _run_cli(
            "uncertainty", str(SMOKE_CFG),
            "--n-samples", "3", "--seed", "42",
            "--engine", "legacy",
        )
        assert cp.returncode == 0, cp.stderr
        assert "engine=legacy" in cp.stdout
        assert "Unified UQ" in cp.stdout
        # Without a calibration store, no posterior kinds are sampled.
        assert "kinds=[material_property]" in cp.stdout

    def test_n_jobs_flag_accepted(self, tmp_path):
        """--n-jobs > 1 must not crash (joblib clamps; loky safe)."""
        if not SMOKE_CFG.exists():
            pytest.skip("fast_smoke.toml missing")
        cp = _run_cli(
            "uncertainty", str(SMOKE_CFG),
            "--n-samples", "2", "--seed", "7",
            "--n-jobs", "2",
            "--engine", "legacy",
        )
        assert cp.returncode == 0, cp.stderr
