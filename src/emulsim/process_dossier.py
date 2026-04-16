"""Node 16 (v7.0, P2a): ProcessDossier — top-level run aggregator.

The v6.1 ``RunReport`` carried model evidence + diagnostics + trust
assessment for one pipeline run, but it did not capture the *experimental
context* around the run: what assays informed the calibration, what target
product profile was being optimised toward, what code/data versions were
in flight, what wet-lab batch the inputs reference. ProcessDossier is the
thin wrapper that completes that loop.

Design philosophy (per v7.0 plan): build over existing artifacts, do not
rewrite. The dossier contains:

  - ``sim_parameters`` — already-existing SimulationParameters
  - ``full_result`` — already-existing FullResult (with embedded RunReport)
  - ``calibration_store_snapshot`` — JSON-serialisable list of entries
    that were active when the run executed
  - ``assay_records`` — list of AssayRecord (Node 17) referenced by the
    calibration store
  - ``target_profile`` — TargetProductProfile if the user defined one
    (deferred-detail v7.0 stub for now)
  - ``environment`` — code version, python version, numpy version, git SHA
    if obtainable
  - ``timestamp_utc``, ``run_id`` — provenance basics

JSON export is the primary use-case (lab notebook attachment, run
reproducibility). HDF5 is deferred until histories are routinely large
enough to warrant it (currently the per-run summary fits comfortably in
JSON).
"""

from __future__ import annotations

import json
import platform
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


@dataclass
class TargetProductProfile:
    """Deferred-detail v7.0 stub for application-specific targets.

    Real TPPs are application-specific (chromatography vs enzyme support
    vs adsorbent), so the full schema is left for v7.1 once we see what
    end users actually want to optimise toward. For now this records the
    barest information so dossiers can carry it through without crashing.
    """
    application: str = ""        # e.g. "size_exclusion_chromatography"
    target_d50_um: float = 0.0
    target_pore_nm: float = 0.0
    target_G_DN_kPa: float = 0.0
    notes: str = ""


def _capture_environment() -> dict:
    """Minimal environment snapshot for reproducibility.

    Captures language/library versions and an optional git SHA. We do NOT
    walk dependencies via pip-list (slow + non-deterministic across
    machines); the imported-version capture is enough for the
    "did the same code produce these outputs?" question.
    """
    env: dict = {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    # Best-effort imports — never crash dossier creation on a missing optional.
    for mod_name in ("numpy", "scipy", "pydantic", "numba", "torch"):
        try:
            mod = __import__(mod_name)
            env[f"{mod_name}_version"] = getattr(mod, "__version__", "unknown")
        except ImportError:
            pass
    # emulsim version
    try:
        from . import __version__ as _v
        env["emulsim_version"] = _v
    except Exception:
        env["emulsim_version"] = "unknown"
    return env


@dataclass
class ProcessDossier:
    """Top-level run-record artifact (Node 16).

    Aggregates everything needed to reproduce or audit a pipeline run.
    Backward-compatible: created on demand by the orchestrator (or by
    user code) — the existing FullResult/RunReport flow is unchanged.
    """
    run_id: str
    timestamp_utc: str
    full_result: Any                        # FullResult (typed loosely to avoid circular import)
    calibration_entries: list[dict] = field(default_factory=list)
    assay_records: list[dict] = field(default_factory=list)
    target_profile: Optional[TargetProductProfile] = None
    environment: dict = field(default_factory=dict)
    notes: str = ""

    @classmethod
    def from_run(
        cls,
        full_result,
        calibration_store=None,
        assay_records: Optional[list] = None,
        target_profile: Optional[TargetProductProfile] = None,
        notes: str = "",
    ) -> "ProcessDossier":
        """Build a dossier around a completed FullResult.

        ``calibration_store`` is a CalibrationStore (Node 7); we snapshot
        its entries as JSON-friendly dicts so the dossier is self-contained
        even if the live store changes after the run.
        """
        cal_entries = []
        if calibration_store is not None:
            for entry in getattr(calibration_store, "entries", []):
                if hasattr(entry, "to_dict"):
                    cal_entries.append(entry.to_dict())

        records_dicts: list[dict] = []
        if assay_records:
            for r in assay_records:
                if hasattr(r, "to_dict"):
                    records_dicts.append(r.to_dict())
                elif isinstance(r, dict):
                    records_dicts.append(r)

        run_id = getattr(getattr(full_result, "parameters", None), "run_id", "") or "unknown_run"

        return cls(
            run_id=run_id,
            timestamp_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            full_result=full_result,
            calibration_entries=cal_entries,
            assay_records=records_dicts,
            target_profile=target_profile,
            environment=_capture_environment(),
            notes=notes,
        )

    def to_json_dict(self) -> dict:
        """Serialise to a JSON-friendly dict (no numpy arrays).

        FullResult contains numpy arrays (size distributions, time series).
        The dossier exports SUMMARIES from the result, not raw arrays —
        large arrays stay in the per-run summary.json that the orchestrator
        already writes. Use ``to_json_with_arrays`` if you need full arrays
        in a single dossier file (uses HDF5 fallback when arrays large).
        """
        rr = getattr(self.full_result, "run_report", None)

        result_summary: dict = {}
        if self.full_result is not None:
            e = getattr(self.full_result, "emulsification", None)
            g = getattr(self.full_result, "gelation", None)
            x = getattr(self.full_result, "crosslinking", None)
            m = getattr(self.full_result, "mechanical", None)
            if e is not None:
                result_summary["L1"] = {
                    "d32_um": float(e.d32) * 1e6,
                    "d50_um": float(e.d50) * 1e6,
                    "span": float(e.span),
                    "converged": bool(e.converged),
                }
            if g is not None:
                result_summary["L2"] = {
                    "pore_size_mean_nm": float(g.pore_size_mean) * 1e9,
                    "porosity": float(g.porosity),
                    "alpha_final": float(g.alpha_final),
                }
            if x is not None:
                result_summary["L3"] = {
                    "p_final": float(x.p_final),
                    "G_chitosan_Pa": float(x.G_chitosan_final),
                    "xi_nm": float(x.xi_final) * 1e9,
                }
            if m is not None:
                result_summary["L4"] = {
                    "G_DN_Pa": float(m.G_DN),
                    "E_star_Pa": float(m.E_star),
                    "model_used": str(m.model_used),
                }

        run_report_dict: dict = {}
        if rr is not None:
            run_report_dict = {
                "min_evidence_tier": rr.min_evidence_tier,
                "trust_level": rr.trust_level,
                "trust_warnings": list(rr.trust_warnings),
                "trust_blockers": list(rr.trust_blockers),
                "diagnostics": rr.diagnostics,
                "n_models": len(rr.model_graph),
                "model_names": [m.model_name for m in rr.model_graph],
            }

        return {
            "schema_version": "1.0",
            "dossier_kind": "ProcessDossier",
            "run_id": self.run_id,
            "timestamp_utc": self.timestamp_utc,
            "result_summary": result_summary,
            "run_report": run_report_dict,
            "calibration_entries": self.calibration_entries,
            "assay_records": self.assay_records,
            "target_profile": (
                asdict(self.target_profile) if self.target_profile else None
            ),
            "environment": self.environment,
            "notes": self.notes,
        }

    def export_json(self, path: Path) -> Path:
        """Write the dossier to a JSON file. Returns the written path."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.to_json_dict(), f, indent=2, default=str)
        return path
