"""Node 20 (v7.0, P1a): calibration fitter scaffold.

Stub fitters that consume ``AssayRecord`` lists from ``data/validation/``
and emit ``CalibrationEntry`` lists ready for ``CalibrationStore.load_json``.

Real fitters (least-squares for empirical models, Bayesian posterior for
calibration-posterior tier) are part of Node 21 — gated on actual wet-lab
data delivery (Study A). The scaffold here lets the data-flow pipeline be
exercised end-to-end with synthetic / placeholder fits so downstream
consumers (RunContext, ProcessDossier) work today.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from ..assay_record import AssayKind, AssayRecord
from .calibration_data import CalibrationEntry

logger = logging.getLogger(__name__)


def load_assay_records(directory: Path) -> list[AssayRecord]:
    """Load every ``*.json`` AssayRecord under a validation subdirectory."""
    directory = Path(directory)
    if not directory.exists():
        return []
    records: list[AssayRecord] = []
    for path in sorted(directory.glob("*.json")):
        try:
            with open(path, encoding="utf-8") as f:
                d = json.load(f)
            records.append(AssayRecord.from_dict(d))
        except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
            logger.warning("Skipping malformed assay JSON %s: %s", path, exc)
    return records


def fit_l1_dsd_to_calibration_entries(
    records: list[AssayRecord],
    profile_key: str = "rotor_stator_legacy",
) -> list[CalibrationEntry]:
    """Stub L1 DSD fitter (Node 20).

    Real implementation (Node 21) will:
      1. Group records by stirrer_type + surfactant_conc + mu_d.
      2. Build (RPM, mean_d32) tuples; minimise residual against
         ``breakage_rate_dispatch + coalescence_rate_dispatch`` over
         (C1, C2, C3, C4, C5) using scipy.optimize.least_squares.
      3. Compute posterior std via Hessian inversion at the optimum and
         emit ``CalibrationEntry(posterior_uncertainty=...)``.

    For v7.0 scaffold this just emits identity entries — value = mean,
    posterior = sample std — so the data-flow pipeline can be tested
    without crashing on an empty fitter.
    """
    out: list[CalibrationEntry] = []
    dsd_records = [r for r in records if r.kind == AssayKind.DROPLET_SIZE_DISTRIBUTION]
    if not dsd_records:
        logger.info("fit_l1_dsd: no DROPLET_SIZE_DISTRIBUTION records found")
        return out
    # Stub: collapse all records to a single d32 reference value.
    all_values = []
    for r in dsd_records:
        all_values.extend(r.values())
    if not all_values:
        return out
    import statistics
    mean_d32 = float(statistics.mean(all_values))
    std_d32 = float(statistics.stdev(all_values)) if len(all_values) >= 2 else 0.0
    out.append(CalibrationEntry(
        profile_key=profile_key,
        parameter_name="d32_reference",
        measured_value=mean_d32,
        units="m",
        confidence="high" if std_d32 / max(mean_d32, 1e-300) < 0.10 else "medium",
        source_reference=f"L1_DSD fit n={len(all_values)} replicates",
        target_module="L1",
        fit_method="stub_mean",
        posterior_uncertainty=std_d32,
    ))
    return out


def write_calibration_json(
    entries: list[CalibrationEntry],
    output_path: Path,
    fit_metadata: Optional[dict] = None,
) -> Path:
    """Write a fits/ JSON loadable by ``CalibrationStore.load_json``."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = [e.to_dict() for e in entries]
    # CalibrationStore.load_json expects a top-level list. Carry fit
    # metadata in a sidecar file alongside the fit JSON.
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    if fit_metadata:
        sidecar = output_path.with_suffix(".meta.json")
        meta = dict(fit_metadata)
        meta["timestamp_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        meta["entry_count"] = len(entries)
        with open(sidecar, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)
    return output_path
