"""Fast end-to-end smoke gate.

Runs the full M1 pipeline (L1 -> L2 -> L3 -> L4) using configs/fast_smoke.toml
and asserts:

  1. The run finishes in well under the budget (target <= 30 s, observed ~0.2 s).
  2. Each level produces finite, in-range outputs (sanity bounds, NOT scientific
     correctness — that is the job of the level-specific test modules).
  3. A RunReport is attached and carries a defined evidence tier.

These are the gates that should run on every commit before more expensive
mechanistic, optimisation, or uncertainty tests. Mark with @pytest.mark.smoke
so they can be selected with ``pytest -m smoke``.

The expected baselines below are hard-coded for the *current* fast_smoke.toml
defaults; if defaults are tuned, update both this file and configs/fast_smoke.toml
in the same commit.
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import pytest


# Per-test budget. fast_smoke.toml ran in ~0.13 s during initial probe; 30 s is
# a generous CI budget that catches catastrophic regressions (e.g., accidental
# n_grid=128 reintroduction) without false-flagging on slow CI runners.
SMOKE_BUDGET_S = 30.0

# Hard sanity bounds. These are NOT calibrated targets — they are guards against
# obviously-broken outputs (e.g., d32 == 0, NaN propagation, negative porosity).
EXPECTED_RANGES = {
    "d32_um":   (0.5, 500.0),    # microns; 500 um upper guards against PBE blow-up
    "pore_nm":  (1.0, 5_000.0),  # nm; 5 um upper guards against gel collapse
    "p_final":  (0.0, 1.0),      # conversion fraction
    "G_DN_Pa":  (1.0, 1.0e7),    # Pa; 10 MPa upper guards against runaway DN model
}


@pytest.fixture(scope="module")
def smoke_config_path() -> Path:
    """Path to the fast smoke TOML, anchored at the repo root."""
    repo_root = Path(__file__).resolve().parents[1]
    cfg = repo_root / "configs" / "fast_smoke.toml"
    if not cfg.exists():
        pytest.skip(f"fast_smoke.toml not found at {cfg}")
    return cfg


@pytest.mark.smoke
def test_pipeline_runs_under_budget(smoke_config_path: Path, tmp_path: Path):
    """End-to-end M1 pipeline finishes under SMOKE_BUDGET_S with finite outputs."""
    from emulsim.config import load_config
    from emulsim.pipeline.orchestrator import PipelineOrchestrator

    params = load_config(smoke_config_path)
    orch = PipelineOrchestrator(output_dir=tmp_path)

    t0 = time.perf_counter()
    result = orch.run_single(params)
    elapsed = time.perf_counter() - t0

    # Budget — catastrophic-regression guard
    assert elapsed < SMOKE_BUDGET_S, (
        f"Smoke pipeline took {elapsed:.2f}s, exceeding {SMOKE_BUDGET_S}s budget. "
        "Did fast_smoke.toml dimensions get widened, or did a solver regress?"
    )

    # Sanity guards on each level's headline output. Use isfinite to catch NaN.
    metrics = {
        "d32_um":  result.emulsification.d32 * 1e6,
        "pore_nm": result.gelation.pore_size_mean * 1e9,
        "p_final": result.crosslinking.p_final,
        "G_DN_Pa": result.mechanical.G_DN,
    }
    for name, value in metrics.items():
        lo, hi = EXPECTED_RANGES[name]
        assert math.isfinite(value), f"{name} is non-finite: {value!r}"
        assert lo <= value <= hi, (
            f"{name}={value:.4g} outside sanity range [{lo}, {hi}]. "
            "This is a smoke guard, not a calibration check — investigate before tightening."
        )


@pytest.mark.smoke
def test_run_report_attached(smoke_config_path: Path, tmp_path: Path):
    """Smoke run must populate FullResult.run_report with a defined evidence tier."""
    from emulsim.config import load_config
    from emulsim.pipeline.orchestrator import PipelineOrchestrator

    params = load_config(smoke_config_path)
    orch = PipelineOrchestrator(output_dir=tmp_path)
    result = orch.run_single(params)

    # RunReport plumbing was added in v6.1; if it disappears, M2/M3 evidence
    # inheritance (consensus plan Sprint 1) silently breaks.
    assert getattr(result, "run_report", None) is not None, (
        "FullResult.run_report missing. The orchestrator must assemble a RunReport "
        "with model manifests and a min_evidence_tier per the consensus plan."
    )
    rr = result.run_report
    assert rr.min_evidence_tier, "run_report.min_evidence_tier is empty"
    assert isinstance(rr.trust_warnings, list)
    assert isinstance(rr.trust_blockers, list)
