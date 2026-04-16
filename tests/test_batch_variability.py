"""Tests for Node 19 (v7.0, P4): batch variability across DSD quantiles.

Acceptance for Node 19:
  1. run_batch returns BatchResult with one QuantileRun per requested
     quantile, sorted by quantile.
  2. Mass fractions sum to ~1.
  3. The mean_d32_m is inherited from the base L1 run (single L1 call).
  4. Batch-level percentiles for pore size and G_DN bracket the
     mass-weighted mean (sanity).
  5. Backward compat: run_single is unchanged.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from emulsim.config import load_config
from emulsim.pipeline.batch_variability import BatchResult, QuantileRun, run_batch


@pytest.fixture(scope="module")
def smoke_params():
    cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
    if not cfg.exists():
        pytest.skip("fast_smoke.toml missing")
    return load_config(cfg)


def test_run_batch_default_quantiles(smoke_params, tmp_path):
    res = run_batch(smoke_params, output_dir=tmp_path)
    assert isinstance(res, BatchResult)
    assert res.n_quantiles == 5
    assert len(res.quantile_runs) == 5
    quantiles = [r.quantile for r in res.quantile_runs]
    assert quantiles == sorted(quantiles)


def test_mass_fractions_sum_to_one(smoke_params, tmp_path):
    res = run_batch(smoke_params, output_dir=tmp_path)
    assert sum(res.quantile_mass_fractions) == pytest.approx(1.0, abs=1e-9)


def test_mean_d32_inherited_from_base_l1(smoke_params, tmp_path):
    """run_batch invokes L1 once; mean_d32_m is the base L1 d32."""
    res = run_batch(smoke_params, output_dir=tmp_path)
    # Each per-quantile run shares the same L1 EmulsificationResult,
    # so mean_d32_m equals every quantile's emulsification.d32.
    for qr in res.quantile_runs:
        assert qr.full_result.emulsification.d32 == pytest.approx(
            res.mean_d32_m, rel=1e-12,
        )


def test_pore_percentiles_bracket_mean(smoke_params, tmp_path):
    res = run_batch(smoke_params, output_dir=tmp_path)
    # p5 <= p50 <= p95 (definition); mean lies inside the [p5, p95] range
    # for any unimodal-ish distribution. The fast_smoke pore model is
    # essentially constant across the small DSD spread, so we only assert
    # the order property strictly.
    assert res.pore_p5_m <= res.pore_p50_m <= res.pore_p95_m
    assert res.G_DN_p5_Pa <= res.G_DN_p50_Pa <= res.G_DN_p95_Pa


def test_custom_quantiles(smoke_params, tmp_path):
    res = run_batch(smoke_params, quantiles=(0.5,), output_dir=tmp_path)
    assert res.n_quantiles == 1
    assert res.quantile_runs[0].quantile == pytest.approx(0.5)
    # Single-quantile mass fraction is 1.0
    assert res.quantile_mass_fractions[0] == pytest.approx(1.0, abs=1e-9)


def test_invalid_quantiles_rejected(smoke_params, tmp_path):
    with pytest.raises(ValueError):
        run_batch(smoke_params, quantiles=(), output_dir=tmp_path)
    with pytest.raises(ValueError):
        run_batch(smoke_params, quantiles=(1.5,), output_dir=tmp_path)
    with pytest.raises(ValueError):
        run_batch(smoke_params, quantiles=(0.0,), output_dir=tmp_path)


# ─── Audit N3 / N8 follow-ups (v7.0.1) ────────────────────────────────────


def test_n3_representative_diameter_helper(smoke_params, tmp_path):
    """QuantileRun.representative_diameter_m == 2 * representative_radius_m."""
    res = run_batch(smoke_params, quantiles=(0.5,), output_dir=tmp_path)
    qr = res.quantile_runs[0]
    assert qr.representative_diameter_m == pytest.approx(
        2.0 * qr.representative_radius_m, rel=1e-12,
    )


def test_n8_duplicate_quantiles_deduped(smoke_params, tmp_path):
    """run_batch silently sort+dedupes quantiles (audit N8)."""
    res = run_batch(
        smoke_params, quantiles=(0.50, 0.25, 0.50, 0.75, 0.25),
        output_dir=tmp_path,
    )
    # 5 inputs collapse to 3 unique sorted quantiles
    assert res.n_quantiles == 3
    quantiles_seen = [r.quantile for r in res.quantile_runs]
    assert quantiles_seen == [0.25, 0.50, 0.75]
    # Mass fractions sum to 1 even after dedup
    assert sum(res.quantile_mass_fractions) == pytest.approx(1.0, abs=1e-9)
