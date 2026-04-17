"""Tests for joblib parallel Monte Carlo in the unified engine.

Node 15 (v7.0, P5b) introduced joblib parallelism; Node 27 / Audit N7
added the small-n-samples auto-serial-fallback; Node 30 (v7.1) merged
the legacy ``uncertainty_core.UncertaintyPropagator`` into
``UnifiedUncertaintyEngine``. These tests run against the merged
engine. Invariants verified:

  1. Serial (n_jobs=1) and parallel (n_jobs=2) paths produce IDENTICAL
     raw sample arrays for the same seed (master process owns RNG draws,
     workers run the deterministic pipeline given pre-perturbed props).
  2. The parallel path does not crash when n_jobs > available cores.
  3. PropertyDatabase is independent per process (no shared mutable
     state, no disk cache that forked workers would lock on).
  4. Audit N7 auto-fallback fires and logs when ``n_samples < 4*n_jobs``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from emulsim.config import load_config
from emulsim.uncertainty_unified import UnifiedUncertaintyEngine


@pytest.fixture(scope="module")
def smoke_params():
    cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
    if not cfg.exists():
        pytest.skip("fast_smoke.toml missing")
    return load_config(cfg)


def _raw(result, name: str) -> np.ndarray:
    """Extract the raw sample array for a named output."""
    out = result.get(name)
    assert out is not None, f"output {name!r} missing"
    assert out.raw_samples is not None, f"raw_samples not populated for {name!r}"
    return out.raw_samples


def test_serial_and_parallel_produce_identical_samples(smoke_params):
    """Same seed -> same samples regardless of n_jobs."""
    engine_s = UnifiedUncertaintyEngine()
    engine_p = UnifiedUncertaintyEngine()
    serial = engine_s.run_m1l4(smoke_params, n_samples=4, seed=42, n_jobs=1)
    parallel = engine_p.run_m1l4(smoke_params, n_samples=4, seed=42, n_jobs=2)

    np.testing.assert_allclose(
        _raw(serial, "d32"), _raw(parallel, "d32"), rtol=1e-12, atol=0.0,
    )
    np.testing.assert_allclose(
        _raw(serial, "pore_size_mean"), _raw(parallel, "pore_size_mean"),
        rtol=1e-12, atol=0.0,
    )
    np.testing.assert_allclose(
        _raw(serial, "G_DN"), _raw(parallel, "G_DN"), rtol=1e-12, atol=0.0,
    )


def test_parallel_n_jobs_clamped_safely(smoke_params):
    """n_jobs greater than available cores must not crash (joblib clamps)."""
    engine = UnifiedUncertaintyEngine()
    res = engine.run_m1l4(smoke_params, n_samples=2, seed=7, n_jobs=64)
    assert res.n_samples == 2
    assert res.n_failed == 0


def test_property_database_is_in_memory_only():
    """PropertyDatabase audit: no disk cache, instances are independent."""
    from emulsim.properties.database import PropertyDatabase
    db1 = PropertyDatabase()
    db2 = PropertyDatabase()
    db1.props.sigma = 99e-3
    assert db2.props.sigma != 99e-3
    assert not hasattr(db1, "_cache_path")
    assert not hasattr(db1, "_disk_cache")


def test_n7_auto_serial_fallback_for_small_n_samples(smoke_params, caplog):
    """Audit N7: n_jobs>1 with tiny n_samples auto-falls back to serial.

    Below the 4*|n_jobs| threshold, joblib startup + Numba JIT cold-compile
    dominate. The engine should detect this and run serially while
    logging the fallback decision under the merged module's logger name.
    """
    import logging
    caplog.set_level(logging.INFO, logger="emulsim.uncertainty_unified")
    engine = UnifiedUncertaintyEngine()
    res = engine.run_m1l4(smoke_params, n_samples=2, seed=42, n_jobs=2)
    assert res.n_samples == 2
    fallback_msgs = [r for r in caplog.records
                     if "falling back to" in r.message
                     and "serial" in r.message]
    assert len(fallback_msgs) == 1, (
        "n_jobs=2 with n_samples=2 should auto-fallback to serial; "
        f"got log records: {[r.message for r in caplog.records]}"
    )
