"""Tests for Node 15 (v7.0, P5b): joblib parallel Monte Carlo.

Acceptance for Node 15:
  1. Serial (n_jobs=1) and parallel (n_jobs=2) paths produce IDENTICAL
     sample arrays for the same seed (the master process owns RNG draws,
     workers do deterministic pipeline runs).
  2. The parallel path doesn't crash when n_jobs > available cores.
  3. PropertyDatabase audit: instances are independent across processes
     (no shared mutable state, no disk cache that needs locking).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from emulsim.config import load_config
from emulsim.uncertainty_core import UncertaintyPropagator


@pytest.fixture(scope="module")
def smoke_params():
    cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
    if not cfg.exists():
        pytest.skip("fast_smoke.toml missing")
    return load_config(cfg)


def test_serial_and_parallel_produce_identical_samples(smoke_params):
    """Same seed -> same samples regardless of n_jobs."""
    serial = UncertaintyPropagator(n_samples=4, seed=42, n_jobs=1).run(smoke_params)
    parallel = UncertaintyPropagator(n_samples=4, seed=42, n_jobs=2).run(smoke_params)

    np.testing.assert_allclose(serial.d32_samples, parallel.d32_samples,
                                rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(serial.pore_samples, parallel.pore_samples,
                                rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(serial.G_DN_samples, parallel.G_DN_samples,
                                rtol=1e-12, atol=0.0)


def test_parallel_n_jobs_clamped_safely(smoke_params):
    """n_jobs greater than available cores must not crash (joblib clamps)."""
    res = UncertaintyPropagator(n_samples=2, seed=7, n_jobs=64).run(smoke_params)
    assert res.n_samples == 2
    assert res.n_failed == 0


def test_property_database_is_in_memory_only():
    """PropertyDatabase audit: no disk cache, instances are independent."""
    from emulsim.properties.database import PropertyDatabase
    db1 = PropertyDatabase()
    db2 = PropertyDatabase()
    # Mutating one instance does not affect the other (no shared state).
    db1.props.sigma = 99e-3
    assert db2.props.sigma != 99e-3
    # No file cache attribute exists -> nothing to lock for forked workers.
    assert not hasattr(db1, "_cache_path")
    assert not hasattr(db1, "_disk_cache")
