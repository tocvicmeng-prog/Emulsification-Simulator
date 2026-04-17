"""Tests for Node F2 Phase 1: digital-twin EnKF replay harness.

10 tests per docs/f2_digital_twin_protocol.md §7.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from emulsim.datatypes import ModelEvidenceTier
from emulsim.digital_twin import (
    DigitalTwinTrace,
    Observation,
    enkf_update,
    load_trace,
    run_replay,
    save_trace,
)


# ─── EnKF linear-Gaussian convergence ───────────────────────────────


class TestEnKFLinearGaussianConvergence:
    def test_persistence_model_mean_converges_to_obs(self):
        """With f(x) = x and h(x) = x, applying enkf_update N times
        with the same observation drives the posterior mean toward
        the observation. Final mean should be within a few % of y_obs.
        """
        rng = np.random.default_rng(0)
        N = 200
        y_obs = 5.0
        # Prior: N(0, 1)
        x = rng.normal(loc=0.0, scale=1.0, size=(N, 1))
        for _ in range(20):
            y_fc = x[:, 0].copy()
            x = enkf_update(x, y_fc, y_obs, R=0.01, rng=rng)
        mean_final = float(x.mean())
        assert abs(mean_final - y_obs) < 0.3

    def test_zero_obs_noise_collapses_mean_to_obs(self):
        """R = 0 + no perturbation: ensemble mean exactly matches
        the observation after one update (for h(x) = x).
        """
        rng = np.random.default_rng(1)
        N = 50
        y_obs = 3.0
        x = rng.normal(loc=0.0, scale=1.0, size=(N, 1))
        y_fc = x[:, 0].copy()
        x_post = enkf_update(x, y_fc, y_obs, R=0.0, rng=rng)
        assert float(x_post.mean()) == pytest.approx(y_obs, abs=1e-6)


class TestEnKFInflation:
    def test_inflation_grows_prior_spread(self):
        """Applying inflation > 1 before a low-information update
        should leave the posterior with strictly more spread than the
        same update with inflation = 1 (for R large relative to P_yy).
        """
        rng_a = np.random.default_rng(42)
        rng_b = np.random.default_rng(42)
        N = 80
        x = np.random.default_rng(0).normal(size=(N, 1))
        y_fc = x[:, 0].copy()
        # Huge R — observation provides almost no information; inflation
        # dominates.
        post_noinf = enkf_update(x, y_fc, 0.0, R=1e6, rng=rng_a, inflation=1.0)
        post_inf = enkf_update(x, y_fc, 0.0, R=1e6, rng=rng_b, inflation=1.5)
        std_noinf = float(post_noinf.std())
        std_inf = float(post_inf.std())
        assert std_inf > std_noinf


class TestEnKFValidation:
    def test_N_below_2_raises(self):
        with pytest.raises(ValueError):
            enkf_update(
                np.array([[1.0]]), np.array([1.0]),
                1.0, R=0.1,
            )

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError):
            enkf_update(
                np.zeros((10, 2)),
                np.zeros(5),      # wrong size
                0.0, R=0.1,
            )

    def test_negative_R_raises(self):
        with pytest.raises(ValueError):
            enkf_update(
                np.zeros((5, 2)), np.zeros(5),
                0.0, R=-1.0,
            )


# ─── DigitalTwinTrace schema ────────────────────────────────────────


class TestTraceSchema:
    def test_roundtrip_preserves_fields(self, tmp_path):
        obs = [
            Observation(t=0.0, name="d32", value=1.5e-5, noise_std=1e-7),
            Observation(t=60.0, name="d32", value=1.3e-5, noise_std=1e-7),
        ]
        trace = DigitalTwinTrace(
            trace_id="run_001",
            process_description="test replay",
            observations=obs,
            metadata={"rpm": 8000.0, "vessel": "glass_beaker"},
        )
        p = tmp_path / "trace.json"
        save_trace(trace, p)
        loaded = load_trace(p)
        assert loaded.trace_id == "run_001"
        assert loaded.process_description == "test replay"
        assert len(loaded.observations) == 2
        assert loaded.observations[0].t == pytest.approx(0.0)
        assert loaded.observations[1].value == pytest.approx(1.3e-5)
        assert loaded.metadata["rpm"] == pytest.approx(8000.0)

    def test_out_of_order_observations_sorted_on_load(self, tmp_path):
        raw = {
            "trace_id": "unsorted",
            "process_description": "test",
            "observations": [
                {"t": 30.0, "name": "d32", "value": 1.0, "noise_std": 0.01},
                {"t": 10.0, "name": "d32", "value": 2.0, "noise_std": 0.01},
                {"t": 20.0, "name": "d32", "value": 1.5, "noise_std": 0.01},
            ],
            "metadata": {},
        }
        p = tmp_path / "unsorted.json"
        p.write_text(json.dumps(raw))
        loaded = load_trace(p)
        times = [o.t for o in loaded.observations]
        assert times == [10.0, 20.0, 30.0]


# ─── Replay harness ──────────────────────────────────────────────────


class TestReplayHarness:
    def _persistence_model(self, x, dt):
        return x  # f(x) = x

    def _identity_observation(self, x, name):
        return x[:, 0]

    def test_returns_trajectory_of_correct_shape(self):
        trace = DigitalTwinTrace(
            trace_id="t", process_description="p",
            observations=[
                Observation(t=10.0, name="d32", value=1.0, noise_std=0.05),
                Observation(t=20.0, name="d32", value=0.8, noise_std=0.05),
                Observation(t=30.0, name="d32", value=0.6, noise_std=0.05),
            ],
        )
        rng = np.random.default_rng(0)
        x0 = rng.normal(size=(30, 1))
        result = run_replay(
            trace, x0,
            state_transition=self._persistence_model,
            observation_operator=self._identity_observation,
        )
        assert len(result.times) == 3
        assert len(result.means) == 3
        assert len(result.stds) == 3
        assert result.final_ensemble.shape == (30, 1)
        assert (
            result.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )

    def test_no_observations_returns_empty_trajectory(self):
        trace = DigitalTwinTrace(
            trace_id="empty", process_description="p",
            observations=[],
        )
        x0 = np.zeros((5, 2))
        result = run_replay(
            trace, x0,
            state_transition=self._persistence_model,
            observation_operator=self._identity_observation,
        )
        assert result.times == []
        assert result.means == []
        assert np.array_equal(result.final_ensemble, x0)

    def test_multi_step_replay_shrinks_spread(self):
        """Repeated informative observations → posterior spread
        shrinks across the trajectory."""
        trace = DigitalTwinTrace(
            trace_id="shrink", process_description="p",
            observations=[
                Observation(
                    t=float(10 * k), name="d32", value=2.0, noise_std=0.05,
                )
                for k in range(1, 11)   # 10 observations at t=10,20,...,100
            ],
        )
        rng = np.random.default_rng(7)
        x0 = rng.normal(loc=0.0, scale=2.0, size=(100, 1))
        result = run_replay(
            trace, x0,
            state_transition=self._persistence_model,
            observation_operator=self._identity_observation,
            inflation=1.0,
            store_ensembles=True,
        )
        std_first = float(result.stds[0][0])
        std_last = float(result.stds[-1][0])
        assert std_last < std_first
        # Final mean should be close to the observed 2.0.
        assert float(result.means[-1][0]) == pytest.approx(2.0, abs=0.2)
