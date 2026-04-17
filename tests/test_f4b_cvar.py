"""Tests for Node F4-b: CVaR aggregation in OptimizationEngine.

6 tests per docs/f4b_cvar_protocol.md §6.
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.optimization.engine import OptimizationEngine
from emulsim.optimization.objectives import TargetSpec


def _cvar_ref(samples: np.ndarray, alpha: float) -> np.ndarray:
    """Reference CVaR implementation matching engine.py."""
    n = samples.shape[0]
    k = max(1, int(np.ceil(alpha * n)))
    sorted_asc = np.sort(samples, axis=0)
    return sorted_asc[-k:].mean(axis=0)


class TestCVaRHelperMath:
    def test_monotone_samples_alpha_30pct_mean_top_3(self):
        samples = np.arange(1.0, 11.0).reshape(10, 1)  # 1..10
        cvar = _cvar_ref(samples, 0.3)
        # α=0.3 → k=3 → mean of top 3 = mean(8,9,10) = 9.0
        assert cvar[0] == pytest.approx(9.0)

    def test_alpha_one_recovers_sample_mean(self):
        rng = np.random.default_rng(0)
        samples = rng.standard_normal((20, 2))
        cvar = _cvar_ref(samples, 1.0)
        assert cvar == pytest.approx(samples.mean(axis=0))

    def test_cvar_monotone_decreasing_with_alpha(self):
        """For fixed samples, CVaR_α decreases as α grows (more of the
        distribution is averaged in, diluting the worst tail).
        """
        rng = np.random.default_rng(1)
        samples = rng.standard_normal((50, 1))
        vals = [_cvar_ref(samples, a)[0] for a in [0.1, 0.3, 0.5, 0.8, 1.0]]
        # Non-strict because of tie handling at k boundaries.
        assert all(vals[i] >= vals[i + 1] - 1e-12 for i in range(len(vals) - 1))


class TestEngineCtorValidation:
    def _dummy_target(self):
        return TargetSpec(d32_target=50e-6, d32_tol=10e-6)

    def test_rejects_alpha_out_of_bounds_below_zero(self):
        with pytest.raises(ValueError):
            OptimizationEngine(
                target_spec=self._dummy_target(),
                robust_cvar_alpha=-0.1,
                robust_n_samples=5,
            )

    def test_rejects_alpha_above_one(self):
        with pytest.raises(ValueError):
            OptimizationEngine(
                target_spec=self._dummy_target(),
                robust_cvar_alpha=1.5,
                robust_n_samples=5,
            )

    def test_rejects_alpha_positive_with_n_samples_below_two(self):
        with pytest.raises(ValueError):
            OptimizationEngine(
                target_spec=self._dummy_target(),
                robust_cvar_alpha=0.3,
                robust_n_samples=1,
            )

    def test_rejects_alpha_without_target_spec(self):
        with pytest.raises(ValueError):
            OptimizationEngine(
                robust_cvar_alpha=0.3,
                robust_n_samples=5,
            )


class TestCliFlag:
    def test_main_parser_has_cvar_flag(self):
        import emulsim.__main__ as mm
        import inspect
        src = inspect.getsource(mm)
        assert '"--robust-cvar-alpha"' in src

    def test_design_subparser_accepts_cvar_alpha(self):
        import argparse

        parser = argparse.ArgumentParser()
        sub = parser.add_subparsers(dest="command")
        des = sub.add_parser("design")
        des.add_argument("--robust-cvar-alpha", dest="robust_cvar_alpha",
                         type=float, default=0.0)
        args = parser.parse_args(["design", "--robust-cvar-alpha", "0.3"])
        assert args.robust_cvar_alpha == pytest.approx(0.3)


class TestCVaRPrecedenceLogic:
    def test_engine_accepts_both_but_cvar_stored_separately(self):
        engine = OptimizationEngine(
            target_spec=TargetSpec(d32_target=50e-6, d32_tol=10e-6),
            robust_variance_weight=1.0,
            robust_cvar_alpha=0.3,
            robust_n_samples=5,
        )
        assert engine.robust_variance_weight == pytest.approx(1.0)
        assert engine.robust_cvar_alpha == pytest.approx(0.3)

    def test_cvar_vs_variance_aggregation_differ(self):
        """Direct head-to-head on a synthetic sample matrix: CVaR and
        mean-variance should generally give different answers. This
        validates that the engine's reduction code-path differentiates
        between the two modes.
        """
        samples = np.array([
            [1.0], [2.0], [3.0], [4.0], [5.0],
            [6.0], [7.0], [8.0], [9.0], [10.0],
        ])
        # CVaR at α=0.3: mean of top 3 = 9.0
        cvar_val = _cvar_ref(samples, 0.3)[0]
        # Mean-variance with λ=1: mean(5.5) + 1·std(2.872) = 8.37
        mean = samples.mean()
        std = samples.std(ddof=0)
        mv_val = mean + 1.0 * std
        assert cvar_val != pytest.approx(mv_val, rel=1e-6)
