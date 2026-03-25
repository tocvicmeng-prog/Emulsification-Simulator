"""Tests for the optimisation engine."""

import numpy as np
import pytest

from emulsim.datatypes import SimulationParameters
from emulsim.optimization.objectives import (
    PARAM_BOUNDS,
    PARAM_NAMES,
    LOG_SCALE_INDICES,
    compute_objectives,
)


class TestObjectives:
    def test_param_names_count(self):
        assert len(PARAM_NAMES) == 7

    def test_bounds_shape(self):
        assert PARAM_BOUNDS.shape == (7, 2)

    def test_bounds_ordered(self):
        for i in range(7):
            assert PARAM_BOUNDS[i, 0] < PARAM_BOUNDS[i, 1]

    def test_log_scale_indices_valid(self):
        for idx in LOG_SCALE_INDICES:
            assert 0 <= idx < 7


class TestSearchSpace:
    def test_roundtrip(self):
        """Transform to search space and back should be identity."""
        from emulsim.optimization.engine import _to_search_space, _from_search_space
        x = np.array([10000, 20, 0.7, 363, 0.167, 2.0, 86400.0])
        x_ss = _to_search_space(x)
        x_back = _from_search_space(x_ss)
        np.testing.assert_allclose(x_back, x, rtol=1e-10)

    def test_log_transforms(self):
        from emulsim.optimization.engine import _to_search_space
        x = np.array([10000, 20, 0.7, 363, 0.167, 2.0, 86400.0])
        x_ss = _to_search_space(x)
        # RPM (index 0) should be log10(10000) = 4
        assert x_ss[0] == pytest.approx(4.0)
        # agarose_frac (index 2) should NOT be log-transformed
        assert x_ss[2] == pytest.approx(0.7)


class TestOptimizationEngine:
    def test_engine_initializes(self):
        from emulsim.optimization.engine import OptimizationEngine
        engine = OptimizationEngine(n_initial=3, max_iterations=2)
        assert engine.n_initial == 3
        assert engine.max_iterations == 2

    @pytest.mark.slow
    def test_short_campaign(self):
        """Run a minimal 3+2 campaign to verify the engine works end-to-end."""
        from emulsim.optimization.engine import OptimizationEngine
        engine = OptimizationEngine(n_initial=3, max_iterations=2)
        state = engine.run()
        assert len(state.X_observed) == 5  # 3 init + 2 BO
        assert len(state.Y_observed) == 5
        assert state.pareto_X.shape[0] >= 1
        assert state.iteration >= 1
