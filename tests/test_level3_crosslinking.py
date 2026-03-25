"""Tests for Level 3 crosslinking kinetics."""

import numpy as np
import pytest

from emulsim.datatypes import MaterialProperties, SimulationParameters
from emulsim.level3_crosslinking.solver import (
    available_amine_concentration,
    genipin_rate_constant,
    solve_crosslinking,
)


class TestCrosslinkingChemistry:
    def test_rate_constant_positive(self):
        k = genipin_rate_constant(310.15, 1.33e4, 52000.0)
        assert k > 0

    def test_rate_constant_increases_with_T(self):
        k1 = genipin_rate_constant(300.0, 1.33e4, 52000.0)
        k2 = genipin_rate_constant(340.0, 1.33e4, 52000.0)
        assert k2 > k1

    def test_amine_concentration(self):
        """18 kg/m³ chitosan with 90% DDA, M_GlcN=161.16 g/mol."""
        nh2 = available_amine_concentration(18.0, 0.9, 161.16)
        assert nh2 > 0
        assert 50 < nh2 < 200  # reasonable mol/m³ range


class TestCrosslinkingSolver:
    def test_solver_runs(self):
        params = SimulationParameters()
        props = MaterialProperties()
        result = solve_crosslinking(params, props)
        assert len(result.t_array) > 0
        assert result.p_final >= 0

    def test_crosslink_density_increases(self):
        params = SimulationParameters()
        props = MaterialProperties()
        result = solve_crosslinking(params, props)
        # X should increase monotonically
        assert np.all(np.diff(result.X_array) >= -1e-15)

    def test_modulus_increases(self):
        params = SimulationParameters()
        props = MaterialProperties()
        result = solve_crosslinking(params, props)
        # G_chitosan should increase with crosslinking
        assert result.G_chitosan_array[-1] >= result.G_chitosan_array[0]

    def test_conversion_bounded(self):
        params = SimulationParameters()
        props = MaterialProperties()
        result = solve_crosslinking(params, props)
        assert 0.0 <= result.p_final <= 1.0

    def test_mesh_size_decreases(self):
        """More crosslinking → smaller mesh."""
        params = SimulationParameters()
        props = MaterialProperties()
        result = solve_crosslinking(params, props)
        assert result.xi_array[-1] <= result.xi_array[0]

    def test_higher_genipin_faster(self):
        """More genipin should give higher conversion in same time."""
        params_low = SimulationParameters()
        params_low.formulation.c_genipin = 1.0
        params_low.formulation.t_crosslink = 3600.0

        params_high = SimulationParameters()
        params_high.formulation.c_genipin = 10.0
        params_high.formulation.t_crosslink = 3600.0

        props = MaterialProperties()
        r_low = solve_crosslinking(params_low, props)
        r_high = solve_crosslinking(params_high, props)
        assert r_high.p_final >= r_low.p_final
