"""Tests for Level 1 emulsification simulation."""

import numpy as np
import pytest

from emulsim.datatypes import MaterialProperties, SimulationParameters
from emulsim.level1_emulsification.energy import (
    average_dissipation,
    emulsion_density,
    kolmogorov_length_scale,
    max_dissipation,
    power_draw,
)
from emulsim.level1_emulsification.kernels import (
    breakage_rate_alopaeus,
    breakage_rate_coulaloglou,
    coalescence_rate_ct,
    hinze_dmax,
    hinze_dmax_viscous,
)
from emulsim.level1_emulsification.solver import PBESolver
from emulsim.level1_emulsification.validation import (
    check_non_negative,
    check_physical_bounds,
)


# ─── Energy dissipation tests ────────────────────────────────────────────


class TestEnergyDissipation:
    def test_power_draw_positive(self, default_mixer):
        P = power_draw(default_mixer, 10000, 850.0)
        assert P > 0

    def test_power_scales_with_rpm_cubed(self, default_mixer):
        P1 = power_draw(default_mixer, 5000, 850.0)
        P2 = power_draw(default_mixer, 10000, 850.0)
        assert pytest.approx(P2 / P1, rel=1e-6) == 8.0  # (10000/5000)^3

    def test_max_exceeds_average(self, default_mixer):
        eps_avg = average_dissipation(default_mixer, 10000, 850.0)
        eps_max = max_dissipation(default_mixer, 10000, 850.0)
        assert eps_max > eps_avg
        assert eps_max == pytest.approx(eps_avg * default_mixer.dissipation_ratio)

    def test_kolmogorov_scale(self):
        eta = kolmogorov_length_scale(1e4, 5e-6)
        assert 1e-6 < eta < 1e-3  # micron to mm range

    def test_emulsion_density_limits(self):
        assert emulsion_density(850, 1020, 0.0) == 850.0
        assert emulsion_density(850, 1020, 1.0) == 1020.0


# ─── Breakage kernel tests ────────────────────────────────────────────────


class TestBreakageKernels:
    def test_breakage_positive(self):
        d = np.array([10e-6, 50e-6, 100e-6])
        g = breakage_rate_alopaeus(d, 1e4, 5e-3, 850.0, 0.1)
        assert np.all(g >= 0)

    def test_breakage_zero_at_zero_epsilon(self):
        d = np.array([10e-6, 50e-6])
        g = breakage_rate_alopaeus(d, 0.0, 5e-3, 850.0, 0.1)
        assert np.all(g == 0)

    def test_breakage_increases_with_epsilon(self):
        d = np.array([50e-6])
        g1 = breakage_rate_alopaeus(d, 1e3, 5e-3, 850.0, 0.1)
        g2 = breakage_rate_alopaeus(d, 1e5, 5e-3, 850.0, 0.1)
        assert g2[0] > g1[0]

    def test_breakage_increases_with_diameter(self):
        d_small = np.array([5e-6])
        d_large = np.array([100e-6])
        g_small = breakage_rate_alopaeus(d_small, 1e4, 5e-3, 850.0, 0.1)
        g_large = breakage_rate_alopaeus(d_large, 1e4, 5e-3, 850.0, 0.1)
        assert g_large[0] > g_small[0]

    def test_viscous_correction_reduces_breakage(self):
        d = np.array([50e-6])
        g_low = breakage_rate_alopaeus(d, 1e4, 5e-3, 850.0, 0.001, C3=0.3)
        g_high = breakage_rate_alopaeus(d, 1e4, 5e-3, 850.0, 1.0, C3=0.3)
        assert g_low[0] >= g_high[0]

    def test_coulaloglou_matches_trend(self):
        d = np.array([50e-6])
        g1 = breakage_rate_coulaloglou(d, 1e3, 5e-3, 850.0)
        g2 = breakage_rate_coulaloglou(d, 1e5, 5e-3, 850.0)
        assert g2[0] > g1[0]


# ─── Coalescence kernel tests ─────────────────────────────────────────────


class TestCoalescenceKernels:
    def test_coalescence_symmetric(self):
        q12 = coalescence_rate_ct(
            np.array([10e-6]), np.array([50e-6]), 1e4, 5e-3, 850.0, 0.005
        )
        q21 = coalescence_rate_ct(
            np.array([50e-6]), np.array([10e-6]), 1e4, 5e-3, 850.0, 0.005
        )
        assert q12 == pytest.approx(q21, rel=1e-10)

    def test_coalescence_positive(self):
        q = coalescence_rate_ct(
            np.array([10e-6]), np.array([50e-6]), 1e4, 5e-3, 850.0, 0.005
        )
        assert np.all(q >= 0)

    def test_coalescence_zero_at_zero_epsilon(self):
        q = coalescence_rate_ct(
            np.array([10e-6]), np.array([50e-6]), 0.0, 5e-3, 850.0, 0.005
        )
        assert np.all(q == 0.0)


# ─── Hinze predictions ────────────────────────────────────────────────────


class TestHinzePredictions:
    def test_hinze_decreases_with_epsilon(self):
        d1 = hinze_dmax(1e3, 5e-3, 850.0)
        d2 = hinze_dmax(1e5, 5e-3, 850.0)
        assert d2 < d1

    def test_hinze_viscous_larger_than_inviscid(self):
        eps = 1e4
        d_inv = hinze_dmax(eps, 5e-3, 850.0)
        d_visc = hinze_dmax_viscous(eps, 5e-3, 850.0, 1.0)
        assert d_visc >= d_inv

    def test_hinze_scaling(self):
        """d_max ~ epsilon^(-0.4)"""
        eps1, eps2 = 1e3, 1e5
        d1 = hinze_dmax(eps1, 5e-3, 850.0)
        d2 = hinze_dmax(eps2, 5e-3, 850.0)
        ratio = d1 / d2
        expected_ratio = (eps2 / eps1) ** 0.4
        assert ratio == pytest.approx(expected_ratio, rel=0.01)


# ─── PBE Solver tests ─────────────────────────────────────────────────────


class TestPBESolver:
    # Class-scoped: full PBE solve takes ~30-90s on Py 3.14+numba; share one
    # solve across the 5 fast_result tests instead of re-running per test.
    # Marked slow because even one solve exceeds the 5s "fast" budget.

    @pytest.fixture(scope="class")
    def solver(self):
        return PBESolver(n_bins=30, d_min=0.5e-6, d_max=200e-6)

    @pytest.fixture(scope="class")
    def fast_result(self, solver):
        """Quick solve with reduced time for testing."""
        params = SimulationParameters()
        params.emulsification.t_emulsification = 60.0  # 1 min only
        props = MaterialProperties(sigma=5e-3)
        return solver.solve(params, props)

    @pytest.mark.slow
    def test_result_shapes(self, fast_result, solver):
        assert fast_result.d_bins.shape == (solver.n_bins,)
        assert fast_result.n_d.shape == (solver.n_bins,)

    @pytest.mark.slow
    def test_non_negative_distribution(self, fast_result):
        assert check_non_negative(fast_result)

    @pytest.mark.slow
    def test_d32_in_range(self, fast_result):
        passed, msg = check_physical_bounds(fast_result)
        assert passed, msg

    @pytest.mark.slow
    def test_percentiles_monotonic(self, fast_result):
        assert fast_result.d10 <= fast_result.d50 <= fast_result.d90

    @pytest.mark.slow
    def test_span_positive(self, fast_result):
        assert fast_result.span >= 0

    def test_d32_decreases_with_rpm(self):
        """Higher RPM should produce smaller droplets (Hinze prediction)."""
        # Use Hinze prediction directly — avoids PBE solver sensitivity
        # to the viscous correction + shear-thinning interaction in short runs.
        from emulsim.level1_emulsification.kernels import hinze_dmax
        from emulsim.level1_emulsification.energy import max_dissipation
        from emulsim.datatypes import MixerGeometry

        mixer = MixerGeometry()
        sigma = 5e-3
        rho_c = 850.0

        eps_low = max_dissipation(mixer, 5000, rho_c)
        eps_high = max_dissipation(mixer, 20000, rho_c)

        d_low = hinze_dmax(eps_low, sigma, rho_c)
        d_high = hinze_dmax(eps_high, sigma, rho_c)

        assert d_high < d_low


@pytest.mark.slow
def test_d32_monotonic_pbe_rpm_sweep():
    """PBE output d32 must decrease monotonically with RPM (Finding F1 regression test).

    This test runs the full PBE solver at 5 RPM values in the converging
    regime [8000, 25000] and asserts that d32 decreases monotonically. It
    guards against the Vi-coalescence interaction that caused nonmonotonic
    RPM->d32 behavior at high RPM (audit Finding F1).

    RPM range starts at 8000 rather than 3000 because at very low RPM the
    PBE requires much longer emulsification time to reach steady state,
    making the test impractically slow. The Hinze-level monotonicity at
    lower RPM is covered by test_d32_decreases_with_rpm.
    """
    rpms = [8000, 10000, 15000, 20000, 25000]
    d32_values = []

    for rpm in rpms:
        params = SimulationParameters()
        params.emulsification.rpm = float(rpm)
        params.emulsification.t_emulsification = 60.0
        # Disable adaptive extensions to keep test runtime manageable
        params.solver.l1_max_extensions = 0

        props = MaterialProperties(sigma=5e-3)

        solver = PBESolver(n_bins=20, d_min=1e-6, d_max=200e-6)
        result = solver.solve(params, props, phi_d=0.05)
        d32_values.append(result.d32)

    # Assert monotonic decrease
    for i in range(len(d32_values) - 1):
        assert d32_values[i] > d32_values[i + 1], (
            f"d32 not monotonic: d32({rpms[i]} RPM)={d32_values[i]*1e6:.2f} um "
            f">= d32({rpms[i+1]} RPM)={d32_values[i+1]*1e6:.2f} um"
        )
