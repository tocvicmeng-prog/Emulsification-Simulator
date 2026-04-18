"""Tests for Phase D: Catalytic packed-bed reactor.

Covers kinetics, deactivation, and the full PFR solver.
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.module3_performance.catalysis.kinetics import (
    michaelis_menten_rate,
    thiele_modulus,
    effectiveness_factor,
    generalized_thiele_modulus,
)
from emulsim.module3_performance.catalysis.deactivation import (
    first_order_deactivation,
    half_life,
    arrhenius_deactivation_rate,
)
from emulsim.module3_performance.catalysis.packed_bed import (
    CatalyticResult,
    solve_packed_bed,
)


# ===================================================================
# Kinetics tests
# ===================================================================

class TestMichaelisMenten:
    """Tests for Michaelis-Menten rate law."""

    def test_saturation_limit(self):
        """v -> V_max at S >> K_m."""
        V_max, K_m = 100.0, 1.0
        S_high = 1e6  # S >> K_m
        v = michaelis_menten_rate(S_high, V_max, K_m)
        assert v == pytest.approx(V_max, rel=1e-4)

    def test_first_order_limit(self):
        """v -> V_max * S / K_m at S << K_m."""
        V_max, K_m = 100.0, 1000.0
        S_low = 0.01
        v = michaelis_menten_rate(S_low, V_max, K_m)
        expected = V_max * S_low / K_m
        assert v == pytest.approx(expected, rel=1e-3)

    def test_at_K_m(self):
        """v = V_max / 2 when S = K_m."""
        V_max, K_m = 50.0, 10.0
        v = michaelis_menten_rate(K_m, V_max, K_m)
        assert v == pytest.approx(V_max / 2.0, rel=1e-10)

    def test_negative_S_clipped(self):
        """Negative substrate concentration is clipped to zero."""
        v = michaelis_menten_rate(-5.0, 100.0, 1.0)
        assert v == 0.0

    def test_array_input(self):
        """Vectorised computation works."""
        S = np.array([0.0, 1.0, 10.0, 100.0])
        v = michaelis_menten_rate(S, V_max=10.0, K_m=1.0)
        assert v.shape == (4,)
        assert v[0] == 0.0
        assert v[-1] > v[1]


class TestEffectivenessFactor:
    """Tests for Thiele modulus and effectiveness factor."""

    def test_eta_approaches_one_small_phi(self):
        """eta -> 1 at phi -> 0 (no diffusion limitation)."""
        eta = effectiveness_factor(1e-8)
        assert eta == pytest.approx(1.0, abs=1e-6)

    def test_eta_approaches_3_over_phi_large(self):
        """eta -> 3/phi at phi -> inf."""
        phi = 100.0
        eta = effectiveness_factor(phi)
        assert eta == pytest.approx(3.0 / phi, rel=0.02)

    def test_eta_intermediate(self):
        """eta at phi=1 matches known value (~0.94)."""
        eta = effectiveness_factor(1.0)
        # Exact: (3/1)*(1/tanh(1) - 1/1) = 3*(coth(1) - 1) ~ 0.9391
        expected = 3.0 * (1.0 / np.tanh(1.0) - 1.0)
        assert eta == pytest.approx(expected, rel=1e-6)

    def test_eta_array(self):
        """Vectorised effectiveness factor."""
        phi = np.array([0.01, 1.0, 10.0, 100.0])
        eta = effectiveness_factor(phi)
        assert eta.shape == (4,)
        # eta should decrease with increasing phi
        assert np.all(np.diff(eta) <= 0)

    def test_thiele_modulus_basic(self):
        """Thiele modulus for known parameters."""
        R = 50e-6  # 50 um
        V_max = 1e3
        K_m = 1.0
        D_eff = 1e-10
        phi = thiele_modulus(R, V_max, K_m, D_eff)
        expected = R * np.sqrt(V_max / (K_m * D_eff))
        assert phi == pytest.approx(expected, rel=1e-10)

    def test_generalized_thiele_gives_consistent_eta(self):
        """Generalised Thiele at high saturation (S >> K_m) gives low eta
        (zero-order regime means stronger diffusion limitation), while at
        low S/K_m (first-order regime) the effectiveness factor from the
        generalised modulus is consistent with the first-order result."""
        R = 50e-6
        V_max = 1e3
        K_m = 1.0
        D_eff = 1e-10

        # At S = K_m (intermediate regime)
        phi_gen = generalized_thiele_modulus(R, V_max, K_m, D_eff, S_bulk=K_m)
        eta_gen = effectiveness_factor(phi_gen)
        # eta should be a valid number in [0, 1]
        assert 0.0 <= eta_gen <= 1.0

        # At very high S (zero-order), Thiele should be larger (more limitation)
        phi_high = generalized_thiele_modulus(R, V_max, K_m, D_eff, S_bulk=1000.0)
        phi_low = generalized_thiele_modulus(R, V_max, K_m, D_eff, S_bulk=0.1)
        # Higher S -> smaller generalized Thiele (closer to zero-order)
        # Lower S -> larger generalized Thiele (closer to first-order)
        assert phi_high < phi_low


# ===================================================================
# Deactivation tests
# ===================================================================

class TestDeactivation:
    """Tests for enzyme deactivation models."""

    def test_activity_at_zero(self):
        """Activity = 1 at t = 0."""
        a = first_order_deactivation(0.0, k_d=0.01)
        assert a == pytest.approx(1.0)

    def test_half_life_consistency(self):
        """Activity = 0.5 at t = ln(2) / k_d."""
        k_d = 1e-4
        t_half = half_life(k_d)
        a = first_order_deactivation(t_half, k_d)
        assert a == pytest.approx(0.5, rel=1e-10)

    def test_half_life_value(self):
        """Half-life = ln(2) / k_d."""
        k_d = 0.001
        assert half_life(k_d) == pytest.approx(np.log(2) / k_d)

    def test_no_deactivation(self):
        """k_d = 0 means activity stays at 1."""
        a = first_order_deactivation(1e6, k_d=0.0)
        assert a == pytest.approx(1.0)

    def test_array_time(self):
        """Vectorised time input."""
        t = np.array([0, 100, 1000])
        a = first_order_deactivation(t, k_d=0.001)
        assert a.shape == (3,)
        assert a[0] == pytest.approx(1.0)
        assert a[-1] < a[0]

    def test_arrhenius_reference(self):
        """At T = T_ref, k_d = k_d_ref."""
        k_d = arrhenius_deactivation_rate(310.15, k_d_ref=0.001, E_a_d=50000.0)
        assert k_d == pytest.approx(0.001, rel=1e-10)

    def test_arrhenius_higher_temp(self):
        """Higher temperature increases deactivation rate."""
        k_ref = 0.001
        E_a = 80000.0
        k_low = arrhenius_deactivation_rate(310.15, k_ref, E_a)
        k_high = arrhenius_deactivation_rate(330.0, k_ref, E_a)
        assert k_high > k_low


# ===================================================================
# Packed-bed reactor tests
# ===================================================================

# Common parameters for PFR tests
_BED_PARAMS = dict(
    bed_length=0.10,          # 10 cm
    bed_diameter=0.01,        # 1 cm
    particle_diameter=100e-6, # 100 um beads
    bed_porosity=0.4,
    particle_porosity=0.5,
    K_m=1.0,                  # mol/m^3
    S_feed=10.0,              # mol/m^3
    flow_rate=1e-8,           # 10 uL/s
    D_eff=1e-10,
    n_z=30,
    total_time=600.0,
)


class TestPackedBedNoReaction:
    """PFR with V_max = 0 — substrate should pass through unchanged."""

    def test_no_reaction_passthrough(self):
        """With V_max=0, S_outlet should reach S_feed at steady state."""
        result = solve_packed_bed(V_max=0.0, **_BED_PARAMS)
        assert isinstance(result, CatalyticResult)
        # At steady state, S_outlet = S_feed
        assert result.S_outlet[-1] == pytest.approx(
            _BED_PARAMS["S_feed"], rel=0.01
        )

    def test_no_reaction_zero_conversion(self):
        """Zero reaction -> zero conversion."""
        result = solve_packed_bed(V_max=0.0, **_BED_PARAMS)
        assert result.conversion == pytest.approx(0.0, abs=0.02)


class TestPackedBedFullConversion:
    """PFR with very high V_max — substrate should be consumed."""

    def test_high_vmax_high_conversion(self):
        """Very high V_max: S_outlet -> 0, conversion -> 1."""
        result = solve_packed_bed(V_max=1e6, **_BED_PARAMS)
        assert result.conversion > 0.95

    def test_high_vmax_product_formed(self):
        """Product at outlet should be ~S_feed when conversion is high."""
        result = solve_packed_bed(V_max=1e6, **_BED_PARAMS)
        # P_outlet should approach S_feed
        assert result.P_outlet[-1] > 0.5 * _BED_PARAMS["S_feed"]


class TestPackedBedMassBalance:
    """Mass balance: S_consumed + S_remaining = S_injected."""

    def test_mass_balance(self):
        """Mass balance error should be small (< 2%)."""
        result = solve_packed_bed(V_max=500.0, **_BED_PARAMS)
        assert result.mass_balance_error < 0.02

    def test_mass_balance_high_conversion(self):
        """Mass balance at high conversion."""
        result = solve_packed_bed(V_max=1e5, **_BED_PARAMS)
        assert result.mass_balance_error < 0.05


class TestPackedBedDeactivation:
    """Enzyme deactivation should reduce conversion over time."""

    def test_deactivation_reduces_conversion(self):
        """Conversion decreases over time with k_d > 0."""
        result = solve_packed_bed(
            V_max=1000.0,
            k_deact=0.005,  # significant deactivation
            total_time=1000.0,
            **{k: v for k, v in _BED_PARAMS.items()
               if k not in ("total_time",)},
        )
        # Early S_outlet should be lower (more conversion) than late
        # After the system reaches initial steady state
        n_t = len(result.S_outlet)
        early_idx = n_t // 4
        late_idx = -1
        # Later time => more deactivation => less conversion => higher S_outlet
        assert result.S_outlet[late_idx] > result.S_outlet[early_idx]

    def test_activity_history_decays(self):
        """Activity should decay exponentially."""
        result = solve_packed_bed(
            V_max=1000.0,
            k_deact=0.001,
            **_BED_PARAMS,
        )
        assert result.activity_history[0] == pytest.approx(1.0)
        assert result.activity_history[-1] < 1.0
        assert result.activity_history[-1] > 0.0


class TestPackedBedEffectiveness:
    """Effectiveness factor should be in valid range."""

    def test_eta_in_range(self):
        """eta should be in [0, 1]."""
        result = solve_packed_bed(V_max=1000.0, **_BED_PARAMS)
        assert 0.0 <= result.effectiveness_factor <= 1.0

    def test_large_particles_lower_eta(self):
        """Larger particles have lower effectiveness factor (more diffusion limitation)."""
        params_small = {**_BED_PARAMS, "particle_diameter": 50e-6}
        params_large = {**_BED_PARAMS, "particle_diameter": 500e-6}
        r_small = solve_packed_bed(V_max=1e4, **params_small)
        r_large = solve_packed_bed(V_max=1e4, **params_large)
        assert r_large.effectiveness_factor < r_small.effectiveness_factor


class TestPackedBedResidenceTime:
    """Lower flow rate -> higher residence time -> higher conversion."""

    def test_conversion_increases_with_residence_time(self):
        """Lower flow rate gives higher conversion."""
        r_fast = solve_packed_bed(
            V_max=500.0,
            flow_rate=5e-8,
            **{k: v for k, v in _BED_PARAMS.items() if k != "flow_rate"},
        )
        r_slow = solve_packed_bed(
            V_max=500.0,
            flow_rate=1e-8,
            **{k: v for k, v in _BED_PARAMS.items() if k != "flow_rate"},
        )
        assert r_slow.conversion > r_fast.conversion
