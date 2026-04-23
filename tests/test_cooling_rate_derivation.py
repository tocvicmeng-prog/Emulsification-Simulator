"""Physics tests for thermal_derivation + cooling_rate suggestion module."""

from __future__ import annotations

import math

import pytest

from emulsim.properties.thermal_derivation import (
    bead_biot_number,
    cooling_rate_for_target_pore,
    dwell_for_target_pore,
    effective_cooling_rate,
    lumped_cooling_time,
    pore_from_dwell,
    spinodal_dwell_time,
)
from emulsim.suggestions import cooling_rate as cr_module
from emulsim.suggestions.types import SuggestionContext


class TestCoolingRateKinematics:
    def test_biot_scales_linearly_with_h(self):
        bi1 = bead_biot_number(d=200e-6, h_oil=500, k_oil=0.15)
        bi2 = bead_biot_number(d=200e-6, h_oil=1000, k_oil=0.15)
        assert math.isclose(bi2 / bi1, 2.0)

    def test_biot_scales_linearly_with_d(self):
        bi1 = bead_biot_number(d=100e-6, h_oil=500, k_oil=0.15)
        bi2 = bead_biot_number(d=200e-6, h_oil=500, k_oil=0.15)
        assert math.isclose(bi2 / bi1, 2.0)

    def test_tau_th_scales_linearly_with_d(self):
        t1 = lumped_cooling_time(d=100e-6, rho_d=1000, cp_d=4180, h_oil=500)
        t2 = lumped_cooling_time(d=200e-6, rho_d=1000, cp_d=4180, h_oil=500)
        assert math.isclose(t2 / t1, 2.0)

    def test_effective_rate_positive(self):
        rate = effective_cooling_rate(T_oil=363.15, T_bath=293.15, tau_th=10.0)
        assert rate > 0

    def test_effective_rate_raises_on_zero_tau(self):
        with pytest.raises(ValueError):
            effective_cooling_rate(T_oil=363.15, T_bath=293.15, tau_th=0.0)


class TestCoolingRateRoundTrip:
    """Forward then inverse should recover the input (within prefactor precision)."""

    def test_pore_forward_inverse(self):
        dwell_in = 10.0
        pore = pore_from_dwell(dwell_in)
        dwell_out = dwell_for_target_pore(pore)
        assert math.isclose(dwell_in, dwell_out, rel_tol=1e-9)

    def test_dwell_forward_inverse(self):
        rate_in = 0.2
        dwell = spinodal_dwell_time(rate_in)
        rate_out = 2.0 * 3.0 / dwell  # dT_band default = 3 K
        assert math.isclose(rate_in, rate_out, rel_tol=1e-9)


class TestCoolingRateTargetSemantics:
    def _target(self, **overrides):
        kw = dict(
            target_pore=100e-9,
            d_bead=200e-6,
            T_oil=363.15,
            T_bath=293.15,
            rho_d=1000.0,
            cp_d=4180.0,
            h_oil=500.0,
            k_oil=0.15,
            l2_mode="ch_2d",
        )
        kw.update(overrides)
        return cooling_rate_for_target_pore(**kw)

    def test_nominal_between_min_and_max(self):
        t = self._target()
        assert t.min <= t.nominal <= t.max

    def test_smaller_pore_needs_faster_rate(self):
        t_small = self._target(target_pore=50e-9)
        t_large = self._target(target_pore=200e-9)
        assert t_small.nominal > t_large.nominal

    def test_l2_mode_empirical_flags_qualitative_tier(self):
        t = self._target(l2_mode="empirical")
        assert t.confidence_tier == "QUALITATIVE_TREND"

    def test_l2_mode_ch_2d_is_semi_quantitative(self):
        t = self._target(l2_mode="ch_2d")
        assert t.confidence_tier == "SEMI_QUANTITATIVE"

    def test_high_biot_gets_flagged(self):
        t = self._target(h_oil=10_000)  # Bi >> 0.1
        assert t.limited_by == "biot_invalid"


class TestCoolingRateModule:
    def _ctx(self, **overrides):
        defaults = dict(
            family="agarose_chitosan",
            d32_actual=2e-6, d50_actual=3e-6,
            pore_actual=300e-9,  # off-target by 3x
            l2_mode="ch_2d",
            cooling_rate_effective=0.17,
            p_final=0.25, G_DN_actual=10_000.0,
            target_d32=2e-6, target_pore=100e-9, target_G=10_000.0,
            rpm=8000.0, T_oil=363.15, cooling_rate_input=0.17,
            c_agarose=42.0, c_chitosan=18.0,
            c_crosslinker_mM=2.0, crosslinker_key="genipin",
            rho_oil=850.0, mu_oil=0.005, rho_d=1000.0, cp_d=4180.0,
            k_oil=0.15, h_coeff=500.0, T_bath=293.15, T_gel=311.15,
            DDA=0.85, M_GlcN=161.16, f_bridge=0.4,
            impeller_D=0.05, phi_d=0.05, run_id="",
        )
        defaults.update(overrides)
        return SuggestionContext(**defaults)  # type: ignore[arg-type]

    def test_generate_emits_when_pore_off_target(self):
        s = cr_module.generate(self._ctx())
        assert s is not None
        assert s.key == "adjust_cooling_rate"

    def test_generate_none_when_pore_on_target(self):
        ctx = self._ctx(pore_actual=100e-9)  # exactly on target
        assert cr_module.generate(ctx) is None

    def test_empirical_l2_gives_qualitative_target_with_reason(self):
        ctx = self._ctx(l2_mode="empirical")
        target = cr_module.derive_target(ctx)
        assert target.is_qualitative_only is True
        assert target.qualitative_reason != ""
        assert target.nominal == 0.0

    def test_ch_2d_l2_gives_numeric_target(self):
        ctx = self._ctx(l2_mode="ch_2d")
        target = cr_module.derive_target(ctx)
        assert target.is_qualitative_only is False
        assert target.nominal > 0.0
        assert target.min <= target.nominal <= target.max
