"""Physics tests for emulsification_derivation + rpm suggestion module."""

from __future__ import annotations

import math


from emulsim.properties.emulsification_derivation import (
    d32_from_weber,
    kolmogorov_scale,
    rpm_for_target_d32,
    weber_for_target_d32,
    weber_number,
)
from emulsim.suggestions import rpm as rpm_module
from emulsim.suggestions.types import SuggestionContext


class TestWeberScaling:
    def test_weber_scales_with_N_squared(self):
        We1 = weber_number(N_rps=10, D=0.05, rho_c=850, sigma=0.03)
        We2 = weber_number(N_rps=20, D=0.05, rho_c=850, sigma=0.03)
        assert math.isclose(We2 / We1, 4.0)

    def test_d32_decreases_with_higher_weber(self):
        We_low = 100
        We_high = 1000
        d_low = d32_from_weber(We_low, D=0.05, phi_d=0.05)
        d_high = d32_from_weber(We_high, D=0.05, phi_d=0.05)
        assert d_high < d_low


class TestRPMRoundTrip:
    def test_weber_forward_inverse_round_trip(self):
        d_in = 3e-6
        We = weber_for_target_d32(d_in, D=0.05, phi_d=0.05)
        d_out = d32_from_weber(We, D=0.05, phi_d=0.05)
        assert math.isclose(d_in, d_out, rel_tol=1e-9)

    def test_target_rpm_reproduces_target_d32(self):
        """Invert to RPM, then run forward -> should recover target d32 within 1 %."""
        target_d = 3e-6
        result = rpm_for_target_d32(
            target_d32=target_d, D=0.05, rho_c=850.0, mu_c=0.005,
            sigma=0.03, phi_d=0.05,
        )
        N_rps = result.nominal / 60.0
        We = weber_number(N_rps, D=0.05, rho_c=850.0, sigma=0.03)
        d_check = d32_from_weber(We, D=0.05, phi_d=0.05)
        assert math.isclose(target_d, d_check, rel_tol=0.01)


class TestRPMTargetSemantics:
    def _target(self, **overrides):
        kw = dict(
            target_d32=3e-6, D=0.05, rho_c=850.0, mu_c=0.005,
            sigma=0.03, phi_d=0.05,
        )
        kw.update(overrides)
        return rpm_for_target_d32(**kw)

    def test_nominal_between_bounds(self):
        t = self._target()
        assert t.min <= t.nominal <= t.max

    def test_smaller_d32_needs_higher_rpm(self):
        t_small = self._target(target_d32=1e-6)
        t_large = self._target(target_d32=10e-6)
        assert t_small.nominal > t_large.nominal

    def test_kolmogorov_scales_inversely_with_rpm(self):
        k_low = kolmogorov_scale(N_rps=10, D=0.05, rho_c=850, mu_c=0.005)
        k_high = kolmogorov_scale(N_rps=100, D=0.05, rho_c=850, mu_c=0.005)
        assert k_high < k_low


class TestRPMModule:
    def _ctx(self, **overrides):
        defaults = dict(
            family="agarose_chitosan",
            d32_actual=6e-6,  # 2x target -> should trigger increase_rpm
            d50_actual=8e-6,
            pore_actual=100e-9, l2_mode="ch_2d",
            cooling_rate_effective=0.17, p_final=0.25,
            G_DN_actual=10_000.0, target_d32=3e-6, target_pore=100e-9, target_G=10_000.0,
            rpm=4000.0, T_oil=363.15, cooling_rate_input=0.17,
            c_agarose=42.0, c_chitosan=18.0,
            c_crosslinker_mM=2.0, crosslinker_key="genipin",
            rho_oil=850.0, mu_oil=0.005, rho_d=1000.0, cp_d=4180.0,
            k_oil=0.15, h_coeff=500.0, T_bath=293.15, T_gel=311.15,
            DDA=0.85, M_GlcN=161.16, f_bridge=0.4,
            impeller_D=0.05, phi_d=0.05, run_id="",
        )
        defaults.update(overrides)
        return SuggestionContext(**defaults)  # type: ignore[arg-type]

    def test_generate_produces_increase_when_too_large(self):
        s = rpm_module.generate(self._ctx())
        assert s is not None
        assert s.key == "increase_rpm"

    def test_generate_produces_decrease_when_too_small(self):
        ctx = self._ctx(d32_actual=1e-6)  # 1/3 of target
        s = rpm_module.generate(ctx)
        assert s is not None
        assert s.key == "decrease_rpm"

    def test_target_unit_is_rpm(self):
        target = rpm_module.derive_target(self._ctx())
        assert target.unit == "rpm"
        assert target.nominal > 0
