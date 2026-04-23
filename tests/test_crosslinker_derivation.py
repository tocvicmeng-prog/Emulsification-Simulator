"""Physics tests for crosslink_derivation + crosslinker + polymer modules."""

from __future__ import annotations

import math

from emulsim.properties.crosslink_derivation import (
    alpha_for_target_G_DN,
    crosslinker_conc_for_target_G,
    p_for_target_G_chit,
)
from emulsim.suggestions import crosslinker as xl_mod
from emulsim.suggestions import polymer as pol_mod
from emulsim.suggestions.types import SuggestionContext


class TestPForTargetG:
    def test_p_scales_with_G(self):
        G_low = 1_000.0
        G_high = 10_000.0
        p_low = p_for_target_G_chit(G_low, c_chitosan=18.0, DDA=0.85,
                                     M_GlcN=161.16, f_bridge=0.4, T=310.15)
        p_high = p_for_target_G_chit(G_high, c_chitosan=18.0, DDA=0.85,
                                      M_GlcN=161.16, f_bridge=0.4, T=310.15)
        assert p_high > p_low
        # Rubber elasticity: G is linear in nu_e which is linear in p ->
        # ratio of p's should equal ratio of G's.
        assert math.isclose(p_high / p_low, G_high / G_low, rel_tol=1e-6)


class TestCrosslinkerTarget:
    def _target(self, **overrides):
        kw = dict(
            target_G_chit=5_000.0, c_chitosan=18.0, DDA=0.85,
            M_GlcN=161.16, f_bridge=0.4, T=310.15,
        )
        kw.update(overrides)
        return crosslinker_conc_for_target_G(**kw)

    def test_nominal_between_bounds(self):
        t = self._target()
        assert t.min <= t.nominal <= t.max

    def test_higher_target_G_needs_more_crosslinker(self):
        low = self._target(target_G_chit=1_000.0)
        high = self._target(target_G_chit=10_000.0)
        assert high.nominal > low.nominal

    def test_stoichiometry_flag_at_extreme(self):
        # G so high that required p > 1 -> stoichiometrically infeasible
        t = self._target(target_G_chit=1e8)
        assert t.limited_by == "stoichiometry"


class TestPolymerScaling:
    def test_alpha_less_than_1_when_G_target_below_current(self):
        r = alpha_for_target_G_DN(
            G_DN_current=20_000.0, target_G_DN=10_000.0,
            c_agarose_current=42.0, c_chitosan_current=18.0,
        )
        assert r.alpha_nominal < 1.0

    def test_alpha_greater_than_1_when_G_target_above_current(self):
        r = alpha_for_target_G_DN(
            G_DN_current=10_000.0, target_G_DN=40_000.0,
            c_agarose_current=42.0, c_chitosan_current=18.0,
        )
        assert r.alpha_nominal > 1.0
        # G ~ c^2 -> alpha = sqrt(G_target/G_current) = sqrt(4) = 2.
        assert math.isclose(r.alpha_nominal, 2.0, rel_tol=1e-9)

    def test_alpha_applies_to_both_polymers(self):
        r = alpha_for_target_G_DN(
            G_DN_current=10_000.0, target_G_DN=2_500.0,
            c_agarose_current=42.0, c_chitosan_current=18.0,
        )
        assert math.isclose(r.new_c_agarose, r.alpha_nominal * 42.0, rel_tol=1e-9)
        assert math.isclose(r.new_c_chitosan, r.alpha_nominal * 18.0, rel_tol=1e-9)

    def test_under_gelation_flag(self):
        r = alpha_for_target_G_DN(
            G_DN_current=1_000_000.0, target_G_DN=100.0,
            c_agarose_current=42.0, c_chitosan_current=18.0,
        )
        # alpha ~ sqrt(1e-4) = 1e-2 -> c_agarose ~ 0.4 kg/m^3
        assert r.limited_by == "under_gelation"


class TestSuggestionModules:
    def _ctx(self, **overrides):
        defaults = dict(
            family="agarose_chitosan",
            d32_actual=2e-6, d50_actual=3e-6,
            pore_actual=100e-9, l2_mode="ch_2d",
            cooling_rate_effective=0.17, p_final=0.25,
            G_DN_actual=10_000.0, target_d32=2e-6, target_pore=100e-9, target_G=10_000.0,
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

    def test_crosslinker_module_only_fires_when_G_below_target(self):
        # G way below target
        s_low = xl_mod.generate(self._ctx(G_DN_actual=1_000.0))
        assert s_low is not None and s_low.key == "increase_crosslinker"
        # G way above target
        s_high = xl_mod.generate(self._ctx(G_DN_actual=1_000_000.0))
        assert s_high is None

    def test_polymer_module_only_fires_when_G_above_target(self):
        s_high = pol_mod.generate(self._ctx(G_DN_actual=1_000_000.0))
        assert s_high is not None and s_high.key == "reduce_polymer"
        s_low = pol_mod.generate(self._ctx(G_DN_actual=1_000.0))
        assert s_low is None

    def test_crosslinker_target_unit(self):
        t = xl_mod.derive_target(self._ctx(G_DN_actual=1_000.0))
        assert "mol/m" in t.unit

    def test_polymer_target_unit(self):
        t = pol_mod.derive_target(self._ctx(G_DN_actual=1_000_000.0))
        assert "multiplier" in t.unit
