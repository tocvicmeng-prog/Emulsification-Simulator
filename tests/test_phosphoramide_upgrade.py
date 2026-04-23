"""Phosphoramide (NH2 co-reaction) upgrade tests.

Covers the v9.2.2 upgrade of SA-002 QUALITATIVE_TREND → SEMI_QUANTITATIVE:
STMP's chitosan-NH2 side-reaction is now modelled as a parallel
second-order ODE alongside the agarose-OH diester track.

Scientific basis: SA-EMULSIM-XL-002 Rev 0.1 + Seal BL (1996) Biomaterials
17:1869 + Salata et al. (2015) Int. J. Biol. Macromol. 81:1009.
"""

from __future__ import annotations

import numpy as np

from emulsim.datatypes import MaterialProperties, SimulationParameters
from emulsim.level3_crosslinking.solver import solve_crosslinking
from emulsim.reagent_library import CROSSLINKERS, NH2CoReaction


class TestNH2CoReactionDataclass:
    def test_stmp_has_nh2_co_reaction_populated(self):
        stmp = CROSSLINKERS["stmp"]
        assert stmp.nh2_co_reaction is not None
        assert isinstance(stmp.nh2_co_reaction, NH2CoReaction)

    def test_non_stmp_hydroxyl_crosslinkers_unchanged(self):
        # ECH, DVS, citric_acid keep nh2_co_reaction = None (their solver
        # path is untouched by v9.2.2; the physics simplification flagged
        # in _solve_second_order_hydroxyl's docstring is a separate issue).
        for key in ("epichlorohydrin", "dvs", "citric_acid"):
            assert CROSSLINKERS[key].nh2_co_reaction is None, (
                f"{key} should not have an NH2 track populated in v9.2.2"
            )

    def test_stmp_nh2_rate_constants_physical(self):
        nh2 = CROSSLINKERS["stmp"].nh2_co_reaction
        assert nh2 is not None
        # k0 calibrated so k_NH2(60C) gives effective rate ~1/5 of k_OH(60C)*[OH]/[NH2]
        assert 1e2 < nh2.k0_nh2 < 1e5
        # Ea lower than diester's 75 kJ/mol: amine deprotonation barrier
        # is lower at pH 11
        assert 40e3 < nh2.E_a_nh2 < 90e3
        # f_bridge less than diester's 0.45 due to pendant adducts
        assert 0.0 < nh2.f_bridge_nh2 < 0.5
        assert nh2.confidence_tier == "SEMI_QUANTITATIVE"


class TestPhosphoramideSolverIntegration:
    def _run(self, crosslinker_key="stmp", **overrides):
        params = SimulationParameters()
        xl = CROSSLINKERS[crosslinker_key]
        params.formulation.T_crosslink = xl.T_crosslink_default
        params.formulation.t_crosslink = xl.t_crosslink_default
        params.formulation.c_genipin = 65.0  # ~2% STMP in mM
        for k, v in overrides.items():
            setattr(params.formulation, k, v)
        props = MaterialProperties()
        return solve_crosslinking(params, props, crosslinker_key=crosslinker_key)

    def test_stmp_produces_nonzero_phosphoramide_contribution(self):
        result = self._run(crosslinker_key="stmp")
        assert result.G_chit_phosphoramide > 0.0, (
            "STMP NH2 track should produce a non-zero phosphoramide modulus"
        )

    def test_stmp_split_sums_to_total(self):
        result = self._run(crosslinker_key="stmp")
        total = result.G_chit_diester + result.G_chit_phosphoramide
        assert np.isclose(total, result.G_chitosan_final, rtol=1e-9), (
            f"G_chit_diester + G_chit_phosphoramide = {total} should equal "
            f"G_chitosan_final = {result.G_chitosan_final}"
        )

    def test_stmp_diester_dominates_phosphoramide(self):
        """NH2 rate is ~1/10 of OH rate, so diester G >> phosphoramide G."""
        result = self._run(crosslinker_key="stmp")
        assert result.G_chit_diester > result.G_chit_phosphoramide, (
            "OH diester should dominate; NH2 phosphoramide is ~1/10 the rate"
        )

    def test_ech_does_not_populate_split_fields(self):
        # ECH (nh2_co_reaction = None) runs the single-track path; split
        # fields remain at their 0.0 defaults.
        result = self._run(crosslinker_key="epichlorohydrin")
        assert result.G_chit_diester == 0.0
        assert result.G_chit_phosphoramide == 0.0
        assert result.p_final_nh2 == 0.0

    def test_stmp_nh2_conversion_in_valid_range(self):
        result = self._run(crosslinker_key="stmp")
        assert 0.0 <= result.p_final_nh2 <= 1.0

    def test_stmp_total_G_greater_than_single_track_would_be(self):
        """The dual-track result must be >= what a hypothetical single-track
        run would give (adding a parallel contribution can only increase G)."""
        stmp_result = self._run(crosslinker_key="stmp")
        # Single-track baseline: the diester component alone, which is what
        # G_chit_diester records.
        assert stmp_result.G_chitosan_final >= stmp_result.G_chit_diester

    def test_arrays_have_consistent_length(self):
        result = self._run(crosslinker_key="stmp")
        assert len(result.G_chitosan_array) == len(result.t_array)
        assert len(result.X_array) == len(result.t_array)


class TestPhosphoramideKineticsFirstPrinciples:
    """Property tests on the NH2 rate constants (not solver integration)."""

    def test_nh2_rate_increases_with_temperature(self):
        from emulsim.level3_crosslinking.solver import arrhenius_rate_constant
        nh2 = CROSSLINKERS["stmp"].nh2_co_reaction
        assert nh2 is not None
        k_low = arrhenius_rate_constant(298.15, nh2.k0_nh2, nh2.E_a_nh2)
        k_high = arrhenius_rate_constant(333.15, nh2.k0_nh2, nh2.E_a_nh2)
        assert k_high > k_low

    def test_effective_rate_ratio_matches_sa_audit(self):
        """Per-molecule, NH2 is ~2x more nucleophilic than OH (alpha effect).
        But in a typical agarose-chitosan bead [NH2]/[OH] ~ 0.1, so the
        EFFECTIVE rate ratio k_NH2*[NH2] / k_OH*[OH] should be ~0.1-0.5.
        This is the measurable quantity that the dual-track solver uses."""
        from emulsim.level3_crosslinking.solver import (
            arrhenius_rate_constant,
            available_amine_concentration,
            available_hydroxyl_concentration,
        )
        stmp = CROSSLINKERS["stmp"]
        nh2 = stmp.nh2_co_reaction
        assert nh2 is not None
        T = stmp.T_crosslink_default  # 60 °C
        k_oh = arrhenius_rate_constant(T, stmp.k_xlink_0, stmp.E_a_xlink)
        k_nh2 = arrhenius_rate_constant(T, nh2.k0_nh2, nh2.E_a_nh2)
        # Reactive-group pools for a typical 4% agarose + 1.8% chitosan bead,
        # 90% DDA. See SA-002 §2.2.
        OH_pool = available_hydroxyl_concentration(42.0)        # c_agarose in kg/m3
        NH2_pool = available_amine_concentration(18.0, 0.9, 161.16)
        rate_oh = k_oh * OH_pool
        rate_nh2 = k_nh2 * NH2_pool
        effective_ratio = rate_nh2 / rate_oh
        assert 0.05 < effective_ratio < 1.0, (
            f"Effective rate ratio k_NH2·[NH2] / k_OH·[OH] = {effective_ratio:.3f}; "
            f"expected 0.05-1.0 per SA audit (NH2 is ~2x more nucleophilic "
            f"per site but ~10x less concentrated)."
        )
