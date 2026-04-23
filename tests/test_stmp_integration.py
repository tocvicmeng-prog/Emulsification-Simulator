"""STMP (Sodium Trimetaphosphate) integration tests.

Covers:
- M-E: primary crosslinker profile presence and parameter sanity
- M-F: end-to-end L3 solver smoke test with stmp
- M-G: secondary Module-2 profile presence and orchestrator validation

Scientific basis: SA-EMULSIM-XL-002 Rev 0.1; implementation plan
JCP-EMULSIM-STMP-001 Rev 0.
"""

from __future__ import annotations

import numpy as np

from emulsim.datatypes import MaterialProperties, SimulationParameters
from emulsim.level3_crosslinking.solver import solve_crosslinking
from emulsim.module2_functionalization.acs import ACSSiteType
from emulsim.module2_functionalization.reagent_profiles import REAGENT_PROFILES
from emulsim.reagent_library import CROSSLINKERS


class TestSTMPPrimaryProfile:
    """M-E: CROSSLINKERS['stmp'] presence and sanity."""

    def test_stmp_registered(self):
        assert "stmp" in CROSSLINKERS

    def test_stmp_cas_is_trimetaphosphate_not_tripolyphosphate(self):
        # 7785-84-4 is STMP (Na3P3O9); 7758-29-4 would be STPP — different compound.
        assert CROSSLINKERS["stmp"].cas == "7785-84-4"
        assert CROSSLINKERS["tpp"].cas == "7758-29-4"

    def test_stmp_routes_through_hydroxyl_mechanism(self):
        # mechanism="hydroxyl" is the whitelist key that dispatches the
        # hydroxyl second-order solver. Changing this breaks routing.
        assert CROSSLINKERS["stmp"].mechanism == "hydroxyl"

    def test_stmp_kinetics_parameters_physical(self):
        xl = CROSSLINKERS["stmp"]
        assert xl.kinetics_model == "second_order"
        assert xl.k_xlink_0 > 0
        assert 30e3 < xl.E_a_xlink < 150e3  # reasonable Arrhenius barrier
        assert 0.0 < xl.f_bridge < 1.0
        assert xl.T_crosslink_default > 273.15
        assert xl.t_crosslink_default > 0

    def test_stmp_network_metadata(self):
        xl = CROSSLINKERS["stmp"]
        assert xl.solver_family == "hydroxyl_covalent"
        assert xl.network_target == "mixed"
        assert xl.reversible is False


class TestSTMPLevel3Smoke:
    """M-F: end-to-end L3 solve via the dispatch path."""

    def _params_props(self):
        params = SimulationParameters()
        # Use STMP defaults so the scenario is physically consistent.
        xl = CROSSLINKERS["stmp"]
        params.formulation.T_crosslink = xl.T_crosslink_default
        params.formulation.t_crosslink = xl.t_crosslink_default
        params.formulation.c_genipin = 65.0  # ~2 % w/v STMP in mM
        props = MaterialProperties()
        return params, props

    def test_stmp_solver_completes(self):
        params, props = self._params_props()
        result = solve_crosslinking(params, props, crosslinker_key="stmp")
        assert len(result.t_array) > 0
        assert 0.0 <= result.p_final <= 1.0

    def test_stmp_produces_nonzero_conversion(self):
        params, props = self._params_props()
        result = solve_crosslinking(params, props, crosslinker_key="stmp")
        assert result.p_final > 0.005, (
            f"STMP default conditions should give measurable conversion, "
            f"got p_final={result.p_final}"
        )

    def test_stmp_routes_to_hydroxyl_covalent_mixed_network(self):
        params, props = self._params_props()
        result = solve_crosslinking(params, props, crosslinker_key="stmp")
        assert result.network_metadata is not None
        assert result.network_metadata.solver_family == "hydroxyl_covalent"
        assert result.network_metadata.network_target == "mixed"

    def test_stmp_conversion_increases_monotonically(self):
        params, props = self._params_props()
        result = solve_crosslinking(params, props, crosslinker_key="stmp")
        assert np.all(np.diff(result.X_array) >= -1e-12)

    def test_stmp_handles_small_and_medium_beads(self):
        """Representative EmulSim bead radii — neither should raise."""
        params, props = self._params_props()
        for R in (50e-6, 250e-6):
            result = solve_crosslinking(
                params, props, crosslinker_key="stmp", R_droplet=R,
            )
            assert result.p_final >= 0.0


class TestSTMPSecondaryProfile:
    """M-G: REAGENT_PROFILES['stmp_secondary'] presence and routing."""

    def test_stmp_secondary_registered(self):
        assert "stmp_secondary" in REAGENT_PROFILES

    def test_stmp_secondary_is_crosslinking_type(self):
        # Orchestrator rule 4: SECONDARY_CROSSLINKING step only accepts
        # reagents whose reaction_type == "crosslinking".
        rp = REAGENT_PROFILES["stmp_secondary"]
        assert rp.reaction_type == "crosslinking"

    def test_stmp_secondary_targets_hydroxyl(self):
        # First HYDROXYL-targeted secondary crosslinker in the library.
        rp = REAGENT_PROFILES["stmp_secondary"]
        assert rp.target_acs == ACSSiteType.HYDROXYL
        assert rp.product_acs is None  # terminal crosslink, no downstream ACS

    def test_stmp_secondary_cas_matches_primary(self):
        # The primary and secondary entries must refer to the same compound.
        assert (
            REAGENT_PROFILES["stmp_secondary"].cas == CROSSLINKERS["stmp"].cas
        )

    def test_stmp_secondary_kinetic_window(self):
        rp = REAGENT_PROFILES["stmp_secondary"]
        assert rp.k_forward > 0
        assert 30e3 < rp.E_a < 150e3
        assert 10.0 <= rp.ph_optimum <= 12.0  # alkaline trigger
        assert rp.temperature_min < rp.temperature_default < rp.temperature_max

    def test_stmp_secondary_chemistry_class_new_string(self):
        # chemistry_class is free-form; STMP introduces "phosphorylation_alkaline".
        rp = REAGENT_PROFILES["stmp_secondary"]
        assert rp.chemistry_class == "phosphorylation_alkaline"
