"""Tests for Level 3 crosslinking kinetics."""

import numpy as np

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


class TestReactionDiffusionPDE:
    """Tests for _solve_reaction_diffusion (diffusion-limited regime, Phi > 1)."""

    def test_pde_solver_runs_and_converges(self):
        """Large droplet + fast kinetics triggers Phi > 1 PDE path."""
        from emulsim.level3_crosslinking.solver import (
            _solve_reaction_diffusion,
            compute_thiele_modulus,
            available_amine_concentration,
        )
        params = SimulationParameters()
        # Large droplet radius to push Thiele > 1
        params.formulation.c_genipin = 50.0   # high crosslinker
        params.formulation.t_crosslink = 3600.0  # 1 hour
        props = MaterialProperties()

        R_droplet = 500e-6  # 500 um — large enough for Phi > 1
        D_xlink = 1e-10

        # Verify Phi > 1
        from emulsim.level3_crosslinking.solver import arrhenius_rate_constant
        k_rate = arrhenius_rate_constant(
            params.formulation.T_crosslink, props.k_xlink_0, props.E_a_xlink
        )
        NH2 = available_amine_concentration(
            params.formulation.c_chitosan, props.DDA, props.M_GlcN,
        )
        Phi = compute_thiele_modulus(R_droplet, k_rate, NH2, D_xlink)
        assert Phi > 1.0, f"Test setup error: Phi={Phi:.2f} should be > 1"

        # Run PDE solver
        result, Phi_out = _solve_reaction_diffusion(params, props, R_droplet, D_xlink)

        # Basic validity checks
        assert result.p_final >= 0.0
        assert result.p_final <= 1.0
        assert len(result.t_array) > 0
        assert result.G_chitosan_final >= 0.0
        assert result.xi_final > 0.0

    def test_pde_no_negative_concentrations(self):
        """BDF solver should not produce negative concentrations after clamp fix."""
        from emulsim.level3_crosslinking.solver import _solve_reaction_diffusion
        params = SimulationParameters()
        params.formulation.c_genipin = 100.0  # very high — stress test
        params.formulation.t_crosslink = 7200.0
        props = MaterialProperties()

        R_droplet = 300e-6
        result, _ = _solve_reaction_diffusion(params, props, R_droplet, 1e-10)

        # Conversion should be valid
        assert 0.0 <= result.p_final <= 1.0
        # X_array (crosslinker consumed) should be non-negative
        assert np.all(result.X_array >= -1e-15), \
            f"X_array has negative values: min={result.X_array.min()}"

    def test_pde_conversion_valid_range(self):
        """PDE solver with different droplet sizes produces valid results."""
        from emulsim.level3_crosslinking.solver import _solve_reaction_diffusion
        params = SimulationParameters()
        params.formulation.c_genipin = 10.0
        params.formulation.t_crosslink = 14400.0  # 4 hours
        props = MaterialProperties()

        # Two different bead sizes — both should produce valid results
        for R in [200e-6, 500e-6]:
            result, Phi = _solve_reaction_diffusion(params, props, R_droplet=R)
            assert 0.0 <= result.p_final <= 1.0, \
                f"p_final out of range for R={R*1e6:.0f} um: {result.p_final}"
            assert result.G_chitosan_final >= 0.0
            assert result.xi_final > 0.0
            assert Phi > 0.0


# ─── Node 9 (v6.1, F9): EDC/NHS chemistry applicability gate ──────────────


class TestNode9EDCNHSApplicability:
    """Verifies that EDC/NHS dispatch downgrades the L3 manifest tier.

    Acceptance for Node 9 (F9):
      - Genipin (matched chemistry, native amines available) -> SEMI_QUANTITATIVE
        and ``approximate_fallback`` False.
      - EDC/NHS (michaelis_menten kinetics_model, no dedicated solver, no M2
        carboxyl context) -> QUALITATIVE_TREND and ``approximate_fallback`` True;
        manifest assumptions explain the fallback so the user/UI/optimizer
        can react.
    """

    def test_genipin_is_semi_quantitative(self):
        from emulsim.datatypes import ModelEvidenceTier
        params = SimulationParameters()
        props = MaterialProperties()
        result = solve_crosslinking(params, props, crosslinker_key="genipin")
        assert result.model_manifest is not None
        assert result.model_manifest.evidence_tier == ModelEvidenceTier.SEMI_QUANTITATIVE
        assert result.model_manifest.diagnostics.get("approximate_fallback") is False

    def test_edc_nhs_downgraded_to_qualitative(self):
        from emulsim.datatypes import ModelEvidenceTier
        params = SimulationParameters()
        props = MaterialProperties()
        result = solve_crosslinking(params, props, crosslinker_key="edc_nhs")
        assert result.model_manifest is not None
        # F9 fix: EDC/NHS without COOH context -> QUALITATIVE_TREND
        assert result.model_manifest.evidence_tier == ModelEvidenceTier.QUALITATIVE_TREND
        assert result.model_manifest.diagnostics.get("approximate_fallback") is True
        # Assumption text must explain the limitation so the trust banner /
        # optimizer Pareto label can carry the warning forward.
        joined_assumptions = " ".join(result.model_manifest.assumptions)
        assert "carboxyl" in joined_assumptions.lower()
        assert "fallback" in joined_assumptions.lower()
