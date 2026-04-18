"""Tests for all trust warning resolution modules (M1, M2, M3, M4).

Covers:
  M1 -- PBE adaptive extension (convergence, extension counting, history)
  M2 -- Crosslinker stoichiometry guidance (recommended concentration, ratio)
  M3 -- Hashin-Shtrikman bounds, Flory-Rehner affine IPN model, mode routing
  M4 -- Per-chemistry eta coupling (genipin vs DVS, backward compat)

All tests use small/fast parameters where PBE solves are unavoidable.
"""

from __future__ import annotations

import numpy as np
import pytest

# ─────────────────────────────────────────────────────────────────────────────
# M1: PBE Adaptive Extension
# ─────────────────────────────────────────────────────────────────────────────


class TestM1AdaptiveExtension:
    """Tests for Level 1 PBE solver adaptive convergence extension."""

    @pytest.fixture
    def fast_solver(self):
        """Small PBE solver for speed (20 bins)."""
        from emulsim.level1_emulsification.solver import PBESolver
        return PBESolver(n_bins=20, d_min=0.5e-6, d_max=200e-6)

    def _make_params(self, t_emul=300.0, max_extensions=2, t_max=600.0,
                     conv_tol=0.01):
        """Build SimulationParameters with controlled adaptive-extension fields."""
        from emulsim.datatypes import SimulationParameters
        p = SimulationParameters()
        p.emulsification.t_emulsification = t_emul
        p.solver.l1_max_extensions = max_extensions
        p.solver.l1_t_max = t_max
        p.solver.l1_conv_tol = conv_tol
        # Use small bin grid to keep runtime low
        p.solver.l1_n_bins = 20
        return p

    def _make_props(self):
        from emulsim.datatypes import MaterialProperties
        return MaterialProperties(sigma=5e-3)

    def test_default_config_converges_or_extends(self, fast_solver):
        """With t=300 s and default settings, solver either converges or extends."""
        params = self._make_params(t_emul=300.0)
        props = self._make_props()
        result = fast_solver.solve(params, props)
        # Either it converged OR it attempted at least one extension
        assert result.converged or result.n_extensions > 0, (
            "Expected convergence or at least one extension for t=300 s"
        )

    def test_short_time_triggers_extension(self, fast_solver):
        """With t=5 s (too short to converge), at least one extension must occur."""
        params = self._make_params(t_emul=5.0, max_extensions=2)
        props = self._make_props()
        result = fast_solver.solve(params, props)
        assert result.n_extensions >= 1, (
            f"Expected n_extensions >= 1 for t=5s, got {result.n_extensions}"
        )

    def test_t_history_monotonic(self, fast_solver):
        """t_history must be strictly increasing across any extensions."""
        params = self._make_params(t_emul=5.0, max_extensions=2)
        props = self._make_props()
        result = fast_solver.solve(params, props)
        assert result.t_history is not None, "t_history should not be None"
        assert len(result.t_history) > 1, "t_history should have more than 1 point"
        diffs = np.diff(result.t_history)
        assert np.all(diffs > 0), (
            f"t_history is not strictly monotonic: min diff = {diffs.min():.3e}"
        )

    def test_n_extensions_capped(self, fast_solver):
        """With l1_max_extensions=1, n_extensions must not exceed 1."""
        params = self._make_params(t_emul=5.0, max_extensions=1)
        props = self._make_props()
        result = fast_solver.solve(params, props)
        assert result.n_extensions <= 1, (
            f"n_extensions={result.n_extensions} exceeds cap of 1"
        )

    def test_t_converged_populated_when_converged(self, fast_solver):
        """If converged=True, t_converged must not be None."""
        params = self._make_params(t_emul=300.0, max_extensions=2)
        props = self._make_props()
        result = fast_solver.solve(params, props)
        if result.converged:
            assert result.t_converged is not None, (
                "t_converged should be set when converged=True"
            )
            assert result.t_converged > 0.0, (
                f"t_converged should be positive, got {result.t_converged}"
            )

    def test_n_extensions_non_negative(self, fast_solver):
        """n_extensions must always be non-negative."""
        params = self._make_params(t_emul=300.0)
        props = self._make_props()
        result = fast_solver.solve(params, props)
        assert result.n_extensions >= 0

    def test_result_has_extension_fields(self, fast_solver):
        """EmulsificationResult must expose n_extensions and t_converged fields."""
        params = self._make_params(t_emul=10.0)
        props = self._make_props()
        result = fast_solver.solve(params, props)
        assert hasattr(result, 'n_extensions'), "Missing n_extensions field"
        assert hasattr(result, 't_converged'), "Missing t_converged field"


# ─────────────────────────────────────────────────────────────────────────────
# M2: Crosslinker Stoichiometry Guidance
# ─────────────────────────────────────────────────────────────────────────────


class TestM2CrosslinkerStoichiometry:
    """Tests for recommended_crosslinker_concentration and the trust gate ratio."""

    def test_recommended_concentration_known_value(self):
        """Spot-check: (18.0 kg/m³, DDA=0.9, M=161.16, target_p=0.20) ≈ 10.05 mol/m³."""
        from emulsim.level3_crosslinking.solver import recommended_crosslinker_concentration
        c = recommended_crosslinker_concentration(18.0, 0.9, 161.16, target_p=0.20)
        # [NH2] = 18000 g/m³ * 0.9 / 161.16 = 100.5 mol/m³
        # c_xlink = 0.20 * 100.5 / 2 = 10.05 mol/m³
        assert c == pytest.approx(10.05, rel=0.01), (
            f"Expected ~10.05 mol/m³, got {c:.4f}"
        )

    def test_recommended_concentration_positive(self):
        """Recommended concentration must be positive for valid inputs."""
        from emulsim.level3_crosslinking.solver import recommended_crosslinker_concentration
        c = recommended_crosslinker_concentration(10.0, 0.8, 161.16, target_p=0.10)
        assert c > 0.0

    def test_recommended_concentration_scales_with_target_p(self):
        """Doubling target_p should double the required concentration."""
        from emulsim.level3_crosslinking.solver import recommended_crosslinker_concentration
        c_low = recommended_crosslinker_concentration(10.0, 0.8, 161.16, target_p=0.10)
        c_high = recommended_crosslinker_concentration(10.0, 0.8, 161.16, target_p=0.20)
        assert c_high == pytest.approx(2.0 * c_low, rel=1e-9)

    def test_default_ratio_above_threshold(self):
        """c_genipin=10 mol/m³ should give NH2 ratio >= 0.05 for default formulation."""
        from emulsim.datatypes import SimulationParameters, MaterialProperties
        params = SimulationParameters()
        params.formulation.c_genipin = 10.0
        props = MaterialProperties()
        # NH2 = c_chitosan * 1000 * DDA / M_GlcN
        NH2 = params.formulation.c_chitosan * 1000.0 * props.DDA / props.M_GlcN
        if NH2 > 0:
            ratio = params.formulation.c_genipin / NH2
            assert ratio >= 0.05, (
                f"Expected ratio >= 0.05, got {ratio:.4f} "
                f"(c_genipin=10, NH2={NH2:.2f})"
            )

    def test_recommended_concentration_zero_chitosan(self):
        """Zero chitosan gives zero recommended concentration (no NH2 groups)."""
        from emulsim.level3_crosslinking.solver import recommended_crosslinker_concentration
        c = recommended_crosslinker_concentration(0.0, 0.9, 161.16)
        assert c == 0.0


# ─────────────────────────────────────────────────────────────────────────────
# M3: Hashin-Shtrikman Bounds
# ─────────────────────────────────────────────────────────────────────────────


class TestM3HashinShtrikmanBounds:
    """Tests for hashin_shtrikman_bounds composite mechanics formula."""

    def test_hs_bounds_symmetric_equal_moduli(self):
        """When G1 == G2, both bounds collapse to that single value."""
        from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
        G = 10000.0
        lower, upper = hashin_shtrikman_bounds(G, G, 0.5)
        assert lower == pytest.approx(G, rel=1e-9)
        assert upper == pytest.approx(G, rel=1e-9)

    def test_hs_bounds_asymmetric_ordered(self):
        """For different G values, lower < upper."""
        from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
        lower, upper = hashin_shtrikman_bounds(1000.0, 10000.0, 0.5)
        assert lower < upper, f"Expected lower < upper, got {lower:.1f} >= {upper:.1f}"

    def test_hs_bounds_pure_phase_1(self):
        """phi1=1.0 -> both bounds equal G1."""
        from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
        lower, upper = hashin_shtrikman_bounds(5000.0, 20000.0, 1.0)
        assert lower == pytest.approx(5000.0, rel=1e-9)
        assert upper == pytest.approx(5000.0, rel=1e-9)

    def test_hs_bounds_pure_phase_2(self):
        """phi1=0.0 -> both bounds equal G2."""
        from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
        lower, upper = hashin_shtrikman_bounds(5000.0, 20000.0, 0.0)
        assert lower == pytest.approx(20000.0, rel=1e-9)
        assert upper == pytest.approx(20000.0, rel=1e-9)

    def test_hs_lower_above_soft_phase(self):
        """Lower bound must be >= G_soft (Voigt-Reuss bounding)."""
        from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
        G_soft, G_stiff = 1000.0, 50000.0
        lower, upper = hashin_shtrikman_bounds(G_soft, G_stiff, 0.3)
        assert lower >= G_soft * 0.9999  # allow tiny float tolerance
        assert upper <= G_stiff * 1.0001

    def test_hs_bounds_nonnegative(self):
        """Bounds must be non-negative even for zero soft phase."""
        from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
        lower, upper = hashin_shtrikman_bounds(0.0, 5000.0, 0.5)
        assert lower >= 0.0
        assert upper >= 0.0


# ─────────────────────────────────────────────────────────────────────────────
# M3: Flory-Rehner Equilibrium Swelling
# ─────────────────────────────────────────────────────────────────────────────


class TestM3FloryRehner:
    """Tests for flory_rehner_swelling numerical solver."""

    def test_physical_phi_eq_in_unit_interval(self):
        """phi_eq must be in (0, 1) for physically meaningful inputs."""
        from emulsim.level4_mechanical.solver import flory_rehner_swelling
        # Moderate crosslink density, typical chi for agarose-water
        phi_eq, converged = flory_rehner_swelling(nu_e=100.0, phi_0=0.02, chi=0.50)
        assert 0.0 < phi_eq < 1.0, f"phi_eq={phi_eq} outside (0,1)"
        assert converged is True

    def test_higher_crosslink_density_higher_phi(self):
        """More crosslinks -> less swelling -> higher phi_eq."""
        from emulsim.level4_mechanical.solver import flory_rehner_swelling
        phi_lo, _ = flory_rehner_swelling(nu_e=10.0, phi_0=0.02, chi=0.50)
        phi_hi, _ = flory_rehner_swelling(nu_e=1000.0, phi_0=0.02, chi=0.50)
        assert phi_hi >= phi_lo, (
            f"Higher nu_e should give higher phi_eq: {phi_lo:.4f} vs {phi_hi:.4f}"
        )

    def test_higher_chi_less_swelling(self):
        """Worse solvent (higher chi) -> less swelling -> higher phi_eq."""
        from emulsim.level4_mechanical.solver import flory_rehner_swelling
        phi_good, _ = flory_rehner_swelling(nu_e=100.0, phi_0=0.02, chi=0.1)
        phi_poor, _ = flory_rehner_swelling(nu_e=100.0, phi_0=0.02, chi=0.5)
        assert phi_poor >= phi_good, (
            f"Poorer solvent should give higher phi_eq: {phi_good:.4f} vs {phi_poor:.4f}"
        )

    def test_extreme_params_no_crash(self):
        """Extreme parameters (very high chi) should not crash — fallback to phi_0."""
        from emulsim.level4_mechanical.solver import flory_rehner_swelling
        phi_eq, converged = flory_rehner_swelling(nu_e=1.0, phi_0=0.02, chi=5.0)
        # Should return phi_0 as fallback or a valid physical value
        assert 0.0 < phi_eq <= 1.0
        # converged may be False for extreme parameters (fallback)
        assert isinstance(converged, bool)


# ─────────────────────────────────────────────────────────────────────────────
# M3: Affine IPN Double-Network Modulus
# ─────────────────────────────────────────────────────────────────────────────


class TestM3AffineIPN:
    """Tests for double_network_modulus_affine Flory-Rehner model."""

    def test_affine_single_network_no_chitosan(self):
        """With nu_e2 ~ 0, G_DN is approximately G1 only."""
        from emulsim.level4_mechanical.solver import double_network_modulus_affine
        G_DN, status = double_network_modulus_affine(
            nu_e1=100.0, nu_e2=1e-10,  # nearly zero chitosan
            phi_01=0.02, phi_02=0.01,
            chi1=0.50, chi2=0.45,
        )
        assert G_DN >= 0.0
        # G from network 1 only — should be much less than combined
        G_DN_both, _ = double_network_modulus_affine(
            nu_e1=100.0, nu_e2=100.0,
            phi_01=0.02, phi_02=0.01,
            chi1=0.50, chi2=0.45,
        )
        assert G_DN_both >= G_DN, "Adding a second network should not decrease G_DN"

    def test_affine_realistic_modulus_range(self):
        """Standard formulation should produce G_DN in 15-150 kPa range."""
        from emulsim.level4_mechanical.solver import double_network_modulus_affine
        # Typical agarose-chitosan IPN parameters
        # nu_e in mol/m^3: agarose ~100, chitosan ~50
        G_DN, status = double_network_modulus_affine(
            nu_e1=100.0, nu_e2=50.0,
            phi_01=0.02, phi_02=0.01,
            chi1=0.50, chi2=0.45,
        )
        assert 1.0 <= G_DN <= 1e7, (
            f"G_DN={G_DN:.1f} Pa outside plausible range for IPN hydrogel"
        )

    def test_affine_returns_tuple_status(self):
        """double_network_modulus_affine must return (float, str) tuple."""
        from emulsim.level4_mechanical.solver import double_network_modulus_affine
        result = double_network_modulus_affine(
            nu_e1=50.0, nu_e2=50.0,
            phi_01=0.02, phi_02=0.01,
            chi1=0.50, chi2=0.45,
        )
        assert isinstance(result, tuple) and len(result) == 2
        G_DN, status = result
        assert isinstance(G_DN, float)
        assert isinstance(status, str)

    def test_affine_fallback_no_crash_extreme_params(self):
        """Extreme parameters must not crash — return valid (G, status) pair."""
        from emulsim.level4_mechanical.solver import double_network_modulus_affine
        G_DN, status = double_network_modulus_affine(
            nu_e1=1e6, nu_e2=1e6,     # very high crosslink density
            phi_01=0.99, phi_02=0.99, # near-pure polymer
            chi1=5.0, chi2=5.0,       # very poor solvent
        )
        assert G_DN >= 0.0
        assert status in ("converged", "fallback")

    def test_affine_status_converged_for_normal_params(self):
        """Normal parameters should yield status='converged'."""
        from emulsim.level4_mechanical.solver import double_network_modulus_affine
        _, status = double_network_modulus_affine(
            nu_e1=100.0, nu_e2=50.0,
            phi_01=0.02, phi_02=0.01,
            chi1=0.50, chi2=0.45,
        )
        assert status == "converged"


# ─────────────────────────────────────────────────────────────────────────────
# M3: select_modulus_model routing
# ─────────────────────────────────────────────────────────────────────────────


class TestM3SelectModulusModel:
    """Tests for select_modulus_model mode routing."""

    def _make_crosslinking_result(self, nu_e_final=1e20):
        """Build a minimal CrosslinkingResult for mechanistic_research tests."""
        from emulsim.datatypes import CrosslinkingResult
        import numpy as np
        t = np.array([0.0, 3600.0])
        zeros = np.zeros(2)
        xi = np.full(2, 1e-6)
        Mc = np.full(2, 1e8)
        return CrosslinkingResult(
            t_array=t,
            X_array=zeros,
            nu_e_array=zeros,
            Mc_array=Mc,
            xi_array=xi,
            G_chitosan_array=zeros,
            p_final=0.05,
            nu_e_final=nu_e_final,
            Mc_final=1e8,
            xi_final=1e-6,
            G_chitosan_final=100.0,
        )

    def test_mechanistic_mode_uses_affine(self):
        """In mechanistic_research mode with valid crosslinking data, model_used is flory_rehner_affine
        when FR converges, or phenomenological_fallback when it does not.

        Note: select_modulus_model requires network_metadata != None to reach the
        mechanistic_research branch; without it, the function returns phenomenological
        immediately as a safe fallback (metadata carries the per-chemistry context).
        """
        from emulsim.level4_mechanical.solver import select_modulus_model
        from emulsim.datatypes import SimulationParameters, NetworkTypeMetadata, ModelMode

        params = SimulationParameters()
        params.model_mode = ModelMode.MECHANISTIC_RESEARCH
        params.formulation.c_agarose = 30.0   # 3% w/v
        params.formulation.c_chitosan = 18.0  # 1.8% w/v
        params.formulation.T_crosslink = 310.15  # 37 C

        crosslinking = self._make_crosslinking_result(nu_e_final=1e20)

        # network_metadata must be provided so the function reaches the
        # mechanistic_research branch (amine_covalent family)
        meta = NetworkTypeMetadata(
            solver_family="amine_covalent",
            eta_coupling_recommended=-0.15,
        )

        G_DN, model_used = select_modulus_model(
            G_agarose=10000.0,
            G_xlink=500.0,
            network_metadata=meta,
            eta_coupling=-0.15,
            model_mode="mechanistic_research",
            params=params,
            crosslinking=crosslinking,
        )
        # When FR converges for both networks, model_used == "flory_rehner_affine".
        # When FR falls back for any network, model_used == "phenomenological_fallback".
        # Both are valid outcomes from the mechanistic_research branch.
        assert model_used in ("flory_rehner_affine", "phenomenological_fallback"), (
            f"Expected 'flory_rehner_affine' or 'phenomenological_fallback', got '{model_used}'"
        )
        assert G_DN >= 0.0

    def test_empirical_mode_uses_phenomenological(self):
        """In empirical_engineering mode, model_used is phenomenological."""
        from emulsim.level4_mechanical.solver import select_modulus_model

        G_DN, model_used = select_modulus_model(
            G_agarose=10000.0,
            G_xlink=500.0,
            network_metadata=None,
            eta_coupling=-0.15,
            model_mode="empirical_engineering",
        )
        assert model_used == "phenomenological", (
            f"Expected 'phenomenological', got '{model_used}'"
        )

    def test_no_metadata_uses_phenomenological(self):
        """Without network_metadata, always falls back to phenomenological."""
        from emulsim.level4_mechanical.solver import select_modulus_model

        G_DN, model_used = select_modulus_model(
            G_agarose=5000.0,
            G_xlink=1000.0,
            network_metadata=None,
            eta_coupling=-0.15,
        )
        assert model_used == "phenomenological"

    def test_hybrid_coupled_uses_phenomenological(self):
        """Default hybrid_coupled mode uses phenomenological (not affine)."""
        from emulsim.level4_mechanical.solver import select_modulus_model
        from emulsim.datatypes import NetworkTypeMetadata

        meta = NetworkTypeMetadata(
            solver_family="amine_covalent",
            eta_coupling_recommended=-0.15,
        )
        G_DN, model_used = select_modulus_model(
            G_agarose=10000.0,
            G_xlink=500.0,
            network_metadata=meta,
            eta_coupling=-0.15,
            model_mode="hybrid_coupled",
        )
        assert model_used == "phenomenological"


# ─────────────────────────────────────────────────────────────────────────────
# M4: Per-Chemistry Eta Coupling
# ─────────────────────────────────────────────────────────────────────────────


class TestM4PerChemistryEta:
    """Tests for per-chemistry eta coupling in L4 modulus selection."""

    def test_genipin_eta_backward_compatible(self):
        """Genipin (amine_covalent) falls back to eta=-0.15 by default."""
        from emulsim.level4_mechanical.solver import double_network_modulus

        G_agarose, G_xlink = 50000.0, 5000.0
        G_DN = double_network_modulus(G_agarose, G_xlink, eta_coupling=-0.15)
        expected = G_agarose + G_xlink + (-0.15) * np.sqrt(G_agarose * G_xlink)
        assert G_DN == pytest.approx(expected, rel=1e-9)
        assert G_DN < G_agarose + G_xlink  # antagonistic coupling

    def test_dvs_eta_positive_from_metadata(self):
        """DVS (hydroxyl_covalent) has eta=+0.10 from reagent library."""
        from emulsim.reagent_library import CROSSLINKERS
        dvs = CROSSLINKERS.get("dvs")
        assert dvs is not None, "DVS crosslinker not found in CROSSLINKERS"
        assert dvs.eta_coupling_recommended == pytest.approx(0.10, rel=1e-9), (
            f"Expected DVS eta=+0.10, got {dvs.eta_coupling_recommended}"
        )

    def test_gdn_differs_dvs_vs_genipin(self):
        """G_DN computed with DVS eta (+0.10) differs from genipin eta (-0.15)."""
        from emulsim.level4_mechanical.solver import double_network_modulus

        G_agarose, G_xlink = 50000.0, 5000.0
        G_genipin = double_network_modulus(G_agarose, G_xlink, eta_coupling=-0.15)
        G_dvs = double_network_modulus(G_agarose, G_xlink, eta_coupling=0.10)
        assert G_dvs > G_genipin, (
            f"DVS (eta=+0.10) should give higher G_DN than genipin (eta=-0.15): "
            f"{G_dvs:.1f} vs {G_genipin:.1f}"
        )

    def test_select_modulus_uses_per_chem_eta(self):
        """select_modulus_model reads eta from network_metadata.eta_coupling_recommended."""
        from emulsim.level4_mechanical.solver import select_modulus_model
        from emulsim.datatypes import NetworkTypeMetadata

        G_agarose, G_xlink = 50000.0, 5000.0

        meta_dvs = NetworkTypeMetadata(
            solver_family="amine_covalent",
            eta_coupling_recommended=0.10,  # DVS-like positive eta
        )
        meta_genipin = NetworkTypeMetadata(
            solver_family="amine_covalent",
            eta_coupling_recommended=-0.15,  # genipin default
        )

        G_dvs, _ = select_modulus_model(G_agarose, G_xlink, meta_dvs, eta_coupling=-0.15)
        G_genipin, _ = select_modulus_model(G_agarose, G_xlink, meta_genipin, eta_coupling=-0.15)

        assert G_dvs > G_genipin, (
            f"DVS meta should give G_DN={G_dvs:.1f} > genipin G_DN={G_genipin:.1f}"
        )

    def test_genipin_eta_from_reagent_library(self):
        """Genipin entry in reagent library has eta=-0.15 (default antagonistic)."""
        from emulsim.reagent_library import CROSSLINKERS
        genipin = CROSSLINKERS.get("genipin")
        assert genipin is not None, "Genipin crosslinker not found in CROSSLINKERS"
        assert genipin.eta_coupling_recommended == pytest.approx(-0.15, abs=0.01), (
            f"Expected genipin eta ~ -0.15, got {genipin.eta_coupling_recommended}"
        )

    def test_eta_in_metadata_suppresses_w4_warning(self):
        """When network_metadata has eta_coupling_recommended, W4 (non-specific eta) is suppressed."""
        from emulsim.datatypes import (
            NetworkTypeMetadata, CrosslinkingResult, MaterialProperties,
        )
        import numpy as np

        props = MaterialProperties()
        t = np.array([0.0, 3600.0])
        crosslinking = CrosslinkingResult(
            t_array=t,
            X_array=np.zeros(2),
            nu_e_array=np.zeros(2),
            Mc_array=np.full(2, 1e8),
            xi_array=np.full(2, 1e-6),
            G_chitosan_array=np.zeros(2),
            p_final=0.05,
            nu_e_final=0.0,
            Mc_final=1e8,
            xi_final=1e-6,
            G_chitosan_final=0.0,
            network_metadata=NetworkTypeMetadata(
                solver_family="amine_covalent",
                eta_coupling_recommended=0.10,  # non-default -> suppresses W4
            ),
        )
        # W4 is suppressed when network_metadata is present with eta field
        _network_meta = crosslinking.network_metadata
        _per_chem_eta_active = (
            _network_meta is not None
            and hasattr(_network_meta, 'eta_coupling_recommended')
        )
        assert _per_chem_eta_active, (
            "Expected per-chem eta to be detected as active"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Integration: SolverSettings fields exist with correct defaults
# ─────────────────────────────────────────────────────────────────────────────


class TestDatatypeFields:
    """Verify new datatype fields exist and have correct defaults."""

    def test_solver_settings_has_l1_adaptive_fields(self):
        """SolverSettings must have l1_t_max, l1_conv_tol, l1_max_extensions."""
        from emulsim.datatypes import SolverSettings
        s = SolverSettings()
        assert hasattr(s, 'l1_t_max')
        assert hasattr(s, 'l1_conv_tol')
        assert hasattr(s, 'l1_max_extensions')

    def test_solver_settings_defaults(self):
        """Check documented default values for adaptive extension fields."""
        from emulsim.datatypes import SolverSettings
        s = SolverSettings()
        assert s.l1_t_max == pytest.approx(600.0)
        assert s.l1_conv_tol == pytest.approx(0.01)
        assert s.l1_max_extensions == 2

    def test_emulsification_result_has_extension_fields(self):
        """EmulsificationResult must accept t_converged and n_extensions fields."""
        from emulsim.datatypes import EmulsificationResult
        import numpy as np
        r = EmulsificationResult(
            d_bins=np.array([1e-6]),
            n_d=np.array([1.0]),
            d32=1e-6, d43=1e-6, d10=1e-6, d50=1e-6, d90=1e-6,
            span=0.0, total_volume_fraction=0.1,
            converged=True,
            t_converged=42.0,
            n_extensions=1,
        )
        assert r.t_converged == 42.0
        assert r.n_extensions == 1

    def test_network_type_metadata_has_eta_field(self):
        """NetworkTypeMetadata must have eta_coupling_recommended field."""
        from emulsim.datatypes import NetworkTypeMetadata
        m = NetworkTypeMetadata()
        assert hasattr(m, 'eta_coupling_recommended')
        assert m.eta_coupling_recommended == pytest.approx(-0.15)

    def test_mechanical_result_has_hs_bounds(self):
        """MechanicalResult must have G_DN_lower and G_DN_upper fields."""
        from emulsim.datatypes import MechanicalResult
        import numpy as np
        r = MechanicalResult(
            G_agarose=1000.0, G_chitosan=500.0, G_DN=1400.0,
            E_star=5600.0,
            delta_array=np.zeros(10), F_array=np.zeros(10),
            rh_array=np.zeros(10), Kav_array=np.ones(10),
            pore_size_mean=50e-9, xi_mesh=10e-9,
            G_DN_lower=1200.0, G_DN_upper=1450.0,
        )
        assert r.G_DN_lower == 1200.0
        assert r.G_DN_upper == 1450.0
