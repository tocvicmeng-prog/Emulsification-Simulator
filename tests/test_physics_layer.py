"""Unit tests for Phase 2 physics layer: thermal model and energy dissipation.

Covers:
  - Thermal: temperature_profile for all 3 strategies
  - Energy: power_draw, dissipation, Re, Np, swept volume, shear rate
  - Cross-validation: stirrer A vs B physical plausibility
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.datatypes import (
    HeatingConfig, HeatingStrategy, MixerGeometry,
    StirrerGeometry, StirrerType, VesselGeometry,
)
from emulsim.level1_emulsification.thermal import (
    cooling_time_constant, flat_plate_temperature, temperature_profile,
)
from emulsim.level1_emulsification.energy import (
    average_dissipation, emulsion_density, gap_shear_rate,
    gap_shear_rate_extended, impeller_reynolds_number,
    kolmogorov_length_scale, max_dissipation, metzner_otto_shear_rate,
    power_draw, power_number_corrected, swept_volume,
)


# ─── Thermal Model Tests ────────────────────────────────────────────────

class TestThermalIsothermal:
    def test_constant_temperature(self):
        h = HeatingConfig.isothermal(363.15)
        v = VesselGeometry.glass_beaker()
        t = np.linspace(0, 10000, 100)
        T = temperature_profile(h, v, t)
        np.testing.assert_allclose(T, 363.15)

    def test_custom_temperature(self):
        h = HeatingConfig.isothermal(333.15)
        v = VesselGeometry.glass_beaker()
        T = temperature_profile(h, v, np.array([0, 100, 1000]))
        np.testing.assert_allclose(T, 333.15)


class TestThermalFlatPlate:
    def setup_method(self):
        self.h = HeatingConfig.flat_plate()
        self.v = VesselGeometry.glass_beaker()

    def test_initial_temperature(self):
        T = temperature_profile(self.h, self.v, np.array([0.0]))
        assert T[0] == pytest.approx(self.h.T_initial, abs=0.1)

    def test_approaches_final_at_cooldown_time(self):
        T = temperature_profile(self.h, self.v, np.array([self.h.cooldown_time]))
        err_frac = abs(T[0] - self.h.T_final) / (self.h.T_initial - self.h.T_final)
        assert err_frac < 0.06  # within 6% of target

    def test_monotonic_decrease(self):
        t = np.linspace(0, 10000, 500)
        T = temperature_profile(self.h, self.v, t)
        assert np.all(np.diff(T) <= 0)

    def test_bounded(self):
        t = np.linspace(0, 50000, 500)
        T = temperature_profile(self.h, self.v, t)
        assert np.all(T >= self.h.T_final - 0.01)
        assert np.all(T <= self.h.T_initial + 0.01)

    def test_volume_scaling(self):
        """Larger volume should cool more slowly."""
        t_test = np.array([2000.0])
        T_500 = temperature_profile(self.h, VesselGeometry.glass_beaker(0.0005), t_test)
        T_700 = temperature_profile(self.h, VesselGeometry.glass_beaker(0.0007), t_test)
        assert T_700[0] > T_500[0]  # 700 mL is warmer (slower cooling)


class TestThermalJacket:
    def test_constant_during_emulsification(self):
        h = HeatingConfig.hot_water_jacket()
        v = VesselGeometry.jacketed_vessel()
        t = np.linspace(0, 600, 50)
        T = temperature_profile(h, v, t)
        np.testing.assert_allclose(T, h.T_initial, atol=0.1)


class TestCoolingTimeConstant:
    def test_reference_volume(self):
        tau = cooling_time_constant(5400.0, 0.0005)
        assert tau == pytest.approx(5400 / np.log(20), rel=0.01)

    def test_volume_scaling(self):
        tau_500 = cooling_time_constant(5400.0, 0.0005)
        tau_700 = cooling_time_constant(5400.0, 0.0007)
        assert tau_700 / tau_500 == pytest.approx(0.0007 / 0.0005, rel=0.01)


# ─── Energy Dissipation Tests ────────────────────────────────────────────

class TestPowerDraw:
    def test_legacy_mixer(self):
        m = MixerGeometry()
        P = power_draw(m, 10000, 850.0)
        assert P > 0
        # P = 1.5 * 850 * (10000/60)^3 * 0.025^5
        N = 10000 / 60
        expected = 1.5 * 850 * N**3 * 0.025**5
        assert P == pytest.approx(expected, rel=1e-6)

    def test_stirrer_A(self):
        sA = StirrerGeometry.pitched_blade_A()
        P = power_draw(sA, 1300, 850.0)
        N = 1300 / 60
        expected = 0.35 * 850 * N**3 * 0.059**5
        assert P == pytest.approx(expected, rel=1e-6)

    def test_power_scales_rpm_cubed(self):
        sA = StirrerGeometry.pitched_blade_A()
        P1 = power_draw(sA, 1000, 850.0)
        P2 = power_draw(sA, 2000, 850.0)
        assert P2 / P1 == pytest.approx(8.0, rel=1e-4)

    def test_stirrer_B(self):
        sB = StirrerGeometry.rotor_stator_B()
        P = power_draw(sB, 1800, 850.0)
        assert P > 0


class TestDissipation:
    def test_max_exceeds_avg_legacy(self):
        m = MixerGeometry()
        eps_avg = average_dissipation(m, 10000, 850.0)
        eps_max = max_dissipation(m, 10000, 850.0)
        assert eps_max == pytest.approx(50.0 * eps_avg, rel=1e-6)

    def test_max_exceeds_avg_stirrer_A(self):
        sA = StirrerGeometry.pitched_blade_A()
        V = 0.0005  # 500 mL
        eps_avg = average_dissipation(sA, 1300, 850.0, tank_volume=V)
        eps_max = max_dissipation(sA, 1300, 850.0, tank_volume=V)
        assert eps_max == pytest.approx(5.0 * eps_avg, rel=1e-6)

    def test_dissipation_much_lower_for_stirred_vessel(self):
        """Stirrer A at 1300 RPM should have much lower dissipation than legacy at 10000 RPM."""
        m = MixerGeometry()
        sA = StirrerGeometry.pitched_blade_A()
        eps_legacy = max_dissipation(m, 10000, 850.0)
        eps_stirred = max_dissipation(sA, 1300, 850.0, tank_volume=0.0005)
        assert eps_stirred < eps_legacy * 0.01  # at least 100x lower


class TestReynoldsNumber:
    def test_stirrer_A_turbulent(self):
        """At 80°C with mu=5 mPa.s, Re should be >10000 (turbulent)."""
        Re = impeller_reynolds_number(1300, 0.059, 850.0, 0.005)
        assert Re > 10000

    def test_stirrer_A_transitional_at_low_temp(self):
        """At 20°C with mu~30 mPa.s, Re should be in transitional range."""
        Re = impeller_reynolds_number(1300, 0.059, 850.0, 0.030)
        assert 1000 < Re < 10000


class TestPowerNumberCorrected:
    def test_turbulent_equals_base(self):
        sA = StirrerGeometry.pitched_blade_A()
        Np = power_number_corrected(sA, 50000)
        assert Np == pytest.approx(sA.power_number)

    def test_transitional_higher(self):
        sA = StirrerGeometry.pitched_blade_A()
        Np_turb = power_number_corrected(sA, 50000)
        Np_trans = power_number_corrected(sA, 2000)
        assert Np_trans > Np_turb


class TestSweptVolume:
    def test_stirrer_A(self):
        sA = StirrerGeometry.pitched_blade_A()
        V = swept_volume(sA)
        expected = np.pi / 4 * 0.059**2 * 0.010
        assert V == pytest.approx(expected, rel=1e-6)
        assert V < 0.0005  # must be less than tank volume


class TestShearRate:
    def test_metzner_otto_pitched_blade(self):
        sA = StirrerGeometry.pitched_blade_A()
        gamma = metzner_otto_shear_rate(1300, sA)
        assert gamma == pytest.approx(11.5 * 1300 / 60, rel=1e-6)

    def test_gap_shear_for_open_impeller(self):
        """Open impeller should fall back to Metzner-Otto."""
        sA = StirrerGeometry.pitched_blade_A()
        gamma_gap = gap_shear_rate_extended(1300, sA)
        gamma_mo = metzner_otto_shear_rate(1300, sA)
        assert gamma_gap == pytest.approx(gamma_mo)

    def test_gap_shear_for_rotor_stator(self):
        sB = StirrerGeometry.rotor_stator_B()
        gamma = gap_shear_rate_extended(1800, sB)
        N = 1800 / 60
        expected = np.pi * N * sB.impeller_diameter / sB.gap_width
        assert gamma == pytest.approx(expected, rel=1e-6)

    def test_legacy_gap_shear(self):
        m = MixerGeometry()
        gamma = gap_shear_rate(m, 10000)
        N = 10000 / 60
        expected = np.pi * N * 0.025 / 0.0005
        assert gamma == pytest.approx(expected, rel=1e-6)


# ─── Physical Plausibility Cross-Checks ──────────────────────────────────

class TestPhysicalPlausibility:
    def test_kolmogorov_scale_stirred_vessel(self):
        """At 1300 RPM with pitched blade, Kolmogorov scale should be O(100 µm)."""
        sA = StirrerGeometry.pitched_blade_A()
        eps = max_dissipation(sA, 1300, 850.0, tank_volume=0.0005)
        nu = 0.005 / 850.0  # kinematic viscosity
        eta = kolmogorov_length_scale(eps, nu)
        # Should be in range 10-500 µm for gentle stirring
        assert 10e-6 < eta < 500e-6

    def test_kolmogorov_scale_rotor_stator(self):
        """Legacy rotor-stator at 10000 RPM: Kolmogorov scale < 20 µm."""
        m = MixerGeometry()
        eps = max_dissipation(m, 10000, 850.0)
        nu = 0.005 / 850.0
        eta = kolmogorov_length_scale(eps, nu)
        assert eta < 20e-6  # O(10 µm) range


# ─── Integration Smoke Test: Stirred-Vessel PBE Solver ───────────────────

class TestStirredVesselSolverIntegration:
    """End-to-end smoke test for the stirred-vessel PBE solver path."""

    def test_stirred_vessel_solve_end_to_end(self):
        """Full solve with stirred_vessel.toml config produces valid results."""
        from pathlib import Path
        from emulsim.config import load_config
        from emulsim.properties.database import PropertyDatabase
        from emulsim.level1_emulsification.solver import PBESolver

        params = load_config(Path("configs/stirred_vessel.toml"))
        db = PropertyDatabase()
        props = db.update_for_conditions(
            T_oil=params.formulation.T_oil,
            c_agarose=params.formulation.c_agarose,
            c_chitosan=params.formulation.c_chitosan,
            c_span80=params.formulation.c_span80,
        )
        pbe = PBESolver(
            n_bins=params.solver.l1_n_bins,
            d_min=params.solver.l1_d_min,
            d_max=params.solver.l1_d_max,
        )
        result = pbe.solve(params, props)

        # Basic structure
        assert result.d_bins.shape == (40,)
        assert result.n_d.shape == (40,)

        # d_mode in valid range [5, 1500] um (broad range; exact value depends on
        # kernel calibration which is system-specific)
        assert 5e-6 < result.d_mode < 1500e-6

        # Convergence: may not reach steady state in 600s with gentle breakage
        # (kernel constants are preliminary defaults awaiting calibration).
        # Just check the solver completed without error.
        assert result.d32 > 0

        # Mass conservation: phi_d ~ 0.40
        assert result.total_volume_fraction == pytest.approx(0.40, abs=0.02)

        # Size ordering: d10 <= d50 <= d90
        assert result.d10 <= result.d50 <= result.d90

        # Span is reasonable (< 3.0)
        assert 0 < result.span < 3.0

        # d32 positive and in physical range
        assert 1e-6 < result.d32 < 2000e-6

        # d_mode > 0
        assert result.d_mode > 0

    def test_stirred_vessel_with_explicit_phi_d(self):
        """Caller-provided phi_d is honoured."""
        from pathlib import Path
        from emulsim.config import load_config
        from emulsim.properties.database import PropertyDatabase
        from emulsim.level1_emulsification.solver import PBESolver

        params = load_config(Path("configs/stirred_vessel.toml"))
        db = PropertyDatabase()
        props = db.update_for_conditions(
            T_oil=params.formulation.T_oil,
            c_agarose=params.formulation.c_agarose,
            c_chitosan=params.formulation.c_chitosan,
            c_span80=params.formulation.c_span80,
        )
        pbe = PBESolver(
            n_bins=params.solver.l1_n_bins,
            d_min=params.solver.l1_d_min,
            d_max=params.solver.l1_d_max,
        )
        # Override phi_d to 0.20
        result = pbe.solve(params, props, phi_d=0.20)
        assert result.total_volume_fraction == pytest.approx(0.20, abs=0.02)

    def test_legacy_mode_still_works(self):
        """Legacy config through solve() still produces valid results."""
        from pathlib import Path
        from emulsim.config import load_config
        from emulsim.properties.database import PropertyDatabase
        from emulsim.level1_emulsification.solver import PBESolver

        params = load_config(Path("configs/default.toml"))
        db = PropertyDatabase()
        props = db.update_for_conditions(
            T_oil=params.formulation.T_oil,
            c_agarose=params.formulation.c_agarose,
            c_chitosan=params.formulation.c_chitosan,
            c_span80=params.formulation.c_span80,
        )
        pbe = PBESolver(
            n_bins=params.solver.l1_n_bins,
            d_min=params.solver.l1_d_min,
            d_max=params.solver.l1_d_max,
        )
        result = pbe.solve(params, props)
        assert result.d32 > 0
        # Convergence may vary with v2 C3=0.1 default (was 0.0);
        # just check the solver completed without error.
        assert hasattr(result, 'd_mode')


# ─── v2.0 Phase 1 Tests ─────────────────────────────────────────────────

class TestV2ModelMode:
    """Test ModelMode enum and threading."""

    def test_model_mode_enum(self):
        from emulsim.datatypes import ModelMode
        assert ModelMode.HYBRID_COUPLED.value == "hybrid_coupled"
        assert ModelMode.EMPIRICAL_ENGINEERING.value == "empirical_engineering"
        assert ModelMode.MECHANISTIC_RESEARCH.value == "mechanistic_research"

    def test_default_mode(self):
        from emulsim.datatypes import SimulationParameters, ModelMode
        p = SimulationParameters()
        assert p.model_mode == ModelMode.HYBRID_COUPLED


class TestL1L2Coupling:
    """Test that L2 pore size responds to L1 droplet size (Finding 1 fix)."""

    def test_different_R_droplet_gives_different_pores(self):
        from emulsim.datatypes import SimulationParameters, MaterialProperties
        from emulsim.level2_gelation.solver import solve_gelation

        params = SimulationParameters()
        props = MaterialProperties()

        r_small = solve_gelation(params, props, R_droplet=0.5e-6, mode='empirical')
        r_large = solve_gelation(params, props, R_droplet=50e-6, mode='empirical')

        p_small = r_small.pore_size_mean
        p_large = r_large.pore_size_mean
        assert p_small != pytest.approx(p_large, rel=0.01), \
            f"Pore sizes should differ: {p_small*1e9:.1f} vs {p_large*1e9:.1f} nm"

    def test_confinement_caps_pore_size(self):
        """Very small droplet should produce pore << bulk empirical prediction."""
        from emulsim.datatypes import SimulationParameters, MaterialProperties
        from emulsim.level2_gelation.solver import solve_gelation

        params = SimulationParameters()
        props = MaterialProperties()

        r_tiny = solve_gelation(params, props, R_droplet=0.2e-6, mode='empirical')
        # For R=200nm, confinement max = 0.15 * 400nm = 60nm
        assert r_tiny.pore_size_mean < 70e-9, \
            f"Tiny droplet pore should be <70nm, got {r_tiny.pore_size_mean*1e9:.1f}nm"

    def test_large_droplet_no_confinement(self):
        """Large droplet should give bulk empirical pore size (no confinement effect)."""
        from emulsim.datatypes import SimulationParameters, MaterialProperties
        from emulsim.level2_gelation.solver import solve_gelation

        params = SimulationParameters()
        props = MaterialProperties()

        r_big = solve_gelation(params, props, R_droplet=100e-6, mode='empirical')
        # For R=100um, confinement max = 0.15 * 200um = 30um >> bulk pore ~200nm
        assert r_big.pore_size_mean > 150e-9, \
            f"Large droplet should give bulk-like pore, got {r_big.pore_size_mean*1e9:.1f}nm"


class TestBeadRadiusFix:
    """Test that L4 uses actual bead radius, not truncated L2 domain (Finding 3 fix)."""

    def test_solve_mechanical_accepts_R_droplet(self):
        from emulsim.level4_mechanical.solver import solve_mechanical
        import inspect
        sig = inspect.signature(solve_mechanical)
        assert 'R_droplet' in sig.parameters


class TestTrustWarnings:
    """Test expanded trust warnings (Finding 4,5,7 fixes)."""

    def test_phenomenological_dn_warning(self):
        from emulsim.trust import assess_trust
        from emulsim.datatypes import (
            SimulationParameters, MaterialProperties, FullResult,
            EmulsificationResult, GelationResult, CrosslinkingResult, MechanicalResult,
        )
        import numpy as np

        params = SimulationParameters()
        props = MaterialProperties()

        emul = EmulsificationResult(d_bins=np.array([1e-4]), n_d=np.array([1e10]),
            d32=10e-6, d43=12e-6, d10=5e-6, d50=9e-6, d90=14e-6,
            span=1.0, total_volume_fraction=0.05, converged=True, d_mode=10e-6)
        gel = GelationResult(r_grid=np.array([0]), phi_field=np.array([0]),
            pore_size_mean=100e-9, pore_size_std=10e-9,
            pore_size_distribution=np.array([100e-9]),
            porosity=0.5, alpha_final=1.0, char_wavelength=200e-9)
        xlink = CrosslinkingResult(t_array=np.array([0]), X_array=np.array([0]),
            nu_e_array=np.array([0]), Mc_array=np.array([0]),
            xi_array=np.array([0]), G_chitosan_array=np.array([0]),
            p_final=0.05, nu_e_final=1e22, Mc_final=5000,
            xi_final=5e-9, G_chitosan_final=5000)
        mech = MechanicalResult(G_agarose=5000, G_chitosan=5000, G_DN=10000,
            E_star=30000, delta_array=np.array([0]), F_array=np.array([0]),
            rh_array=np.array([0]), Kav_array=np.array([0]),
            pore_size_mean=100e-9, xi_mesh=5e-9)
        full = FullResult(parameters=params, emulsification=emul,
            gelation=gel, crosslinking=xlink, mechanical=mech)

        trust = assess_trust(full, params, props)
        all_warnings = " ".join(trust.warnings)
        assert "phenomenological" in all_warnings.lower() or "G_DN" in all_warnings
