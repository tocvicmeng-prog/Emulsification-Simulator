"""Unit tests for Phase 1 data layer: equipment types, kernel config, validation, config loading.

Tests cover:
  - Enums (VesselType, StirrerType, HeatingStrategy, BreakageModel)
  - VesselGeometry factory methods and properties
  - StirrerGeometry factory methods and measured dimensions
  - HeatingConfig factory methods
  - KernelConfig per-stirrer dispatch
  - EmulsificationParameters dual-mode (legacy vs stirred-vessel)
  - FormulationParameters volumetric fields
  - SimulationParameters.validate() for both modes
  - Config loading from TOML (legacy and stirred-vessel)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from emulsim.datatypes import (
    BreakageModel,
    EmulsificationParameters,
    FormulationParameters,
    HeatingConfig,
    HeatingStrategy,
    KernelConfig,
    MixerGeometry,
    SimulationParameters,
    StirrerGeometry,
    StirrerType,
    VesselGeometry,
    VesselType,
)
from emulsim.config import load_config


# ─── Enum Tests ──────────────────────────────────────────────────────────

class TestEnums:
    def test_vessel_type_values(self):
        assert VesselType.GLASS_BEAKER.value == "glass_beaker"
        assert VesselType.JACKETED_VESSEL.value == "jacketed_vessel"

    def test_stirrer_type_values(self):
        assert StirrerType.PITCHED_BLADE.value == "pitched_blade"
        assert StirrerType.ROTOR_STATOR_SMALL.value == "rotor_stator_small"
        assert StirrerType.ROTOR_STATOR_LEGACY.value == "rotor_stator_legacy"

    def test_heating_strategy_values(self):
        assert HeatingStrategy.FLAT_PLATE.value == "flat_plate"
        assert HeatingStrategy.HOT_WATER_JACKET.value == "hot_water_jacket"
        assert HeatingStrategy.ISOTHERMAL.value == "isothermal"

    def test_breakage_model_values(self):
        assert BreakageModel.ALOPAEUS.value == "alopaeus"
        assert BreakageModel.COULALOGLOU_TAVLARIDES.value == "coulaloglou_tavlarides"


# ─── VesselGeometry Tests ────────────────────────────────────────────────

class TestVesselGeometry:
    def test_glass_beaker_defaults(self):
        v = VesselGeometry.glass_beaker()
        assert v.vessel_type == VesselType.GLASS_BEAKER
        assert v.inner_diameter == pytest.approx(0.100, abs=1e-6)
        assert v.wall_thickness == pytest.approx(0.0015, abs=1e-6)
        assert v.height == pytest.approx(0.130, abs=1e-6)
        assert v.working_volume == pytest.approx(0.0005, abs=1e-9)

    def test_jacketed_vessel_defaults(self):
        v = VesselGeometry.jacketed_vessel()
        assert v.vessel_type == VesselType.JACKETED_VESSEL
        assert v.inner_diameter == pytest.approx(0.092, abs=1e-6)
        assert v.wall_thickness == pytest.approx(0.002, abs=1e-6)
        assert v.height == pytest.approx(0.160, abs=1e-6)

    def test_custom_working_volume(self):
        v = VesselGeometry.glass_beaker(working_volume=0.0003)
        assert v.working_volume == pytest.approx(0.0003, abs=1e-9)

    def test_cross_section_area(self):
        v = VesselGeometry.glass_beaker()
        expected = np.pi / 4 * 0.100**2  # ~7.854e-3 m²
        assert v.cross_section_area == pytest.approx(expected, rel=1e-6)

    def test_liquid_height(self):
        v = VesselGeometry.glass_beaker()
        # 500 mL / (π/4 × 0.1²) ≈ 0.0637 m
        assert v.liquid_height == pytest.approx(0.0637, rel=0.01)
        assert v.liquid_height < v.height  # must fit in vessel

    def test_volume_range(self):
        v = VesselGeometry.glass_beaker()
        assert v.working_volume_min == pytest.approx(0.00025)
        assert v.working_volume_max == pytest.approx(0.0007)


# ─── StirrerGeometry Tests ───────────────────────────────────────────────

class TestStirrerGeometry:
    def test_pitched_blade_A_dimensions(self):
        """Verify measured dimensions from photographs."""
        s = StirrerGeometry.pitched_blade_A()
        assert s.stirrer_type == StirrerType.PITCHED_BLADE
        assert s.impeller_diameter == pytest.approx(0.059, abs=1e-6)  # 59 mm
        assert s.shaft_diameter == pytest.approx(0.008, abs=1e-6)     # 8 mm
        assert s.blade_thickness == pytest.approx(0.001, abs=1e-6)    # 1 mm
        assert s.blade_height == pytest.approx(0.010, abs=1e-6)       # 10 mm
        assert s.blade_length == pytest.approx(0.009, abs=1e-6)       # 9 mm
        assert s.blade_angle == pytest.approx(10.0)
        assert s.max_rpm == pytest.approx(2000.0)
        assert s.has_stator is False

    def test_rotor_stator_B_dimensions(self):
        """Verify measured dimensions from photographs."""
        s = StirrerGeometry.rotor_stator_B()
        assert s.stirrer_type == StirrerType.ROTOR_STATOR_SMALL
        assert s.impeller_diameter == pytest.approx(0.0257, abs=1e-6)  # 25.7 mm
        assert s.stator_diameter == pytest.approx(0.03203, abs=1e-6)   # 32.03 mm
        assert s.blade_thickness == pytest.approx(0.002, abs=1e-6)     # 2 mm
        assert s.wall_thickness == pytest.approx(0.0022, abs=1e-6)     # 2.2 mm
        assert s.wall_height == pytest.approx(0.018, abs=1e-6)         # 18 mm
        assert s.perforation_diameter == pytest.approx(0.003, abs=1e-6) # 3 mm
        assert s.max_rpm == pytest.approx(9000.0)
        assert s.has_stator is True

    def test_rotor_stator_legacy(self):
        s = StirrerGeometry.rotor_stator_legacy()
        assert s.stirrer_type == StirrerType.ROTOR_STATOR_LEGACY
        assert s.impeller_diameter == pytest.approx(0.025)
        assert s.gap_width == pytest.approx(0.0005)
        assert s.dissipation_ratio == pytest.approx(50.0)
        assert s.max_rpm == pytest.approx(25000.0)

    def test_tip_speed(self):
        s = StirrerGeometry.pitched_blade_A()
        # v_tip = π × 0.059 × 2000/60 ≈ 6.18 m/s
        assert s.tip_speed == pytest.approx(6.18, rel=0.01)

    def test_gap_width_stirrer_B(self):
        s = StirrerGeometry.rotor_stator_B()
        # gap = (32.03 - 25.7) / 2 ≈ 3.165 mm
        assert s.gap_width == pytest.approx(0.003165, rel=0.01)


# ─── HeatingConfig Tests ─────────────────────────────────────────────────

class TestHeatingConfig:
    def test_flat_plate(self):
        h = HeatingConfig.flat_plate()
        assert h.strategy == HeatingStrategy.FLAT_PLATE
        assert h.T_initial == pytest.approx(353.15)    # 80°C
        assert h.T_final == pytest.approx(293.15)      # 20°C
        assert h.cooldown_time == pytest.approx(5400.0) # 1.5 h

    def test_hot_water_jacket(self):
        h = HeatingConfig.hot_water_jacket()
        assert h.strategy == HeatingStrategy.HOT_WATER_JACKET
        assert h.jacket_water_temperature == pytest.approx(358.15)  # 85°C

    def test_isothermal(self):
        h = HeatingConfig.isothermal(363.15)
        assert h.strategy == HeatingStrategy.ISOTHERMAL
        assert h.T_initial == pytest.approx(363.15)
        assert h.T_final == pytest.approx(363.15)
        assert h.cooldown_time == pytest.approx(0.0)


# ─── KernelConfig Tests ──────────────────────────────────────────────────

class TestKernelConfig:
    def test_for_pitched_blade(self):
        k = KernelConfig.for_pitched_blade()
        assert k.breakage_model == BreakageModel.ALOPAEUS  # switched from CT (viscous sub-range)
        assert k.breakage_C1 == pytest.approx(0.04)
        assert k.breakage_C3 == pytest.approx(2.0)  # strong viscous correction
        assert k.phi_d_correction is True
        assert k.coalescence_exponent == 2

    def test_for_rotor_stator_legacy(self):
        # F1 fix (2026-04-17) enabled phi_d_correction + exponent=2 on the
        # legacy preset to match for_rotor_stator_small; see KernelConfig
        # docstring for the high-RPM nonphysical-d32 rationale.
        k = KernelConfig.for_rotor_stator_legacy()
        assert k.breakage_model == BreakageModel.ALOPAEUS
        assert k.breakage_C1 == pytest.approx(0.986)
        assert k.phi_d_correction is True
        assert k.coalescence_exponent == 2

    def test_for_rotor_stator_small(self):
        k = KernelConfig.for_rotor_stator_small()
        assert k.breakage_model == BreakageModel.ALOPAEUS
        assert k.phi_d_correction is True

    def test_dispatch(self):
        for st in StirrerType:
            k = KernelConfig.for_stirrer_type(st)
            assert isinstance(k, KernelConfig)


# ─── EmulsificationParameters Tests ──────────────────────────────────────

class TestEmulsificationParameters:
    def test_legacy_defaults(self):
        e = EmulsificationParameters()
        assert e.mode == "rotor_stator_legacy"
        assert e.rpm == pytest.approx(10000.0)
        assert e.vessel is None
        assert e.stirrer is None

    def test_stirred_vessel_auto_populate(self):
        e = EmulsificationParameters(mode="stirred_vessel", rpm=1300)
        assert e.vessel is not None
        assert e.stirrer is not None
        assert e.heating is not None
        assert e.kernels is not None
        assert e.stirrer.stirrer_type == StirrerType.PITCHED_BLADE

    def test_effective_properties_legacy(self):
        e = EmulsificationParameters()
        assert e.effective_tank_volume == pytest.approx(0.0005)
        assert e.effective_impeller_diameter == pytest.approx(0.025)
        assert e.effective_power_number == pytest.approx(1.5)

    def test_effective_properties_stirred(self):
        e = EmulsificationParameters(mode="stirred_vessel", rpm=1300)
        assert e.effective_tank_volume == pytest.approx(0.0005)
        assert e.effective_impeller_diameter == pytest.approx(0.059)
        assert e.effective_power_number == pytest.approx(0.35)


# ─── FormulationParameters Tests ─────────────────────────────────────────

class TestFormulationParameters:
    def test_phi_d_from_volumes(self):
        f = FormulationParameters(v_oil_span80_mL=300, v_polysaccharide_mL=200)
        assert f.phi_d_from_volumes == pytest.approx(0.40, abs=0.01)

    def test_total_working_volume(self):
        f = FormulationParameters(v_oil_span80_mL=300, v_polysaccharide_mL=200)
        assert f.total_working_volume_mL == pytest.approx(500.0)


# ─── Validation Tests ────────────────────────────────────────────────────

class TestValidation:
    def test_legacy_validates_clean(self):
        p = SimulationParameters()
        assert p.validate() == []

    def test_stirred_vessel_validates_clean(self):
        p = SimulationParameters(
            emulsification=EmulsificationParameters(
                mode="stirred_vessel", rpm=1300, t_emulsification=600
            ),
            formulation=FormulationParameters(phi_d=0.40),
        )
        assert p.validate() == []

    def test_rpm_exceeds_max(self):
        p = SimulationParameters(
            emulsification=EmulsificationParameters(
                mode="stirred_vessel", rpm=3000
            ),
        )
        errors = p.validate()
        assert any("exceeds stirrer max" in e for e in errors)

    def test_flat_plate_with_jacket_vessel(self):
        p = SimulationParameters(
            emulsification=EmulsificationParameters(
                mode="stirred_vessel", rpm=1300,
                vessel=VesselGeometry.jacketed_vessel(),
                heating=HeatingConfig.flat_plate(),
            ),
        )
        errors = p.validate()
        assert any("Flat-plate" in e for e in errors)

    def test_jacket_with_beaker(self):
        p = SimulationParameters(
            emulsification=EmulsificationParameters(
                mode="stirred_vessel", rpm=1300,
                vessel=VesselGeometry.glass_beaker(),
                heating=HeatingConfig.hot_water_jacket(),
            ),
        )
        errors = p.validate()
        assert any("jacketed vessel" in e for e in errors)

    def test_volume_mismatch(self):
        p = SimulationParameters(
            emulsification=EmulsificationParameters(
                mode="stirred_vessel", rpm=1300,
                vessel=VesselGeometry.glass_beaker(working_volume=0.0006),
            ),
            formulation=FormulationParameters(
                v_oil_span80_mL=300, v_polysaccharide_mL=200,
            ),
        )
        errors = p.validate()
        assert any("volumes" in e.lower() for e in errors)

    def test_unknown_mode(self):
        p = SimulationParameters(
            emulsification=EmulsificationParameters(mode="unknown"),
        )
        errors = p.validate()
        assert any("Unknown" in e for e in errors)


# ─── Config Loading Tests ────────────────────────────────────────────────

class TestConfigLoading:
    def test_load_legacy_config(self):
        p = load_config(Path("configs/default.toml"))
        assert p.emulsification.mode == "rotor_stator_legacy"
        assert p.emulsification.rpm == pytest.approx(10000)
        assert p.emulsification.mixer.rotor_diameter == pytest.approx(0.025)
        assert p.validate() == []

    def test_load_stirred_vessel_config(self):
        p = load_config(Path("configs/stirred_vessel.toml"))
        assert p.emulsification.mode == "stirred_vessel"
        assert p.emulsification.rpm == pytest.approx(1300)
        assert p.emulsification.vessel.vessel_type == VesselType.GLASS_BEAKER
        assert p.emulsification.stirrer.stirrer_type == StirrerType.PITCHED_BLADE
        assert p.emulsification.heating.strategy == HeatingStrategy.FLAT_PLATE
        assert p.emulsification.kernels.breakage_model == BreakageModel.ALOPAEUS
        assert p.formulation.phi_d == pytest.approx(0.40)
        assert p.solver.l1_n_bins == 40
        assert p.solver.l1_d_max == pytest.approx(2000e-6)
        assert p.validate() == []


# ─── Backward Compatibility Tests ────────────────────────────────────────

class TestBackwardCompatibility:
    def test_mixer_geometry_still_works(self):
        """MixerGeometry should still be importable and functional."""
        m = MixerGeometry()
        assert m.rotor_diameter == pytest.approx(0.025)
        assert m.power_number == pytest.approx(1.5)

    def test_legacy_emul_params_default(self):
        """Default EmulsificationParameters should behave as before."""
        e = EmulsificationParameters()
        assert e.rpm == pytest.approx(10000.0)
        assert e.mixer.tank_volume == pytest.approx(0.0005)
