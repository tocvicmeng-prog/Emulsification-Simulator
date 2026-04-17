"""Core data structures for the emulsification simulation pipeline.

Supports two hardware modes:
  - Legacy rotor-stator homogeniser (original, 2 µm target)
  - Stirred-vessel double-emulsification (new, 100 µm target)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional

import numpy as np


# ─── Platform / polymer family ───────────────────────────────────────────


class PolymerFamily(Enum):
    """Microsphere polymer family (Node F1-a, v8.0).

    Drives L2 gelation-solver dispatch. AGAROSE_CHITOSAN is the legacy
    v7.x platform (thermal TIPS + optional covalent crosslinking);
    ALGINATE uses diffusion-limited Ca²⁺ ionic gelation; CELLULOSE and
    PLGA are stubs for future F1-b / F1-c work.
    """
    AGAROSE_CHITOSAN = "agarose_chitosan"
    ALGINATE = "alginate"
    CELLULOSE = "cellulose"
    PLGA = "plga"


# ─── Equipment Enums ─────────────────────────────────────────────────────

class VesselType(Enum):
    """Reaction vessel type for the double-emulsification process."""
    GLASS_BEAKER = "glass_beaker"
    JACKETED_VESSEL = "jacketed_vessel"


class StirrerType(Enum):
    """Stirrer/impeller type."""
    PITCHED_BLADE = "pitched_blade"             # Stirrer A — open impeller, ≤2000 RPM
    ROTOR_STATOR_SMALL = "rotor_stator_small"   # Stirrer B — small homogeniser, ≤9000 RPM
    ROTOR_STATOR_LEGACY = "rotor_stator_legacy" # Original 25 mm rotor-stator


class HeatingStrategy(Enum):
    """Heating/cooling strategy during emulsification."""
    FLAT_PLATE = "flat_plate"           # External hot plate (beaker only)
    HOT_WATER_JACKET = "hot_water_jacket"  # Circulating jacket (jacketed vessel only)
    ISOTHERMAL = "isothermal"           # Constant T_oil (legacy mode)


class ModelEvidenceTier(Enum):
    """Evidence quality tier for a model output.

    Every numeric result should carry an evidence tier indicating how
    much the prediction can be trusted for decision-making.
    Ordered from strongest to weakest.
    """
    VALIDATED_QUANTITATIVE = "validated_quantitative"
    """Calibrated against experimental data for this specific system."""

    CALIBRATED_LOCAL = "calibrated_local"
    """Calibrated against data from an analogous system."""

    SEMI_QUANTITATIVE = "semi_quantitative"
    """Empirical model, not locally calibrated. Trends reliable, magnitudes approximate."""

    QUALITATIVE_TREND = "qualitative_trend"
    """Directional predictions only. Numeric values are order-of-magnitude estimates."""

    UNSUPPORTED = "unsupported"
    """Model not applicable to this chemistry or regime. Output should not be used."""


class ModelMode(Enum):
    """Scientific operating mode for the simulation pipeline.

    Controls which model types are permitted and what claims are defensible.
    """
    EMPIRICAL_ENGINEERING = "empirical_engineering"
    """Fast screening. Calibration-driven. Only trend/ranking claims allowed."""

    HYBRID_COUPLED = "hybrid_coupled"
    """Main production mode. Empirical where evidence exists, mechanistic
    where coupling is defensible. Default."""

    MECHANISTIC_RESEARCH = "mechanistic_research"
    """Exploratory hypothesis-testing. Slower. Not for production decisions."""


# ─── Equipment Geometry Dataclasses ──────────────────────────────────────

@dataclass
class VesselGeometry:
    """Reaction vessel geometry.

    Working liquid = paraffin oil + polysaccharide solution + surfactant.
    """
    vessel_type: VesselType = VesselType.GLASS_BEAKER
    inner_diameter: float = 0.100       # [m]
    wall_thickness: float = 0.0015      # [m]
    height: float = 0.130               # [m]
    material: str = "borosilicate_glass"
    working_volume: float = 0.0005      # [m³] (500 mL default)
    working_volume_min: float = 0.00025 # [m³] (250 mL)
    working_volume_max: float = 0.0007  # [m³] (700 mL)

    @property
    def cross_section_area(self) -> float:
        """Internal cross-section area [m²]."""
        return np.pi / 4 * self.inner_diameter ** 2

    @property
    def liquid_height(self) -> float:
        """Height of working liquid in the vessel [m]."""
        return self.working_volume / self.cross_section_area

    @classmethod
    def glass_beaker(cls, working_volume: float = 0.0005) -> VesselGeometry:
        """Factory: standard glass beaker (Ø100 mm × 130 mm)."""
        return cls(
            vessel_type=VesselType.GLASS_BEAKER,
            inner_diameter=0.100,       # 10 cm
            wall_thickness=0.0015,      # 1.5 mm
            height=0.130,               # 13 cm
            material="borosilicate_glass",
            working_volume=working_volume,
            working_volume_min=0.00025, # 250 mL
            working_volume_max=0.0007,  # 700 mL
        )

    @classmethod
    def jacketed_vessel(cls, working_volume: float = 0.0005) -> VesselGeometry:
        """Factory: jacketed glass vessel (Ø92 mm × 160 mm)."""
        return cls(
            vessel_type=VesselType.JACKETED_VESSEL,
            inner_diameter=0.092,       # 9.2 cm
            wall_thickness=0.002,       # 2 mm (inner + outer wall)
            height=0.160,               # 16 cm
            material="borosilicate_glass",
            working_volume=working_volume,
            working_volume_min=0.0002,  # 200 mL
            working_volume_max=0.0007,  # 700 mL
        )


@dataclass
class StirrerGeometry:
    """Stirrer/impeller geometry for the emulsification process.

    Dimensions from measurement photographs (2026-03-27).
    Supports pitched-blade (Stirrer A), small rotor-stator (Stirrer B),
    and legacy rotor-stator (original 25 mm) configurations.
    """
    stirrer_type: StirrerType = StirrerType.PITCHED_BLADE
    impeller_diameter: float = 0.059    # [m] outer diameter of active element
    shaft_diameter: float = 0.008       # [m]
    blade_count: int = 6                # [-] number of blades/fins
    blade_angle: float = 10.0           # [°] angle relative to tangent
    blade_thickness: float = 0.001      # [m]
    blade_height: float = 0.010         # [m] parallel to axis
    blade_length: float = 0.009         # [m] perpendicular to axis (fin length)
    power_number: float = 0.35          # [-] estimated for 10° pitched, alternately bent
    max_rpm: float = 2000.0             # [rev/min]
    # Rotor-stator specific fields (0 for open impellers)
    has_stator: bool = False
    stator_diameter: float = 0.0        # [m]
    gap_width: float = 0.0              # [m]
    wall_height: float = 0.0            # [m] outer wall height (Stirrer B)
    wall_thickness: float = 0.0         # [m] outer wall thickness
    perforation_diameter: float = 0.0   # [m] hole diameter in perforated stator
    # Computed from stirrer_type, not user-set
    dissipation_ratio: float = 5.0      # ε_max / ε_avg

    @property
    def tip_speed(self) -> float:
        """Maximum tip speed [m/s] at max_rpm."""
        return np.pi * self.impeller_diameter * self.max_rpm / 60.0

    @classmethod
    def pitched_blade_A(cls) -> StirrerGeometry:
        """Factory: Stirrer A — pitched-blade impeller.

        Measured dimensions:
          - Shaft Ø8 mm, fin outer Ø59 mm (root and top)
          - Fin length 9 mm (⊥ axis), fin height 10 mm (∥ axis)
          - Blade thickness 1 mm, 10° angle, alternately bent
          - Opposite fin tip distance 18 mm
        """
        return cls(
            stirrer_type=StirrerType.PITCHED_BLADE,
            impeller_diameter=0.059,        # Ø59 mm outer diameter
            shaft_diameter=0.008,           # Ø8 mm
            blade_count=6,                  # estimated from photos (radial fins)
            blade_angle=10.0,               # 10° from tangent
            blade_thickness=0.001,          # 1 mm
            blade_height=0.010,             # 10 mm parallel to axis
            blade_length=0.009,             # 9 mm perpendicular to axis
            # ASSUMPTION: Np ≈ 0.35 for small-angle (10°) alternately bent blades.
            # Literature Np for 45° PBT is ~1.3; 10° with alternating bend is
            # closer to a flat disc turbine with minimal axial pumping.
            power_number=0.35,
            max_rpm=2000.0,
            has_stator=False,
            dissipation_ratio=5.0,          # typical for open impellers
        )

    @classmethod
    def rotor_stator_B(cls) -> StirrerGeometry:
        """Factory: Stirrer B — small rotor-stator homogeniser.

        Measured dimensions:
          - Outer wall Ø32.03 mm, wall thickness 2.2 mm, wall height 18 mm
          - Cross blade root Ø8.5 mm, blade edge Ø25.7 mm, blade thickness 2 mm
          - Perforated outer wall with Ø3 mm holes
          - Closed top, open bottom — centrifugal ejection through peripheral holes
        """
        return cls(
            stirrer_type=StirrerType.ROTOR_STATOR_SMALL,
            impeller_diameter=0.0257,       # 25.7 mm cross-blade edge (active rotor)
            shaft_diameter=0.0085,          # 8.5 mm root diameter
            blade_count=4,                  # cross-shaped (4 arms)
            blade_angle=0.0,                # radial (cross blade)
            blade_thickness=0.002,          # 2 mm
            blade_height=0.018,             # 18 mm (outer wall height)
            blade_length=0.0,               # not applicable for cross-blade
            power_number=2.0,               # rotor-stator, estimated from Padron (2005)
            max_rpm=9000.0,
            has_stator=True,
            stator_diameter=0.03203,        # 32.03 mm outer wall
            gap_width=(0.03203 - 0.0257) / 2,  # ~3.2 mm radial gap
            wall_height=0.018,              # 18 mm
            wall_thickness=0.0022,          # 2.2 mm
            perforation_diameter=0.003,     # 3 mm holes
            dissipation_ratio=25.0,         # intermediate: open perforations reduce confinement
        )

    @classmethod
    def rotor_stator_legacy(cls) -> StirrerGeometry:
        """Factory: original 25 mm rotor-stator (backward compatibility)."""
        return cls(
            stirrer_type=StirrerType.ROTOR_STATOR_LEGACY,
            impeller_diameter=0.025,
            shaft_diameter=0.008,
            blade_count=4,
            blade_angle=0.0,
            blade_thickness=0.001,
            blade_height=0.010,
            blade_length=0.0,
            power_number=1.5,
            max_rpm=25000.0,
            has_stator=True,
            stator_diameter=0.026,
            gap_width=0.0005,
            wall_height=0.010,
            wall_thickness=0.001,
            perforation_diameter=0.0,
            dissipation_ratio=50.0,
        )


@dataclass
class HeatingConfig:
    """Heating/cooling configuration during emulsification.

    Empirical calibration data:
      - Flat plate at 150°C → 300 mL oil/Span80 reaches 80°C steady state
      - After off: 500 mL cools 80°C → 20°C in ~1.5 hours
      - Jacket at 85°C circulating → 300 mL oil/Span80 reaches 80°C steady state
    """
    strategy: HeatingStrategy = HeatingStrategy.FLAT_PLATE
    heater_temperature: float = 423.15  # [K] (150°C for flat plate)
    T_initial: float = 353.15           # [K] (80°C — oil temp at mixing start)
    T_final: float = 293.15             # [K] (20°C — ambient)
    cooldown_time: float = 5400.0       # [s] (1.5 hours for 500 mL with flat plate)
    jacket_water_temperature: float = 358.15  # [K] (85°C circulating water)

    @classmethod
    def flat_plate(cls) -> HeatingConfig:
        """Factory: flat-plate heater (for glass beaker only)."""
        return cls(
            strategy=HeatingStrategy.FLAT_PLATE,
            heater_temperature=423.15,      # 150°C
            T_initial=353.15,               # 80°C
            T_final=293.15,                 # 20°C
            cooldown_time=5400.0,           # 1.5 h
        )

    @classmethod
    def hot_water_jacket(cls) -> HeatingConfig:
        """Factory: hot-water jacket (for jacketed vessel only)."""
        return cls(
            strategy=HeatingStrategy.HOT_WATER_JACKET,
            heater_temperature=0.0,         # not used
            T_initial=353.15,               # 80°C
            T_final=293.15,                 # 20°C
            cooldown_time=7200.0,           # ASSUMPTION: ~2 h for jacket cooling (slower)
            jacket_water_temperature=358.15,  # 85°C
        )

    @classmethod
    def isothermal(cls, T: float = 363.15) -> HeatingConfig:
        """Factory: constant temperature (legacy mode)."""
        return cls(
            strategy=HeatingStrategy.ISOTHERMAL,
            heater_temperature=0.0,
            T_initial=T,
            T_final=T,
            cooldown_time=0.0,
        )


# ─── Kernel Configuration ────────────────────────────────────────────────

class BreakageModel(Enum):
    """Breakage kernel model selection."""
    ALOPAEUS = "alopaeus"               # Alopaeus et al. — suited for rotor-stator
    COULALOGLOU_TAVLARIDES = "coulaloglou_tavlarides"  # CT (1977) original — suited for stirred vessels


class CoalescenceModel(Enum):
    """Coalescence kernel model selection."""
    COULALOGLOU_TAVLARIDES = "coulaloglou_tavlarides"  # Standard CT coalescence


@dataclass
class KernelConfig:
    """Breakage and coalescence kernel configuration.

    Constants are system-dependent.  Factory methods provide calibrated
    defaults for each stirrer type.

    Uncertainty bands (Liao & Lucas, 2009 review):
      - Alopaeus C1: [0.5, 2.0] (default 0.986)
      - Alopaeus C2: [0.005, 0.03] (default 0.0115)
      - CT breakage C1: [0.003, 0.01] (default 0.00481)
      - CT breakage C2: [0.05, 0.15] (default 0.08)
      - CT coalescence C4: [1e-4, 5e-4] (default 2.17e-4)
      - CT coalescence C5: [1e13, 5e13] (default 2.28e13)
    These ranges span published values across different systems.
    Calibration against experimental DSD data is recommended.
    """
    breakage_model: BreakageModel = BreakageModel.ALOPAEUS
    coalescence_model: CoalescenceModel = CoalescenceModel.COULALOGLOU_TAVLARIDES
    # Breakage constants
    breakage_C1: float = 0.986          # [-] Alopaeus default; range [0.5, 2.0]
    breakage_C2: float = 0.0115         # [-] Alopaeus default; range [0.005, 0.03]
    breakage_C3: float = 0.0            # [-] Alopaeus viscous correction. Default 0.0: for aqueous
                                         # dispersed phases (mu_d ~ 0.01-0.1 Pa·s at emulsification T),
                                         # the viscous correction is physically negligible. The Vi cap
                                         # in breakage_rate_alopaeus() prevents exponential blowout for
                                         # any C3 > 0 value. Set C3 > 0 only for highly viscous
                                         # dispersed phases (mu_d > 1 Pa·s) after calibration.
    # Coalescence constants (Coulaloglou-Tavlarides)
    coalescence_C4: float = 2.17e-4     # [-]; range [1e-4, 5e-4]
    coalescence_C5: float = 2.28e13     # [-]; range [1e13, 5e13]
    # Concentrated-emulsion correction
    phi_d_correction: bool = False       # enable coalescence damping
    coalescence_exponent: int = 1        # coalescence ~ 1/(1+phi_d)^n; n=1 legacy, n=2 recommended for phi_d>0.3

    @classmethod
    def for_rotor_stator_legacy(cls) -> KernelConfig:
        """Factory: original rotor-stator calibration (Alopaeus breakage).

        F1 fix (2026-04-17): phi_d_correction and coalescence_exponent=2 are
        now enabled to match the for_rotor_stator_small preset. With them
        disabled (the prior default), high-RPM runs at viscous-dispersed-phase
        conditions exhibited a nonphysical d32 increase with RPM because
        coalescence damping `exp(-C5*mu_c*rho*eps/sigma^2 * d_h^4)` collapses
        for small d_h and the bare CT collision frequency dominates. The
        crowding correction `1/(1+phi_d)^2` adds the missing damping at small
        scales without changing the well-behaved behaviour at low RPM.
        """
        return cls(
            breakage_model=BreakageModel.ALOPAEUS,
            breakage_C1=0.986,
            breakage_C2=0.0115,
            breakage_C3=0.0,                # disabled for aqueous dispersed phases (Vi cap handles high-mu_d)
            coalescence_C4=2.17e-4,
            coalescence_C5=2.28e13,
            phi_d_correction=True,
            coalescence_exponent=2,
        )

    @classmethod
    def for_pitched_blade(cls) -> KernelConfig:
        """Factory: stirred-vessel with pitched-blade impeller.

        Uses Alopaeus breakage kernel (with viscous sub-range correction)
        instead of CT, because at typical stirred-vessel conditions
        d_mode ~ eta_K (Kolmogorov scale), placing droplets in the
        viscous/transitional breakage regime where the CT inertial
        assumption (d >> eta_K) breaks down.

        Alopaeus C1/C2 are reduced from the rotor-stator defaults to
        reflect the lower turbulent intensity in stirred vessels.

        Coalescence exponent n=2 for concentrated emulsions (phi_d~0.40).
        """
        return cls(
            breakage_model=BreakageModel.ALOPAEUS,
            breakage_C1=0.04,               # reduced ~25x from rotor-stator (0.986)
            breakage_C2=0.0115,             # Alopaeus default surface-tension term
            breakage_C3=2.0,                # strong viscous correction: suppresses breakage of small droplets
            coalescence_C4=2.17e-4,
            coalescence_C5=2.28e13,
            phi_d_correction=True,
            coalescence_exponent=2,         # 1/(1+phi_d)^2 for phi_d~0.40
        )

    @classmethod
    def for_rotor_stator_small(cls) -> KernelConfig:
        """Factory: small rotor-stator (Stirrer B) — intermediate calibration."""
        return cls(
            breakage_model=BreakageModel.ALOPAEUS,
            breakage_C1=0.986,
            breakage_C2=0.0115,
            breakage_C3=0.0,
            coalescence_C4=2.17e-4,
            coalescence_C5=2.28e13,
            phi_d_correction=True,
            coalescence_exponent=2,
        )

    @classmethod
    def for_stirrer_type(cls, stirrer_type: StirrerType) -> KernelConfig:
        """Select kernel config based on stirrer type."""
        _dispatch = {
            StirrerType.PITCHED_BLADE: cls.for_pitched_blade,
            StirrerType.ROTOR_STATOR_SMALL: cls.for_rotor_stator_small,
            StirrerType.ROTOR_STATOR_LEGACY: cls.for_rotor_stator_legacy,
        }
        return _dispatch[stirrer_type]()


# ─── Simulation Parameters ───────────────────────────────────────────────

@dataclass
class MixerGeometry:
    """Rotor-stator mixer geometry (legacy — kept for backward compatibility).

    New code should use StirrerGeometry + VesselGeometry instead.
    """
    rotor_diameter: float = 0.025       # [m] (25 mm default)
    stator_diameter: float = 0.026      # [m]
    gap_width: float = 0.0005           # [m] (0.5 mm)
    tank_volume: float = 0.0005         # [m³] (500 mL)
    power_number: float = 1.5           # [-]
    dissipation_ratio: float = 50.0     # ε_max / ε_avg


@dataclass
class EmulsificationParameters:
    """Process parameters for Level 1.

    Supports two modes:
      - Legacy: uses ``mixer`` (MixerGeometry) — rotor-stator at 10,000 RPM
      - Stirred-vessel: uses ``vessel``, ``stirrer``, ``heating``, ``kernels``
    The ``mode`` field selects which set of fields the solver reads.
    """
    mode: str = "rotor_stator_legacy"   # "rotor_stator_legacy" or "stirred_vessel"
    rpm: float = 10000.0                # [rev/min]
    t_emulsification: float = 60.0      # [s]
    # Legacy fields
    mixer: MixerGeometry = field(default_factory=MixerGeometry)
    # New stirred-vessel fields
    vessel: Optional[VesselGeometry] = None
    stirrer: Optional[StirrerGeometry] = None
    heating: Optional[HeatingConfig] = None
    kernels: Optional[KernelConfig] = None

    def __post_init__(self) -> None:
        """Auto-populate stirred-vessel fields if mode is set."""
        if self.mode == "stirred_vessel":
            if self.vessel is None:
                self.vessel = VesselGeometry.glass_beaker()
            if self.stirrer is None:
                self.stirrer = StirrerGeometry.pitched_blade_A()
            if self.heating is None:
                # Pick vessel-compatible heating strategy
                if (self.vessel is not None
                        and getattr(self.vessel.vessel_type, 'value', '') == "jacketed_vessel"):
                    self.heating = HeatingConfig.hot_water_jacket()
                else:
                    self.heating = HeatingConfig.flat_plate()
            if self.kernels is None:
                self.kernels = KernelConfig.for_stirrer_type(
                    self.stirrer.stirrer_type
                )

    @property
    def effective_tank_volume(self) -> float:
        """Tank volume [m³] from either legacy mixer or new vessel."""
        if self.mode == "stirred_vessel" and self.vessel is not None:
            return self.vessel.working_volume
        return self.mixer.tank_volume

    @property
    def effective_impeller_diameter(self) -> float:
        """Impeller/rotor diameter [m] from either legacy or new stirrer."""
        if self.mode == "stirred_vessel" and self.stirrer is not None:
            return self.stirrer.impeller_diameter
        return self.mixer.rotor_diameter

    @property
    def effective_power_number(self) -> float:
        """Power number [-] from either legacy or new stirrer."""
        if self.mode == "stirred_vessel" and self.stirrer is not None:
            return self.stirrer.power_number
        return self.mixer.power_number

    @property
    def effective_dissipation_ratio(self) -> float:
        """ε_max / ε_avg from either legacy or new stirrer."""
        if self.mode == "stirred_vessel" and self.stirrer is not None:
            return self.stirrer.dissipation_ratio
        return self.mixer.dissipation_ratio


@dataclass
class FormulationParameters:
    """Chemical formulation parameters.

    In stirred-vessel mode, the surfactant concentration and dispersed-phase
    volume fraction are computed from volumetric quantities (v_oil_span80_mL,
    v_polysaccharide_mL, c_span80_vol_pct).
    """
    c_agarose: float = 42.0             # [kg/m³] (4.2% w/v)
    c_chitosan: float = 18.0            # [kg/m³] (1.8% w/v)
    c_alginate: float = 0.0             # [kg/m³] Node F1-a; 0 = not-alginate
    phi_cellulose_0: float = 0.0        # [-] Node F1-b; initial cellulose vol. fraction
    solvent_system: str = ""            # [Node F1-b] cellulose solvent preset key; "" = skip
    phi_PLGA_0: float = 0.0             # [-] Node F1-c; initial PLGA vol. fraction in droplet
    plga_grade: str = ""                # [Node F1-c] PLGA grade preset key ("" = skip)
    c_span80: float = 20.0              # [kg/m³] (~2% w/v)
    c_genipin: float = 2.0              # [mol/m³] (~2 mM)
    c_Ca_bath: float = 100.0            # [mol/m³] Node F1-a; external CaCl₂ bath
    T_oil: float = 363.15              # [K] (90°C)
    cooling_rate: float = 0.167         # [K/s] (~10°C/min)
    T_crosslink: float = 310.15         # [K] (37°C)
    t_crosslink: float = 86400.0        # [s] (24 hours)
    phi_d: float = 0.05                 # [-] dispersed phase volume fraction
    pH: float = 7.0                     # [-] reaction pH (Node 31; EDC/NHS + future pH-dependent chemistries)
    # ── Stirred-vessel volumetric fields ──
    c_span80_vol_pct: float = 1.5       # [% v/v] Span-80 in paraffin oil
    v_oil_span80_mL: float = 300.0      # [mL] volume of oil + Span-80 mixture
    v_polysaccharide_mL: float = 200.0  # [mL] max volume of polysaccharide solution

    @property
    def agarose_fraction(self) -> float:
        """Agarose mass fraction in polymer blend."""
        total = self.c_agarose + self.c_chitosan
        if total == 0:
            return 0.0
        return self.c_agarose / total

    @property
    def total_polymer(self) -> float:
        """Total polymer concentration [kg/m³]."""
        return self.c_agarose + self.c_chitosan

    @property
    def c_span80_from_vol_pct(self) -> float:
        """Span-80 concentration [kg/m^3] derived from volumetric %.

        c_span80 = rho_span80 * (vol_pct / 100)
        where rho_span80 ~ 986 kg/m^3 (Sigma-Aldrich).
        """
        RHO_SPAN80 = 986.0  # [kg/m^3]
        return RHO_SPAN80 * self.c_span80_vol_pct / 100.0

    @property
    def total_working_volume_mL(self) -> float:
        """Total working liquid volume [mL]."""
        return self.v_oil_span80_mL + self.v_polysaccharide_mL

    @property
    def phi_d_from_volumes(self) -> float:
        """Dispersed-phase volume fraction computed from volumetric quantities.

        In W/O emulsion: polysaccharide solution is the dispersed phase.
        """
        total = self.total_working_volume_mL
        if total <= 0:
            return 0.0
        return self.v_polysaccharide_mL / total


@dataclass
class SolverSettings:
    """Numerical solver settings."""
    # Level 1
    l1_n_bins: int = 20
    l1_d_min: float = 1e-6             # [m]
    l1_d_max: float = 500e-6           # [m] (must exceed premix d32 by ≥3σ)
    l1_rtol: float = 1e-6
    l1_atol: float = 1e-8
    # Level 2
    l2_n_r: int = 1000
    l2_n_grid: int = 128               # 2D grid side length (N×N)
    l2_dt_initial: float = 1e-4        # [s]
    l2_dt_max: float = 1.0             # [s]
    l2_arrest_exponent: float = 2.5
    # Level 3
    l3_method: str = "Radau"
    l3_rtol: float = 1e-8
    l3_atol: float = 1e-10
    # Level 1 adaptive convergence
    l1_t_max: float = 600.0            # [s] absolute max emulsification time for adaptive extensions
    l1_conv_tol: float = 0.01          # [-] relative d32 variation threshold for convergence
    l1_max_extensions: int = 2         # [-] max number of adaptive time extensions


@dataclass
class SimulationParameters:
    """Top-level parameter container for the full pipeline."""
    model_mode: ModelMode = ModelMode.HYBRID_COUPLED
    emulsification: EmulsificationParameters = field(default_factory=EmulsificationParameters)
    formulation: FormulationParameters = field(default_factory=FormulationParameters)
    solver: SolverSettings = field(default_factory=SolverSettings)
    run_id: str = ""
    notes: str = ""

    def validate(self) -> list[str]:
        """Validate parameters. Returns list of error messages (empty = valid)."""
        errors = []
        e = self.emulsification
        f = self.formulation
        s = self.solver

        # ── Common checks ──
        if e.rpm <= 0:
            errors.append(f"RPM must be positive, got {e.rpm}")
        if e.t_emulsification <= 0:
            errors.append("Emulsification time must be positive")
        if f.c_agarose < 0:
            errors.append("c_agarose must be non-negative")
        if f.c_chitosan < 0:
            errors.append("c_chitosan must be non-negative")
        if f.c_span80 < 0:
            errors.append("c_span80 must be non-negative")
        if f.c_genipin < 0:
            errors.append("c_genipin must be non-negative")
        if f.T_oil <= 0:
            errors.append(f"T_oil must be positive (Kelvin), got {f.T_oil}")
        if f.T_crosslink <= 0:
            errors.append(f"T_crosslink must be positive (Kelvin), got {f.T_crosslink}")
        if f.cooling_rate <= 0:
            errors.append("cooling_rate must be positive")
        if f.t_crosslink <= 0:
            errors.append("t_crosslink must be positive")
        if not (0 < f.phi_d < 1):
            errors.append(f"phi_d must be in (0, 1), got {f.phi_d}")
        if s.l1_n_bins < 5:
            errors.append("l1_n_bins must be >= 5")
        if s.l1_d_min <= 0:
            errors.append("l1_d_min must be positive")
        if s.l1_d_max <= s.l1_d_min:
            errors.append("l1_d_max must exceed l1_d_min")

        # ── Mode-specific checks ──
        if e.mode == "rotor_stator_legacy":
            m = e.mixer
            if m.gap_width <= 0:
                errors.append("Mixer gap_width must be positive")
            if m.tank_volume <= 0:
                errors.append("Mixer tank_volume must be positive")
            if m.rotor_diameter <= 0:
                errors.append("Mixer rotor_diameter must be positive")

        elif e.mode == "stirred_vessel":
            # Stirrer must be populated
            if e.stirrer is None:
                errors.append("Stirred-vessel mode requires stirrer")
            else:
                if e.rpm > e.stirrer.max_rpm:
                    errors.append(
                        f"RPM {e.rpm} exceeds stirrer max {e.stirrer.max_rpm}"
                    )
                if e.stirrer.impeller_diameter <= 0:
                    errors.append("Stirrer impeller_diameter must be positive")
            # Vessel must be populated
            if e.vessel is None:
                errors.append("Stirred-vessel mode requires vessel")
            else:
                v = e.vessel
                if not (v.working_volume_min <= v.working_volume <= v.working_volume_max):
                    errors.append(
                        f"Working volume {v.working_volume*1e6:.0f} mL outside "
                        f"range [{v.working_volume_min*1e6:.0f}, "
                        f"{v.working_volume_max*1e6:.0f}] mL"
                    )
                if v.liquid_height > v.height:
                    errors.append(
                        f"Liquid height {v.liquid_height*1000:.1f} mm exceeds "
                        f"vessel height {v.height*1000:.0f} mm"
                    )
            # Heating-vessel compatibility
            if e.heating is not None and e.vessel is not None:
                h_val = getattr(e.heating.strategy, 'value', str(e.heating.strategy))
                v_val = getattr(e.vessel.vessel_type, 'value', str(e.vessel.vessel_type))
                if h_val == "flat_plate" and v_val != "glass_beaker":
                    errors.append(
                        "Flat-plate heating is only compatible with glass beaker"
                    )
                if h_val == "hot_water_jacket" and v_val != "jacketed_vessel":
                    errors.append(
                        "Hot-water jacket requires jacketed vessel"
                    )
            # Volumetric consistency
            total_mL = f.total_working_volume_mL
            if e.vessel is not None:
                vessel_mL = e.vessel.working_volume * 1e6
                if abs(total_mL - vessel_mL) > 1.0:  # tolerance 1 mL
                    errors.append(
                        f"Formulation volumes ({total_mL:.0f} mL) != "
                        f"vessel working volume ({vessel_mL:.0f} mL)"
                    )
        else:
            errors.append(
                f"Unknown emulsification mode: '{e.mode}'. "
                f"Use 'rotor_stator_legacy' or 'stirred_vessel'."
            )

        return errors

    def to_optimization_vector(self) -> np.ndarray:
        """Flatten to 7D vector for BoTorch."""
        return np.array([
            self.emulsification.rpm,
            self.formulation.c_span80,
            self.formulation.agarose_fraction,
            self.formulation.T_oil,
            self.formulation.cooling_rate,
            self.formulation.c_genipin,
            self.formulation.t_crosslink,
        ])

    @classmethod
    def from_optimization_vector(
        cls, x: np.ndarray, template: SimulationParameters
    ) -> SimulationParameters:
        """Reconstruct from 7D vector + template for fixed params."""
        import copy
        params = copy.deepcopy(template)
        params.emulsification.rpm = float(x[0])
        params.formulation.c_span80 = float(x[1])
        # Reconstruct agarose/chitosan from fraction and total
        frac = float(x[2])
        total = params.formulation.total_polymer
        params.formulation.c_agarose = frac * total
        params.formulation.c_chitosan = (1.0 - frac) * total
        params.formulation.T_oil = float(x[3])
        params.formulation.cooling_rate = float(x[4])
        params.formulation.c_genipin = float(x[5])
        params.formulation.t_crosslink = float(x[6])
        return params


# ─── Material Properties ─────────────────────────────────────────────────

@dataclass
class PropertyValue:
    """A single material property with metadata."""
    value: float
    unit: str
    uncertainty: float = 0.0
    source: str = ""
    T_ref: float = 298.15              # [K]


@dataclass
class MaterialProperties:
    """Aggregated material properties for the simulation."""
    # Oil phase
    rho_oil: float = 850.0              # [kg/m³] at 20°C reference (interpolated to T_oil by PropertyDatabase)
    mu_oil: float = 0.005               # [Pa·s] at 90°C reference (interpolated to T_oil by PropertyDatabase)

    # Aqueous / dispersed phase
    rho_aq: float = 1020.0              # [kg/m³]
    mu_d: float = 1.0                   # [Pa·s] dispersed phase viscosity at T_oil

    # Interfacial
    sigma: float = 5.0e-3               # [N/m] interfacial tension with Span-80

    # Thermodynamic
    chi_0: float = 0.497                # Flory-Huggins χ at reference T
    chi_T_coeffs: tuple = (515.5, -0.720)  # (A, B) for χ(T) = A/T + B; spinodal ~325 K for Np=10
    kappa_CH: float = 5.0e-12           # [J/m] Cahn-Hilliard gradient coefficient
    M_0: float = 1.0e-9                # [m⁵/(J·s)] bare mobility (calibrated for 50-100 nm coarsening)

    # Gelation
    T_gel: float = 311.15               # [K] (~38°C)
    k_gel_0: float = 1.0                # [1/s] Avrami rate prefactor
    n_avrami: float = 2.5               # [-]
    gel_arrest_exponent: float = 2.5    # β

    # Crosslinking
    k_xlink_0: float = 2806.0          # [m³/(mol·s)] calibrated: k(37°C)=5e-3 L/(mol·s)
    E_a_xlink: float = 52000.0         # [J/mol]
    DDA: float = 0.90                   # degree of deacetylation
    M_GlcN: float = 161.16             # [g/mol] glucosamine molar mass
    M_genipin: float = 226.23          # [g/mol]

    # Agarose gel
    G_agarose_prefactor: float = 3000.0  # [Pa] at 1% w/v
    G_agarose_exponent: float = 2.2      # power law exponent

    # Crosslinking bridge efficiency
    f_bridge: float = 0.4              # fraction of genipin reactions producing elastically active crosslinks

    # IPN coupling
    eta_coupling: float = -0.15        # IPN coupling coefficient

    # Chitosan viscosity
    eta_intr_chit: float = 800.0       # [mL/g] intrinsic viscosity of chitosan

    # Shear-thinning (Cross model) — set to enable non-Newtonian dispersed phase
    cross_mu_0: float = 0.0       # [Pa.s] zero-shear viscosity (0 = disabled, use constant mu_d)
    cross_mu_inf: float = 0.001   # [Pa.s] infinite-shear viscosity
    cross_K: float = 0.01         # [s] relaxation time
    cross_n: float = 0.6          # [-] power-law index

    # Breakage kernel
    breakage_C3: float = 0.0           # [-] Alopaeus viscous correction. Default 0.0: for aqueous
                                        # dispersed phases (mu_d ~ 0.01-0.1 Pa·s at emulsification T),
                                        # the viscous correction is physically negligible. The Vi cap
                                        # in breakage_rate_alopaeus() prevents exponential blowout for
                                        # any C3 > 0. Set C3 > 0 only for highly viscous dispersed
                                        # phases (mu_d > 1 Pa·s) after calibration.

    # Network / pore
    r_fiber: float = 1.5e-9            # [m] agarose fiber radius

    # Surface carboxyl (Node 31): introduced by prior M2 steps (e.g.
    # succinylation of chitosan NH2 -> NH-CO-CH2-CH2-COOH). Zero on
    # native chitosan/agarose; non-zero gates L3 EDC/NHS to run the
    # mechanistic path instead of falling back to QUALITATIVE_TREND.
    surface_cooh_concentration: float = 0.0  # [mol/m^3] grafted COOH sites

    # Platform / polymer family (Node F1-a, v8.0): drives L2 solver
    # dispatch. Default keeps legacy chitosan/agarose behaviour.
    polymer_family: PolymerFamily = PolymerFamily.AGAROSE_CHITOSAN

    # Alginate-specific defaults (F1-a). Harmless for other families.
    f_guluronate: float = 0.5          # [-] guluronate fraction (G-block)
    D_Ca: float = 1.0e-9               # [m²/s] Ca²⁺ diffusivity in alginate gel
    k_bind_Ca: float = 1.0e3           # [M⁻²·s⁻¹] Ca²⁺ + 2 COO⁻ binding rate
    K_alg_modulus: float = 30.0e3      # [Pa] modulus prefactor (Kong 2004)
    n_alg_modulus: float = 2.0         # [-] modulus exponent

    # Cellulose-specific defaults (F1-b Phase 1, NaOH/urea). Harmless for
    # other families. See docs/f1b_cellulose_nips_protocol.md §5.
    N_p_cellulose: float = 370.0           # [-] degree of polymerisation (M_n/M_AGU)
    chi_PS_cellulose: float = 0.45         # [-] χ polymer-solvent
    chi_PN_cellulose: float = 0.85         # [-] χ polymer-nonsolvent (water)
    chi_SN_cellulose: float = 0.30         # [-] χ solvent-nonsolvent
    D_solvent_cellulose: float = 5.0e-11   # [m²/s] solvent self-diffusion in gel
    D_nonsolvent_cellulose: float = 1.0e-10  # [m²/s] water self-diffusion in gel
    kappa_CH_cellulose: float = 1.0e-17    # [J·m⁻¹] Cahn-Hilliard gradient energy coef
    K_cell_modulus: float = 5.0e5          # [Pa] modulus prefactor (Zhang 2020)
    alpha_cell_modulus: float = 2.25       # [-] modulus exponent (entangled regime)

    # PLGA-specific defaults (F1-c Phase 1, PLGA 50:50). Harmless for
    # other families. See docs/f1c_plga_protocol.md §5.
    D_DCM_plga: float = 1.0e-9             # [m²/s] DCM effective diffusivity in PLGA/DCM
    phi_DCM_eq: float = 0.005              # [-] bath-equilibrium DCM fraction (Henry)
    G_glassy_plga: float = 7.0e8           # [Pa] PLGA glassy-state shear modulus
    n_plga_modulus: float = 2.0            # [-] Gibson-Ashby exponent


# ─── Model Evidence and Provenance ────────────────────────────────────────

@dataclass
class ModelManifest:
    """Provenance record for a single model/solver used in a pipeline level.

    Attached to each result object to make evidence tier, assumptions, and
    validity domain traceable without reading source code.
    """
    model_name: str                         # e.g. "L1.PBE.FixedPivot.AlopaeusCT"
    evidence_tier: ModelEvidenceTier = ModelEvidenceTier.SEMI_QUANTITATIVE
    valid_domain: dict = field(default_factory=dict)  # {"Re": (100, 1e6), ...}
    calibration_ref: str = ""               # "" or CalibrationEntry ID
    assumptions: list[str] = field(default_factory=list)
    diagnostics: dict = field(default_factory=dict)   # {"thiele": 0.3, "mass_conservation": 0.999}


@dataclass
class RunReport:
    """Structured report for a complete pipeline run.

    Collects model manifests from all levels plus trust assessment and
    solver diagnostics. Suitable for JSON export and lab notebook attachment.
    """
    model_graph: list[ModelManifest] = field(default_factory=list)
    trust_level: str = ""                   # "TRUSTWORTHY" | "CAUTION" | "UNRELIABLE"
    trust_warnings: list[str] = field(default_factory=list)
    trust_blockers: list[str] = field(default_factory=list)
    diagnostics: dict = field(default_factory=dict)   # aggregated diagnostics
    min_evidence_tier: str = ""             # weakest tier across all levels

    def compute_min_tier(self) -> ModelEvidenceTier:
        """Return the weakest evidence tier across all models in the graph."""
        if not self.model_graph:
            return ModelEvidenceTier.UNSUPPORTED
        # Tier ordering: validated > calibrated > semi > qualitative > unsupported
        _ORDER = list(ModelEvidenceTier)
        # Compare by .value (string), not identity, so manifests created
        # against a stale class object (e.g. after importlib.reload of this
        # module by the Streamlit app) still match the current enum members.
        _order_values = [t.value for t in _ORDER]
        worst_idx = 0
        for m in self.model_graph:
            tier = m.evidence_tier
            tier_value = getattr(tier, "value", tier)
            try:
                idx = _order_values.index(tier_value)
            except ValueError:
                idx = _order_values.index(ModelEvidenceTier.UNSUPPORTED.value)
            if idx > worst_idx:
                worst_idx = idx
        return _ORDER[worst_idx]


# ─── Cross-cutting Run Context (Node 7, v6.1) ─────────────────────────────


@dataclass
class RunContext:
    """Cross-cutting inputs to a pipeline run.

    Carries optional services that the orchestrator should consult before
    invoking the per-level solvers. Designed to be backward-compatible:
    ``run_context=None`` reproduces the v6.0 behaviour exactly.

    v6.1 scope (Node 7): only ``calibration_store`` is wired. Future
    additions (random_seed routing, output_policy, uncertainty_spec) follow
    the same opt-in pattern — pass them in to enable; omit to keep current
    behaviour. The richer ``ParameterProvider`` with full provenance
    tracking (ResolvedParameter source/uncertainty) is deferred to v7.0.
    """
    calibration_store: Optional[Any] = None
    run_id: str = ""
    notes: str = ""


# ─── Result Structures ────────────────────────────────────────────────────

@dataclass
class EmulsificationResult:
    """Output of Level 1: PBE solver."""
    d_bins: np.ndarray                  # [m] pivot diameters (N_bins,)
    n_d: np.ndarray                     # [#/m³] number density (N_bins,)
    d32: float                          # [m] Sauter mean diameter
    d43: float                          # [m] volume-weighted mean
    d10: float                          # [m] 10th percentile
    d50: float                          # [m] median
    d90: float                          # [m] 90th percentile
    span: float                         # [-] (d90 - d10) / d50
    total_volume_fraction: float        # [-]
    converged: bool
    d_mode: float = 0.0                 # [m] modal diameter (volume-weighted)
    t_history: Optional[np.ndarray] = None
    n_d_history: Optional[np.ndarray] = None
    t_converged: Optional[float] = None # [s] time at which d32 convergence was first achieved (None if never)
    n_extensions: int = 0               # [-] number of adaptive extensions performed
    model_manifest: Optional[ModelManifest] = None  # v6.1: evidence provenance


@dataclass
class GelationTimingResult:
    """Output of Level 2a: Gelation timing and arrest.

    Describes WHEN and HOW FAST gelation occurs, separate from
    what microstructure forms (Level 2b).
    """
    T_history: np.ndarray           # [K] temperature history (N_t,)
    t_gel_onset: float              # [s] time at which gelation starts
    alpha_final: float              # [-] final gelation fraction (0-1)
    mobility_arrest_factor: float   # [-] how much mobility is reduced at arrest
    cooling_rate_effective: float   # [K/s] effective cooling rate (may be size-dependent)


@dataclass
class GelationResult:
    """Output of Level 2: Phase-field solver.

    Supports both 1D radial (legacy) and 2D Cartesian solvers.
    For 1D: r_grid is (N_r,), phi_field is (N_r,).
    For 2D: r_grid is (N,) coordinate array, phi_field is (N, N).
    """
    r_grid: np.ndarray                  # [m] (N_r,) or (N,) for 2D coords
    phi_field: np.ndarray               # [-] (N_r,) or (N, N)
    pore_size_mean: float               # [m]
    pore_size_std: float                # [m]
    pore_size_distribution: np.ndarray  # [m]
    porosity: float                     # [-]
    alpha_final: float                  # [-]
    char_wavelength: float              # [m]
    T_history: Optional[np.ndarray] = None
    phi_snapshots: Optional[np.ndarray] = None
    L_domain: float = 0.0              # [m] domain side length (2D solver)
    grid_spacing: float = 0.0          # [m] uniform grid spacing
    # Morphology descriptors (v3.0)
    bicontinuous_score: float = 0.5   # [-] 0-1
    anisotropy: float = 0.0           # [-] 0-1
    connectivity: float = 1.0         # [-] 0-1, fraction of pore space connected
    chord_skewness: float = 0.0       # [-] skewness of chord length distribution
    model_tier: str = "unknown"       # "empirical_calibrated" or "mechanistic" — identifies L2 model type
    model_manifest: Optional[ModelManifest] = None  # v6.1: evidence provenance


@dataclass
class NetworkTypeMetadata:
    """Describes what network was formed by crosslinking."""
    solver_family: str = "amine_covalent"
    network_target: str = "chitosan"      # "chitosan", "agarose", "independent", "mixed"
    bond_type: str = "covalent"           # "covalent", "ionic", "reversible"
    is_true_second_network: bool = True   # True for IPN, False for reinforcement
    eta_coupling_recommended: float = -0.15  # [-] per-chemistry IPN coupling coefficient from CrosslinkerProfile


@dataclass
class CrosslinkingResult:
    """Output of Level 3: ODE kinetics."""
    t_array: np.ndarray                 # [s] (N_t,)
    X_array: np.ndarray                 # [mol/m³] (N_t,)
    nu_e_array: np.ndarray              # [1/m³] (N_t,)
    Mc_array: np.ndarray                # [g/mol] (N_t,)
    xi_array: np.ndarray                # [m] (N_t,)
    G_chitosan_array: np.ndarray        # [Pa] (N_t,)
    p_final: float                      # [-]
    nu_e_final: float                   # [1/m³]
    Mc_final: float                     # [g/mol]
    xi_final: float                     # [m]
    G_chitosan_final: float             # [Pa]
    network_metadata: Optional[NetworkTypeMetadata] = None
    model_manifest: Optional[ModelManifest] = None  # v6.1: evidence provenance
    # v6.1: solver diagnostics
    thiele_modulus: float = 0.0         # [-] Phi (0 = not computed)
    regime: str = "unknown"             # "reaction_limited", "borderline", "diffusion_limited", "unknown"
    stoichiometric_ceiling: float = 1.0 # [-] max conversion given reagent/reactive-group ratio
    residual_reactive_groups: float = 0.0  # [mol/m3] unreacted -NH2 or -OH after crosslinking


@dataclass
class MechanicalResult:
    """Output of Level 4: Property prediction."""
    G_agarose: float                    # [Pa]
    G_chitosan: float                   # [Pa]
    G_DN: float                         # [Pa]
    E_star: float                       # [Pa]
    delta_array: np.ndarray             # [m]
    F_array: np.ndarray                 # [N]
    rh_array: np.ndarray                # [m]
    Kav_array: np.ndarray               # [-]
    pore_size_mean: float               # [m]
    xi_mesh: float                      # [m]
    model_used: str = "phenomenological"  # which modulus model produced G_DN
    G_DN_lower: float = 0.0             # [Pa] Single-phase composite reference (HS bounds, not applicable to IPN)
    G_DN_upper: float = 0.0             # [Pa] Single-phase composite reference (HS bounds, not applicable to IPN)
    model_manifest: Optional[ModelManifest] = None  # v6.1: evidence provenance
    network_type: str = "unknown"       # v6.1: "true_IPN", "semi_IPN", "independent_network", "ionic_reinforced", "unknown"


@dataclass
class FullResult:
    """Complete pipeline output."""
    parameters: SimulationParameters
    emulsification: EmulsificationResult
    gelation: GelationResult
    crosslinking: CrosslinkingResult
    mechanical: MechanicalResult
    gelation_timing: Optional[GelationTimingResult] = None
    run_report: Optional[RunReport] = None  # v6.1: model evidence + diagnostics

    def objective_vector(self) -> np.ndarray:
        """Compute the 3 objective values for optimization.

        Delegates to the canonical implementation in
        ``emulsim.optimization.objectives.compute_objectives``.
        """
        from emulsim.optimization.objectives import compute_objectives
        return compute_objectives(self)


@dataclass
class M1ExportContract:
    """Stable interface between Module 1 (fabrication) and Module 2 (functionalization).

    Contains all Module 1 outputs needed by Module 2, with explicit units,
    types, and trust metadata. This is the ONLY object Module 2 should read
    from Module 1 — it decouples the two modules.
    """
    # ── Geometry (from L1) ──
    bead_radius: float              # [m] from d50/2
    bead_d32: float                 # [m] Sauter mean diameter
    bead_d50: float                 # [m] median diameter

    # ── Pore structure (from L2) ──
    pore_size_mean: float           # [m] mean pore diameter
    pore_size_std: float            # [m] standard deviation of pore sizes
    porosity: float                 # [-] particle porosity (fraction void)
    l2_model_tier: str              # "empirical_calibrated" or "mechanistic"

    # ── Network structure (from L3) ──
    mesh_size_xi: float             # [m] crosslinked network mesh size
    p_final: float                  # [-] crosslinking conversion fraction
    primary_crosslinker: str        # e.g., "genipin", "dvs"

    # ── Residual reactive groups (from L3 + formulation) ──
    nh2_bulk_concentration: float   # [mol/m^3] residual NH2 after primary crosslinking
    oh_bulk_concentration: float    # [mol/m^3] total agarose OH

    # ── Mechanical properties (from L4) ──
    G_DN: float                     # [Pa] double-network shear modulus
    E_star: float                   # [Pa] effective Young's modulus
    model_used: str                 # "phenomenological" or "flory_rehner_affine"

    # ── Formulation (source parameters) ──
    c_agarose: float                # [kg/m^3] agarose concentration
    c_chitosan: float               # [kg/m^3] chitosan concentration
    DDA: float                      # [-] degree of deacetylation

    # ── Trust and uncertainty ──
    trust_level: str                # "TRUSTWORTHY", "CAUTION", "UNRELIABLE"
    trust_warnings: list[str] = field(default_factory=list)
    uncertainty_notes: str = ""     # Human-readable uncertainty description

    def validate_units(self) -> list[str]:
        """Node 10 (F11): boundary-level unit/range sanity checks.

        Returns a list of violation messages (empty = pass). Catches the
        common silent-failure modes of someone shoving a value in the wrong
        unit through the M1->M2 contract: e.g. a bead_radius in microns
        instead of meters, a pH in mol/m^3, a concentration in g/L.

        These are *order-of-magnitude* guards, not scientific calibration —
        a reasonable polysaccharide microsphere lives entirely inside every
        bound checked here. Failing one of these means the value is
        physically impossible for our system, not just unusual.
        """
        violations: list[str] = []

        # Geometry: SI metres. A µm number passed as m would land at 1e-6;
        # a cm value would land at 1e-2.
        if not (1e-7 <= self.bead_radius <= 5e-3):
            violations.append(
                f"bead_radius={self.bead_radius:g} m outside [1e-7, 5e-3]; "
                "wrong unit (µm vs m)?"
            )
        if not (1e-7 <= self.bead_d50 <= 1e-2):
            violations.append(
                f"bead_d50={self.bead_d50:g} m outside [1e-7, 1e-2]; wrong unit?"
            )
        if not (1e-9 <= self.pore_size_mean <= 1e-5):
            violations.append(
                f"pore_size_mean={self.pore_size_mean:g} m outside [1e-9, 1e-5]; "
                "did you pass nm as m?"
            )
        if self.pore_size_std < 0:
            violations.append(
                f"pore_size_std={self.pore_size_std:g} must be non-negative."
            )
        if not (0.0 <= self.porosity <= 1.0):
            violations.append(
                f"porosity={self.porosity:g} outside [0, 1]; must be a fraction."
            )
        if not (1e-10 <= self.mesh_size_xi <= 1e-5):
            violations.append(
                f"mesh_size_xi={self.mesh_size_xi:g} m outside [1e-10, 1e-5]."
            )
        if not (0.0 <= self.p_final <= 1.0):
            violations.append(
                f"p_final={self.p_final:g} outside [0, 1]; conversion fraction."
            )

        # Concentrations: mol/m^3. Native chitosan amine sits ~50-150 mol/m^3
        # for a 1-3% w/v solution. Pushing 1e5 means kg/m^3 was passed
        # as mol/m^3.
        if not (0.0 <= self.nh2_bulk_concentration <= 1e4):
            violations.append(
                f"nh2_bulk_concentration={self.nh2_bulk_concentration:g} mol/m^3 "
                "outside [0, 1e4]; wrong unit (mM vs mol/m^3)?"
            )
        if not (0.0 <= self.oh_bulk_concentration <= 1e5):
            violations.append(
                f"oh_bulk_concentration={self.oh_bulk_concentration:g} mol/m^3 "
                "outside [0, 1e5]."
            )

        # Mechanical: Pa. Gel modulus 1e3-1e6 Pa. 1 GPa would be a steel
        # bead; 1e-3 Pa would be water.
        if not (1.0 <= self.G_DN <= 1e9):
            violations.append(
                f"G_DN={self.G_DN:g} Pa outside [1, 1e9]; wrong unit (kPa vs Pa)?"
            )
        if not (1.0 <= self.E_star <= 1e10):
            violations.append(
                f"E_star={self.E_star:g} Pa outside [1, 1e10]."
            )

        # Formulation: kg/m^3 (i.e. g/L). 0-200 kg/m^3 covers every
        # processable polysaccharide gel. A value passed in % w/v would
        # land in [0, 20] which we cannot distinguish; only flag obvious
        # nonsense (negative, >300 kg/m^3 = 30% which is unprocessable).
        for name, val in (("c_agarose", self.c_agarose),
                          ("c_chitosan", self.c_chitosan)):
            if not (0.0 <= val <= 300.0):
                violations.append(
                    f"{name}={val:g} kg/m^3 outside [0, 300]."
                )

        if not (0.0 <= self.DDA <= 1.0):
            violations.append(
                f"DDA={self.DDA:g} outside [0, 1]; degree of deacetylation."
            )

        return violations


@dataclass
class OptimizationState:
    """State of the Bayesian optimization campaign."""
    X_observed: np.ndarray              # (N_eval, 7)
    Y_observed: np.ndarray              # (N_eval, 3)
    pareto_X: np.ndarray
    pareto_Y: np.ndarray
    iteration: int
    hypervolume: float
    hypervolume_history: list = field(default_factory=list)
    converged: bool = False
    gp_state: Optional[dict] = None
    # Node 6 (v6.1): per-Pareto-candidate weakest evidence tier (string form,
    # JSON-friendly). Empty when no run_report was attached. Length matches
    # pareto_X.shape[0].
    pareto_evidence_tiers: list[str] = field(default_factory=list)
