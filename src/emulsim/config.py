"""Configuration loading and validation for the simulation.

Supports two config formats:
  - Legacy: [emulsification.mixer] with rotor-stator fields
  - Stirred-vessel: [vessel], [stirrer], [heating] top-level sections
The ``mode`` field in [emulsification] selects the code path.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

from .datatypes import (
    EmulsificationParameters,
    FormulationParameters,
    HeatingConfig,
    HeatingStrategy,
    KernelConfig,
    MaterialProperties,
    MixerGeometry,
    SimulationParameters,
    SolverSettings,
    StirrerGeometry,
    StirrerType,
    VesselGeometry,
    VesselType,
)


def _warn_unknown_keys(raw: dict, known: set[str], section: str) -> None:
    """Warn about unrecognized keys in a TOML section."""
    unknown = set(raw.keys()) - known
    if unknown:
        logger.warning(
            "[%s] unrecognized config keys (ignored): %s",
            section, ", ".join(sorted(unknown)),
        )


_VESSEL_KEYS = {
    "type", "working_volume", "inner_diameter", "wall_thickness",
    "height", "material", "working_volume_min", "working_volume_max",
}
_STIRRER_KEYS = {
    "type", "impeller_diameter", "shaft_diameter", "blade_count",
    "blade_angle", "blade_thickness", "blade_height", "blade_length",
    "power_number", "max_rpm", "stator_diameter", "gap_width",
    "wall_height", "wall_thickness", "perforation_diameter",
    "dissipation_ratio",
}
_HEATING_KEYS = {
    "strategy", "heater_temperature", "T_initial", "T_final",
    "cooldown_time", "jacket_water_temperature",
}


def _build_vessel(raw: dict) -> VesselGeometry:
    """Construct VesselGeometry from TOML [vessel] section."""
    _warn_unknown_keys(raw, _VESSEL_KEYS, "vessel")
    vtype_str = raw.get("type", "glass_beaker")
    vtype = VesselType(vtype_str)
    wv = raw.get("working_volume", 0.0005)
    if vtype == VesselType.GLASS_BEAKER:
        vessel = VesselGeometry.glass_beaker(working_volume=wv)
    elif vtype == VesselType.JACKETED_VESSEL:
        vessel = VesselGeometry.jacketed_vessel(working_volume=wv)
    else:
        vessel = VesselGeometry(working_volume=wv)
    # Override any explicitly set fields
    for key in ("inner_diameter", "wall_thickness", "height", "material",
                "working_volume_min", "working_volume_max"):
        if key in raw:
            setattr(vessel, key, raw[key])
    return vessel


def _build_stirrer(raw: dict) -> StirrerGeometry:
    """Construct StirrerGeometry from TOML [stirrer] section."""
    _warn_unknown_keys(raw, _STIRRER_KEYS, "stirrer")
    stype_str = raw.get("type", "pitched_blade")
    stype = StirrerType(stype_str)
    if stype == StirrerType.PITCHED_BLADE:
        stirrer = StirrerGeometry.pitched_blade_A()
    elif stype == StirrerType.ROTOR_STATOR_SMALL:
        stirrer = StirrerGeometry.rotor_stator_B()
    elif stype == StirrerType.ROTOR_STATOR_LEGACY:
        stirrer = StirrerGeometry.rotor_stator_legacy()
    else:
        stirrer = StirrerGeometry()
    # Override any explicitly set fields
    for key in ("impeller_diameter", "shaft_diameter", "blade_count",
                "blade_angle", "blade_thickness", "blade_height",
                "blade_length", "power_number", "max_rpm",
                "stator_diameter", "gap_width", "wall_height",
                "wall_thickness", "perforation_diameter", "dissipation_ratio"):
        if key in raw:
            setattr(stirrer, key, raw[key])
    return stirrer


def _build_heating(raw: dict) -> HeatingConfig:
    """Construct HeatingConfig from TOML [heating] section."""
    _warn_unknown_keys(raw, _HEATING_KEYS, "heating")
    strategy_str = raw.get("strategy", "flat_plate")
    strategy = HeatingStrategy(strategy_str)
    if strategy == HeatingStrategy.FLAT_PLATE:
        heating = HeatingConfig.flat_plate()
    elif strategy == HeatingStrategy.HOT_WATER_JACKET:
        heating = HeatingConfig.hot_water_jacket()
    elif strategy == HeatingStrategy.ISOTHERMAL:
        T = raw.get("T_initial", 363.15)
        heating = HeatingConfig.isothermal(T)
    else:
        heating = HeatingConfig()
    # Override any explicitly set fields
    for key in ("heater_temperature", "T_initial", "T_final",
                "cooldown_time", "jacket_water_temperature"):
        if key in raw:
            setattr(heating, key, raw[key])
    return heating


def load_config(path: Path) -> SimulationParameters:
    """Load simulation parameters from a TOML file.

    Detects mode from the ``[emulsification].mode`` field:
      - ``"rotor_stator_legacy"`` (default): reads [emulsification.mixer]
      - ``"stirred_vessel"``: reads [vessel], [stirrer], [heating]
    """
    import tomllib

    with open(path, "rb") as f:
        raw = tomllib.load(f)

    emul_raw = raw.get("emulsification", {})
    mode = emul_raw.get("mode", "rotor_stator_legacy")

    _VALID_MODES = {"rotor_stator_legacy", "stirred_vessel"}
    if mode not in _VALID_MODES:
        raise ValueError(
            f"Unknown emulsification mode '{mode}' in config. "
            f"Valid modes: {sorted(_VALID_MODES)}"
        )

    if mode == "stirred_vessel":
        # ── New stirred-vessel config path ──
        vessel = _build_vessel(raw.get("vessel", {}))
        stirrer = _build_stirrer(raw.get("stirrer", {}))
        heating = _build_heating(raw.get("heating", {}))
        kernels = KernelConfig.for_stirrer_type(stirrer.stirrer_type)

        emul_fields = {k: v for k, v in emul_raw.items()
                       if k not in ("mixer", "mode")}
        emul = EmulsificationParameters(
            mode="stirred_vessel",
            vessel=vessel,
            stirrer=stirrer,
            heating=heating,
            kernels=kernels,
            **emul_fields,
        )
    else:
        # ── Legacy rotor-stator config path ──
        mixer_data = emul_raw.get("mixer", {})
        mixer = MixerGeometry(**mixer_data) if mixer_data else MixerGeometry()

        emul_fields = {k: v for k, v in emul_raw.items()
                       if k not in ("mixer", "mode")}
        emul = EmulsificationParameters(
            mode="rotor_stator_legacy",
            mixer=mixer,
            **emul_fields,
        )

    # ── Formulation ──
    form_data = raw.get("formulation", {})
    form = FormulationParameters(**form_data) if form_data else FormulationParameters()

    # ── Solver ──
    solver_raw = raw.get("solver", {})
    solver_flat = {}
    for section_key, section_val in solver_raw.items():
        if isinstance(section_val, dict):
            prefix = section_key.replace("level", "l")
            for k, v in section_val.items():
                solver_flat[f"{prefix}_{k}"] = v
        else:
            solver_flat[section_key] = section_val
    solver = SolverSettings(**solver_flat) if solver_flat else SolverSettings()

    # ── Simulation metadata ──
    sim_data = raw.get("simulation", {})

    # Model mode
    from .datatypes import ModelMode
    mode_str = sim_data.get("model_mode", "hybrid_coupled")
    try:
        model_mode = ModelMode(mode_str)
    except ValueError:
        model_mode = ModelMode.HYBRID_COUPLED

    return SimulationParameters(
        model_mode=model_mode,
        emulsification=emul,
        formulation=form,
        solver=solver,
        run_id=sim_data.get("run_id", ""),
        notes=sim_data.get("notes", ""),
    )


def load_properties(path: Optional[Path] = None) -> MaterialProperties:
    """Load material properties from TOML or return defaults."""
    if path is None:
        return MaterialProperties()

    import tomllib

    with open(path, "rb") as f:
        raw = tomllib.load(f)

    props = {}
    # Flatten nested structure into MaterialProperties fields
    for section in raw.values():
        if isinstance(section, dict):
            for key, val in section.items():
                if isinstance(val, dict) and "value" in val:
                    props[key] = val["value"]
                elif not isinstance(val, dict):
                    props[key] = val

    return MaterialProperties(**{k: v for k, v in props.items() if hasattr(MaterialProperties, k)})
