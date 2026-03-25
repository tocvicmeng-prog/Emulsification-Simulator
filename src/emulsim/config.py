"""Configuration loading and validation for the simulation."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .datatypes import (
    EmulsificationParameters,
    FormulationParameters,
    MaterialProperties,
    MixerGeometry,
    SimulationParameters,
    SolverSettings,
)


def load_config(path: Path) -> SimulationParameters:
    """Load simulation parameters from a TOML file."""
    import tomllib

    with open(path, "rb") as f:
        raw = tomllib.load(f)

    mixer_data = raw.get("emulsification", {}).get("mixer", {})
    mixer = MixerGeometry(**mixer_data) if mixer_data else MixerGeometry()

    emul_data = {k: v for k, v in raw.get("emulsification", {}).items() if k != "mixer"}
    emul = EmulsificationParameters(mixer=mixer, **emul_data) if emul_data else EmulsificationParameters(mixer=mixer)

    form_data = raw.get("formulation", {})
    form = FormulationParameters(**form_data) if form_data else FormulationParameters()

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

    sim_data = raw.get("simulation", {})

    return SimulationParameters(
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
