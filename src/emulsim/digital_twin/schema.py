"""Node F2 Phase 1: trace schema for the digital-twin replay harness.

A ``DigitalTwinTrace`` is a time-ordered sequence of scalar
``Observation`` records plus free-form metadata. Persistence is JSON
(trivial forward-compat via ``extra`` fields).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class Observation:
    """A single scalar process observation.

    Attributes
    ----------
    t : float
        Time since experiment start [s].
    name : str
        Observable identifier — the replay harness's
        ``observation_operator`` callback receives this to decide what
        slice of the ensemble state to project.
    value : float
        Measured scalar value.
    noise_std : float
        One-sigma observation noise. Used by ``enkf_update`` as
        ``R = noise_std**2``.
    """

    t: float
    name: str
    value: float
    noise_std: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "t": float(self.t),
            "name": str(self.name),
            "value": float(self.value),
            "noise_std": float(self.noise_std),
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> "Observation":
        return cls(
            t=float(raw["t"]),
            name=str(raw["name"]),
            value=float(raw["value"]),
            noise_std=float(raw["noise_std"]),
        )


@dataclass
class DigitalTwinTrace:
    """Recorded process trajectory plus metadata.

    ``observations`` is always time-sorted on load (defensive; the
    EnKF replay relies on monotone time).
    """

    trace_id: str
    process_description: str
    observations: list[Observation]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Defensive re-sort — replays of corrupted / manually-edited
        # traces still work.
        self.observations = sorted(self.observations, key=lambda o: o.t)

    def to_dict(self) -> dict[str, Any]:
        return {
            "trace_id": self.trace_id,
            "process_description": self.process_description,
            "observations": [o.to_dict() for o in self.observations],
            "metadata": dict(self.metadata),
        }


def load_trace(path: Path | str) -> DigitalTwinTrace:
    """Parse a JSON trace file into a :class:`DigitalTwinTrace`.

    Raises
    ------
    ValueError
        If the JSON is missing required keys.
    FileNotFoundError
        If ``path`` does not exist.
    """
    p = Path(path)
    with open(p, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return _from_dict(raw)


def _from_dict(raw: dict[str, Any]) -> DigitalTwinTrace:
    for key in ("trace_id", "process_description", "observations"):
        if key not in raw:
            raise ValueError(
                f"DigitalTwinTrace missing required key {key!r}"
            )
    obs_list = [Observation.from_dict(o) for o in raw["observations"]]
    return DigitalTwinTrace(
        trace_id=str(raw["trace_id"]),
        process_description=str(raw["process_description"]),
        observations=obs_list,
        metadata=dict(raw.get("metadata", {})),
    )


def save_trace(trace: DigitalTwinTrace, path: Path | str) -> None:
    """Serialise a :class:`DigitalTwinTrace` to JSON."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(trace.to_dict(), f, indent=2)
