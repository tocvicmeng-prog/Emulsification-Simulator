"""Calibration data types for user-supplied measured parameters.

v6.0-alpha: Typed schema with units, target molecule, validity domain,
and source reference (audit F2 requirement).
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CalibrationEntry:
    """A single calibration measurement for a reagent profile parameter.

    Captures not just the value, but the conditions under which it was
    measured, so the simulator can validate applicability.

    Attributes:
        profile_key: Reagent profile key (e.g., "protein_a_coupling").
        parameter_name: Parameter being calibrated (e.g., "q_max", "K_L",
            "activity_retention", "estimated_q_max").
        measured_value: The measured value in the specified units.
        units: SI or common units (e.g., "mg/mL", "mol/m3", "fraction").
        target_molecule: Identity of the target protein/analyte.
        temperature_C: Measurement temperature [Celsius].
        ph: Measurement pH.
        salt_concentration_M: Salt concentration [M] during measurement.
        salt_type: Salt identity (e.g., "NaCl", "(NH4)2SO4").
        measurement_type: How the value was measured.
        confidence: Data quality tier.
        source_reference: Literature reference or lab notebook ID.
        replicates: Number of independent replicates.
    """
    profile_key: str
    parameter_name: str
    measured_value: float
    units: str
    target_molecule: str = ""
    temperature_C: float = 25.0
    ph: float = 7.0
    salt_concentration_M: float = 0.0
    salt_type: str = ""
    measurement_type: str = ""    # "static_binding", "DBC10", "DBC5", "batch_uptake"
    confidence: str = "measured"  # "measured", "literature", "estimated"
    source_reference: str = ""
    replicates: int = 1

    def to_dict(self) -> dict:
        """Serialize to dict for JSON storage."""
        return {
            "profile_key": self.profile_key,
            "parameter_name": self.parameter_name,
            "measured_value": self.measured_value,
            "units": self.units,
            "target_molecule": self.target_molecule,
            "temperature_C": self.temperature_C,
            "ph": self.ph,
            "salt_concentration_M": self.salt_concentration_M,
            "salt_type": self.salt_type,
            "measurement_type": self.measurement_type,
            "confidence": self.confidence,
            "source_reference": self.source_reference,
            "replicates": self.replicates,
        }

    @classmethod
    def from_dict(cls, d: dict) -> CalibrationEntry:
        """Deserialize from dict."""
        return cls(
            profile_key=d["profile_key"],
            parameter_name=d["parameter_name"],
            measured_value=float(d["measured_value"]),
            units=d.get("units", ""),
            target_molecule=d.get("target_molecule", ""),
            temperature_C=float(d.get("temperature_C", 25.0)),
            ph=float(d.get("ph", 7.0)),
            salt_concentration_M=float(d.get("salt_concentration_M", 0.0)),
            salt_type=d.get("salt_type", ""),
            measurement_type=d.get("measurement_type", ""),
            confidence=d.get("confidence", "measured"),
            source_reference=d.get("source_reference", ""),
            replicates=int(d.get("replicates", 1)),
        )
