"""Node 17 (v7.0, P2b): AssayRecord — public data model for wet-lab ingest.

This is the contract that lets a wet-lab assay (a triplicate IFT, a DSD-vs-RPM
sweep, a static binding capacity) flow into the platform as machine-readable
data and become the source of a CalibrationEntry. It closes consensus-plan
Rule 9 ("every assay should be ingestible as data") which was deferred from
v6.1.

Design goals (per v7.0 plan):
  - **Stable schema for v7.0**: this is now a public data contract.
    Wet-lab tooling and external scripts will reference it by field name.
    Breaking changes need a schema_version bump and migration notes.
  - **Self-describing units**: every measurement carries an explicit unit
    string. The downstream calibration fitter is responsible for SI
    coercion at the boundary (see ``CalibrationEntry.units`` -- same
    convention).
  - **Replicate-aware**: most useful assays produce a small N of replicate
    measurements with attached metadata. AssayRecord stores them as a list
    so the fitter can compute mean / std / SEM as needed.
  - **Provenance**: every record traces back to an instrument, an
    operator, and a notebook entry. Without this, a calibration run is
    not auditable.

Out of scope for v7.0 (push to v7.1+):
  - Validated JSON-Schema with a vendored validator
  - UI upload widgets (the data flow is JSON-on-disk for v7.0)
  - Multi-batch aggregation (one AssayRecord = one assay invocation)
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any


# Schema version. Bump on incompatible AssayRecord field changes.
ASSAY_SCHEMA_VERSION = "1.0"


class AssayKind(Enum):
    """Categorical type of measurement.

    Drives downstream interpretation: the fitter dispatches on
    ``AssayKind`` to decide which CalibrationEntry to emit. Add new
    kinds by extending this enum + the fit-function dispatch tables in
    ``data/validation/`` (Node 20).
    """

    DROPLET_SIZE_DISTRIBUTION = "droplet_size_distribution"
    """Laser diffraction or microscopy DSD; informs L1 kernel fits."""

    INTERFACIAL_TENSION = "interfacial_tension"
    """Pendant drop or Wilhelmy plate; informs IFT model."""

    DISPERSED_VISCOSITY = "dispersed_viscosity"
    """Rheometer flow curve; informs L1 mu_d input + Cross model."""

    PORE_SIZE = "pore_size"
    """Cryo-SEM, tracer exclusion, SEC probe, porosimetry; informs L2 fit."""

    POROSITY = "porosity"
    """Wet/dry swelling or tracer exclusion; informs L2 porosity."""

    GELATION_ONSET = "gelation_onset"
    """Rheology G'-G'' crossover; informs L2 timing constants."""

    CROSSLINK_CONVERSION = "crosslink_conversion"
    """Ninhydrin/TNBS residual amine assay; informs L3 kinetics."""

    BULK_MODULUS = "bulk_modulus"
    """Compression / oscillatory rheometry; informs L4 G_DN fit."""

    STATIC_BINDING_CAPACITY = "static_binding_capacity"
    """Isotherm point or full curve; informs M3 q_max."""

    DYNAMIC_BINDING_CAPACITY = "dynamic_binding_capacity"
    """Breakthrough curve at given residence time; informs M3 DBC10/DBC50."""

    LIGAND_DENSITY = "ligand_density"
    """Colorimetric or elemental assay; informs M2 functional density."""

    ACTIVITY_RETENTION = "activity_retention"
    """Binding-assay vs known concentration; informs M2 activity_retention."""


@dataclass
class Replicate:
    """One numerical replicate of a measurement.

    Replicates carry their own (optional) standard deviation when the
    assay reports per-replicate uncertainty (e.g. instrument ± from a
    laser-diffraction histogram width). Most assays leave ``std`` at
    zero and let the calibration fitter compute scatter across the
    Replicate list.
    """
    value: float
    std: float = 0.0
    flag: str = ""        # "" | "outlier" | "censored_low" | "censored_high"
    notes: str = ""


@dataclass
class AssayRecord:
    """One wet-lab assay invocation, ready for fitter ingest.

    Attributes:
        record_id: Unique identifier (lab-notebook entry, LIMS ID, etc.).
        kind: Categorical assay type (drives fitter dispatch).
        units: Physical unit string for ``replicates[*].value``.
            Use SI where possible (m, K, Pa, mol/m^3); document deviations
            in ``notes`` and let the fitter coerce. Same convention as
            ``CalibrationEntry.units``.
        replicates: Per-replicate measurements.
        process_conditions: Free-form dict of conditions held constant
            across replicates (RPM, temperature, surfactant_conc, pH,
            salt_conc, etc.). Keys SHOULD be SI-named; the fitter uses
            this dict to build the calibration-domain matrix.
        sample_id: Material lot identifier — what was actually measured.
        instrument: Instrument model + serial / settings string.
        operator: Operator initials (auditability).
        notebook_ref: Lab-notebook page or LIMS link.
        timestamp_utc: ISO-8601 UTC stamp of measurement.
        target_module: Optional hint for the fitter ("L1"/"L2"/"L3"/"L4"/"M2"/"M3"
            or ""). When set, the fitter only emits CalibrationEntries
            against this module. When empty, the kind is used to infer.
        notes: Free-form notes (anomalies, calibration cycles, etc.).

    The dataclass is JSON-serialisable via ``to_dict`` / ``from_dict``;
    same round-trip discipline as CalibrationEntry. Carry one record per
    assay invocation; aggregate across runs by storing many records in
    the same JSON file (`data/validation/<kind>/*.json`).
    """
    record_id: str
    kind: AssayKind
    units: str
    replicates: list[Replicate] = field(default_factory=list)
    process_conditions: dict[str, Any] = field(default_factory=dict)
    sample_id: str = ""
    instrument: str = ""
    operator: str = ""
    notebook_ref: str = ""
    timestamp_utc: str = ""
    target_module: str = ""
    notes: str = ""

    def __post_init__(self):
        if not self.timestamp_utc:
            self.timestamp_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")

    # ── Statistics on the replicate set ─────────────────────────────────

    def values(self) -> list[float]:
        return [r.value for r in self.replicates if r.flag != "outlier"]

    def mean(self) -> float:
        v = self.values()
        if not v:
            return float("nan")
        return sum(v) / len(v)

    def std(self) -> float:
        """Sample standard deviation across non-outlier replicates."""
        v = self.values()
        if len(v) < 2:
            return 0.0
        m = self.mean()
        return (sum((x - m) ** 2 for x in v) / (len(v) - 1)) ** 0.5

    def cv(self) -> float:
        """Coefficient of variation (std/|mean|); NaN if mean is zero."""
        m = self.mean()
        if abs(m) < 1e-300:
            return float("nan")
        return self.std() / abs(m)

    def n_replicates(self) -> int:
        return len(self.values())

    # ── JSON round-trip ─────────────────────────────────────────────────

    def to_dict(self) -> dict:
        d = asdict(self)
        d["kind"] = self.kind.value           # enum -> string
        d["schema_version"] = ASSAY_SCHEMA_VERSION
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "AssayRecord":
        # Tolerate older schemas by ignoring unknown keys; bump
        # ASSAY_SCHEMA_VERSION when an incompatible change requires
        # a real migration.
        d = dict(d)
        version = d.pop("schema_version", "1.0")
        if version != ASSAY_SCHEMA_VERSION:
            # v7.0 has only one schema version. v7.1+ should add
            # field translation here when bumping.
            pass

        replicates = [
            Replicate(**r) if isinstance(r, dict) else r
            for r in d.pop("replicates", [])
        ]
        kind_raw = d.pop("kind")
        kind = AssayKind(kind_raw) if isinstance(kind_raw, str) else kind_raw
        return cls(replicates=replicates, kind=kind, **d)
