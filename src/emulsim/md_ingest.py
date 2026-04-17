"""Node F5 Phase 1 (v8.3-alpha): ingest MARTINI coarse-grained MD
records into EmulSim MaterialProperties.

EmulSim itself does **not** run MARTINI simulations — we accept
user-supplied output as a structured JSON record. This matches the
Node 32 roadmap's "ingest-only default" scoping: MD is the
authoritative source of χ parameters, and we provide the plumbing to
move those numbers into a continuum Flory-Huggins model (currently
used by the cellulose NIPS L2 solver).

See :file:`docs/f5_md_ingest_protocol.md` for the scientific basis,
schema definition, and API examples.

JSON schema (Phase 1)
---------------------
Required top-level keys:

- ``source``            — str, e.g. ``"martini_3.0"``
- ``system_description`` — str, free-form
- ``beads``             — dict mapping EmulSim role → MARTINI bead
                          type, e.g. ``{"polymer": "P3", ...}``
- ``chi``               — dict of effective χ parameters, keys from
                          ``{"polymer_solvent", "polymer_nonsolvent",
                          "solvent_nonsolvent"}``; all three optional
                          (missing = "not measured")
- ``diagnostics``       — dict, free-form; conventionally includes
                          ``temperature_K``, ``simulation_ns``, etc.

Optional:

- ``paper_doi``         — str
- ``notes``             — str
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


# Required top-level keys in a MARTINI record JSON.
_REQUIRED_KEYS = frozenset({
    "source", "system_description", "beads", "chi", "diagnostics",
})

# Recognised χ sub-keys (extras are preserved but not applied).
_KNOWN_CHI_KEYS = frozenset({
    "polymer_solvent", "polymer_nonsolvent", "solvent_nonsolvent",
})


@dataclass
class MartiniRecord:
    """Parsed MARTINI MD record.

    Attributes
    ----------
    source : str
        Identifier of the MD force field version / package.
    system_description : str
        Free-form description of what was simulated.
    beads : dict
        Mapping EmulSim role → MARTINI bead type.
    chi : dict
        Flory-Huggins χ parameters. Keys from
        :data:`_KNOWN_CHI_KEYS`; any missing key means "not measured"
        and should be left alone downstream.
    diagnostics : dict
        Free-form dictionary of simulation metadata
        (temperature, duration, chain length, etc.).
    paper_doi : str
        Optional DOI of the accompanying publication.
    notes : str
        Free-form notes.
    extra : dict
        Any additional keys present in the JSON that aren't part of
        the schema — preserved for forward compatibility.
    """

    source: str
    system_description: str
    beads: dict[str, str]
    chi: dict[str, float]
    diagnostics: dict[str, Any]
    paper_doi: str = ""
    notes: str = ""
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialise back to the JSON shape (round-trip)."""
        payload = {
            "source": self.source,
            "system_description": self.system_description,
            "beads": dict(self.beads),
            "chi": dict(self.chi),
            "diagnostics": dict(self.diagnostics),
        }
        if self.paper_doi:
            payload["paper_doi"] = self.paper_doi
        if self.notes:
            payload["notes"] = self.notes
        payload.update(self.extra)
        return payload


def load_martini_record(path: Path | str) -> MartiniRecord:
    """Parse a MARTINI MD record JSON file into a :class:`MartiniRecord`.

    Raises
    ------
    ValueError
        If any required top-level key is missing, or if χ contains a
        non-numeric or non-finite value.
    FileNotFoundError
        If ``path`` does not exist.
    json.JSONDecodeError
        If the file is not valid JSON.
    """
    p = Path(path)
    with open(p, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return _from_dict(raw)


def _from_dict(raw: dict[str, Any]) -> MartiniRecord:
    """Shared parser for ``load_martini_record`` and round-trip tests."""
    missing = _REQUIRED_KEYS - set(raw.keys())
    if missing:
        raise ValueError(
            f"MARTINI record missing required keys: {sorted(missing)}. "
            f"Required: {sorted(_REQUIRED_KEYS)}"
        )

    chi_raw = raw["chi"]
    if not isinstance(chi_raw, dict):
        raise ValueError(
            f"`chi` must be a dict, got {type(chi_raw).__name__}"
        )
    chi_validated: dict[str, float] = {}
    for k, v in chi_raw.items():
        if k in _KNOWN_CHI_KEYS:
            try:
                fv = float(v)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"chi[{k!r}] is not numeric: {v!r}"
                ) from exc
            if not math.isfinite(fv):
                raise ValueError(
                    f"chi[{k!r}] is not finite: {fv}"
                )
            chi_validated[k] = fv
        else:
            # Unknown chi key — keep it in extra-chi for future use
            chi_validated[k] = v

    extra = {
        k: v for k, v in raw.items()
        if k not in _REQUIRED_KEYS
        and k not in ("paper_doi", "notes")
    }

    return MartiniRecord(
        source=str(raw["source"]),
        system_description=str(raw["system_description"]),
        beads=dict(raw["beads"]),
        chi=chi_validated,
        diagnostics=dict(raw["diagnostics"]),
        paper_doi=str(raw.get("paper_doi", "")),
        notes=str(raw.get("notes", "")),
        extra=extra,
    )


def apply_chi_to_props(props, record: MartiniRecord) -> list[str]:
    """Apply the record's χ parameters to a ``MaterialProperties``.

    Current mapping (Phase 1, cellulose NIPS only):

    - ``polymer_solvent``    → ``props.chi_PS_cellulose``
    - ``polymer_nonsolvent`` → ``props.chi_PN_cellulose``
    - ``solvent_nonsolvent`` → ``props.chi_SN_cellulose``

    Missing χ sub-keys leave the corresponding field untouched.
    Non-cellulose fields on ``MaterialProperties`` are never
    modified.

    Returns
    -------
    list[str]
        Names of the MaterialProperties fields that were written,
        in application order.

    Raises
    ------
    ValueError
        If any applied χ is non-finite. (``_from_dict`` already
        validates this at load time; this guard catches records
        constructed in-memory with corrupted values.)
    """
    mapping = {
        "polymer_solvent": "chi_PS_cellulose",
        "polymer_nonsolvent": "chi_PN_cellulose",
        "solvent_nonsolvent": "chi_SN_cellulose",
    }
    written: list[str] = []
    for chi_key, field_name in mapping.items():
        if chi_key not in record.chi:
            continue
        value = record.chi[chi_key]
        fv = float(value)
        if not math.isfinite(fv):
            raise ValueError(
                f"chi[{chi_key!r}] is non-finite ({fv}); refusing to apply"
            )
        setattr(props, field_name, fv)
        written.append(field_name)
    return written


def save_martini_record(record: MartiniRecord, path: Path | str) -> None:
    """Serialise a :class:`MartiniRecord` to JSON (round-trip helper)."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(record.to_dict(), f, indent=2)
