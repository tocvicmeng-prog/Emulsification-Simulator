"""Calibration store: load, save, query, and apply calibration data.

v6.0-alpha: In-memory store backed by JSON file. Applies measured values
to FunctionalMediaContract, logging every override for transparency
(audit F2 requirement).
"""

from __future__ import annotations

import copy
import json
import logging
from pathlib import Path
from typing import Any

from .calibration_data import CalibrationEntry

logger = logging.getLogger(__name__)


class CalibrationStore:
    """Load, save, query, and apply user-supplied calibration data.

    Usage::

        store = CalibrationStore()
        store.load_json("my_calibration.json")
        fmc_calibrated, overrides = store.apply_to_fmc(fmc)
        # overrides is a list of human-readable override descriptions
    """

    def __init__(self) -> None:
        self._entries: list[CalibrationEntry] = []

    def add(self, entry: CalibrationEntry) -> None:
        """Add a calibration entry."""
        self._entries.append(entry)

    def load_json(self, path: str | Path) -> int:
        """Load calibration entries from a JSON file.

        JSON format: list of dicts, each matching CalibrationEntry fields.

        Args:
            path: Path to JSON file.

        Returns:
            Number of entries loaded.
        """
        path = Path(path)
        with open(path, encoding="utf-8") as f:
            data = json.load(f)

        if not isinstance(data, list):
            raise ValueError(f"Expected JSON array, got {type(data).__name__}")

        count = 0
        for item in data:
            try:
                entry = CalibrationEntry.from_dict(item)
                self._entries.append(entry)
                count += 1
            except (KeyError, TypeError, ValueError) as e:
                logger.warning("Skipping invalid calibration entry: %s", e)

        logger.info("Loaded %d calibration entries from %s", count, path)
        return count

    def save_json(self, path: str | Path) -> None:
        """Save all entries to a JSON file."""
        path = Path(path)
        data = [e.to_dict() for e in self._entries]
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        logger.info("Saved %d calibration entries to %s", len(data), path)

    def query(
        self,
        profile_key: str,
        parameter_name: str = "",
    ) -> list[CalibrationEntry]:
        """Query entries by profile key and optionally parameter name.

        Args:
            profile_key: Reagent profile key to match.
            parameter_name: Optional parameter name filter.

        Returns:
            List of matching CalibrationEntry objects.
        """
        results = [e for e in self._entries if e.profile_key == profile_key]
        if parameter_name:
            results = [e for e in results if e.parameter_name == parameter_name]
        return results

    @property
    def entries(self) -> list[CalibrationEntry]:
        """All stored entries (read-only view)."""
        return list(self._entries)

    def __len__(self) -> int:
        return len(self._entries)

    def apply_to_fmc(self, fmc: Any) -> tuple[Any, list[str]]:
        """Apply calibration overrides to a FunctionalMediaContract.

        For each calibration entry that matches the FMC's profile context,
        overrides the corresponding FMC field. Returns a modified copy and
        a log of all overrides applied.

        Supported parameter_name -> FMC field mappings:
            "estimated_q_max" -> fmc.estimated_q_max
            "q_max" -> fmc.estimated_q_max (alias)
            "K_L" -> stored in calibration_overrides
            "activity_retention" -> stored in calibration_overrides

        Args:
            fmc: FunctionalMediaContract instance.

        Returns:
            Tuple of (modified_fmc_copy, list_of_override_descriptions).
        """
        fmc_out = copy.deepcopy(fmc)
        overrides: list[str] = []

        # Map parameter names to FMC fields
        _FMC_FIELD_MAP = {
            "estimated_q_max": "estimated_q_max",
            "q_max": "estimated_q_max",
            "charge_density": "charge_density",
            "functional_ligand_density": "functional_ligand_density",
        }

        for entry in self._entries:
            field_name = _FMC_FIELD_MAP.get(entry.parameter_name)
            if field_name and hasattr(fmc_out, field_name):
                old_val = getattr(fmc_out, field_name)
                setattr(fmc_out, field_name, entry.measured_value)
                desc = (
                    f"OVERRIDE: {field_name} {old_val:.4g} -> {entry.measured_value:.4g} "
                    f"({entry.units}, {entry.confidence}, {entry.source_reference})"
                )
                overrides.append(desc)
                logger.info(desc)

                # Update confidence to reflect calibration
                if hasattr(fmc_out, 'q_max_confidence') and "q_max" in entry.parameter_name:
                    fmc_out.q_max_confidence = "calibrated"
                if hasattr(fmc_out, 'confidence_tier'):
                    fmc_out.confidence_tier = "calibrated"

        if overrides:
            logger.info("Applied %d calibration overrides to FMC", len(overrides))
        else:
            logger.debug("No calibration entries matched FMC context")

        return fmc_out, overrides
