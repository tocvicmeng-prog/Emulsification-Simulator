"""Tests for emulsim.protocols.protocol_document.

Covers:
  - to_markdown() section headings
  - Temperature K -> Celsius conversion in Markdown output
  - Time formatting: hours and minutes
  - Safety warning prefix normalisation
  - Reagent table row contents
  - Pipe character escaping in Markdown table cells
  - Empty document (minimal data) does not raise
"""

from __future__ import annotations

import pytest

from emulsim.protocols.protocol_document import (
    ProtocolDocument,
    ProtocolStep,
    ReagentRequirement,
)


# ═══════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════


def _minimal_doc(**overrides) -> ProtocolDocument:
    """Return a ProtocolDocument with sensible defaults for testing."""
    defaults = dict(
        title="Test Protocol",
        reagent_key="test_reagent",
        generation_timestamp="2026-04-16T12:00:00",
        user_temperature_K=298.15,
        user_time_s=3600.0,
        user_concentration_mM=1.0,
        user_pH=7.4,
        user_bead_volume_mL=1.0,
        reagent_table=[],
        procedure_steps=[],
        qc_targets=[],
        safety_warnings=[],
        mechanism_summary="Test mechanism.",
    )
    defaults.update(overrides)
    return ProtocolDocument(**defaults)


# ═══════════════════════════════════════════════════════════════════════════
#  TestMarkdownSections
# ═══════════════════════════════════════════════════════════════════════════


class TestMarkdownSections:
    """Tests that all expected section headings appear in to_markdown() output."""

    def _md(self) -> str:
        return _minimal_doc().to_markdown()

    def test_contains_protocol_title_heading(self):
        """Markdown must start with '# Protocol:' heading."""
        assert "# Protocol:" in self._md()

    def test_contains_reagents_required_heading(self):
        assert "## Reagents Required" in self._md()

    def test_contains_procedure_heading(self):
        assert "## Procedure" in self._md()

    def test_contains_qc_targets_heading(self):
        assert "## QC Targets" in self._md()

    def test_contains_safety_warnings_heading(self):
        assert "## Safety Warnings" in self._md()


# ═══════════════════════════════════════════════════════════════════════════
#  TestTemperatureConversion
# ═══════════════════════════════════════════════════════════════════════════


class TestTemperatureConversion:
    """Tests Kelvin -> Celsius conversion displayed in Markdown."""

    def test_temperature_37_celsius(self):
        """310.15 K must appear as '37.0' (Celsius) in the Markdown output."""
        md = _minimal_doc(user_temperature_K=310.15).to_markdown()
        assert "37.0" in md


# ═══════════════════════════════════════════════════════════════════════════
#  TestTimeFormatting
# ═══════════════════════════════════════════════════════════════════════════


class TestTimeFormatting:
    """Tests human-readable time formatting in Markdown."""

    def test_time_formatting_hours(self):
        """86400 s must appear as '24.0 h' in the Markdown output."""
        md = _minimal_doc(user_time_s=86400).to_markdown()
        assert "24.0 h" in md

    def test_time_formatting_minutes(self):
        """1800 s must appear as '30 min' in the Markdown output."""
        md = _minimal_doc(user_time_s=1800).to_markdown()
        assert "30 min" in md


# ═══════════════════════════════════════════════════════════════════════════
#  TestSafetyWarnings
# ═══════════════════════════════════════════════════════════════════════════


class TestSafetyWarnings:
    """Tests safety warning formatting in Markdown."""

    def test_safety_warning_prefix_added(self):
        """Each safety warning in the Markdown must be prefixed with 'WARNING:'."""
        md = _minimal_doc(safety_warnings=["Toxic reagent"]).to_markdown()
        assert "WARNING: Toxic reagent" in md

    def test_safety_warning_prefix_not_doubled(self):
        """Warnings that already carry 'WARNING:' prefix must not be doubled."""
        md = _minimal_doc(safety_warnings=["WARNING: Pre-prefixed"]).to_markdown()
        assert "WARNING: WARNING:" not in md


# ═══════════════════════════════════════════════════════════════════════════
#  TestReagentTable
# ═══════════════════════════════════════════════════════════════════════════


class TestReagentTable:
    """Tests that reagent table rows are rendered with name and CAS."""

    def _doc_with_reagent(self) -> ProtocolDocument:
        reagent = ReagentRequirement(
            name="Genipin",
            cas="6902-77-8",
            amount="4.5 mg",
            grade=">=98%",
            role="crosslinker",
        )
        return _minimal_doc(reagent_table=[reagent])

    def test_reagent_name_in_markdown(self):
        """Reagent name must appear in the Markdown output."""
        md = self._doc_with_reagent().to_markdown()
        assert "Genipin" in md

    def test_reagent_cas_in_markdown(self):
        """Reagent CAS number must appear in the Markdown output."""
        md = self._doc_with_reagent().to_markdown()
        assert "6902-77-8" in md


# ═══════════════════════════════════════════════════════════════════════════
#  TestPipeEscaping
# ═══════════════════════════════════════════════════════════════════════════


class TestPipeEscaping:
    """Tests that pipe characters inside reagent names are escaped for Markdown tables."""

    def test_pipe_in_name_is_escaped(self):
        """A reagent name containing '|' must have it escaped as '\\|' in Markdown."""
        reagent = ReagentRequirement(
            name="Reagent A | Variant B",
            cas="N/A",
            amount="1.0 mg",
            grade="AR grade",
            role="test",
        )
        md = _minimal_doc(reagent_table=[reagent]).to_markdown()
        assert r"\|" in md


# ═══════════════════════════════════════════════════════════════════════════
#  TestEmptyDocument
# ═══════════════════════════════════════════════════════════════════════════


class TestEmptyDocument:
    """Tests robustness with minimal / empty data."""

    def test_empty_document_does_not_raise(self):
        """to_markdown() on a minimal ProtocolDocument must not raise any exception."""
        doc = _minimal_doc()
        md = doc.to_markdown()
        assert isinstance(md, str)
        assert len(md) > 0
