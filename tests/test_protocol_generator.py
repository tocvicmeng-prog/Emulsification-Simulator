"""Tests for emulsim.protocols.protocol_generator.

Covers:
  - Protocol generation for standard crosslinkers (genipin, edc_nhs, pegda_uv)
  - Protocol generation for reagent_profiles (ech_activation)
  - Reagent mass calculation correctness
  - User parameter propagation into ProtocolDocument
  - ValueError on unknown key
  - Broad coverage: every CROSSLINKER and REAGENT_PROFILE key generates without error
"""

from __future__ import annotations

import pytest

from emulsim.protocols.protocol_generator import ProtocolGenerator
from emulsim.protocols.protocol_document import ProtocolDocument
from emulsim.reagent_library import CROSSLINKERS
from emulsim.module2_functionalization.reagent_profiles import REAGENT_PROFILES


# ═══════════════════════════════════════════════════════════════════════════
#  TestProtocolGeneratorCrosslinkers
# ═══════════════════════════════════════════════════════════════════════════


class TestProtocolGeneratorCrosslinkers:
    """Tests covering generate() with source='crosslinkers'."""

    def test_generate_crosslinker_genipin_title_not_empty(self):
        """Genipin protocol must have a non-empty title."""
        doc = ProtocolGenerator().generate(
            "genipin",
            source="crosslinkers",
            temperature_K=310.15,
            time_s=86400,
            concentration_mM=2.0,
            pH=7.4,
        )
        assert doc.title != ""

    def test_generate_crosslinker_genipin_reagent_table_min_length(self):
        """Genipin protocol reagent_table must have at least 2 entries."""
        doc = ProtocolGenerator().generate(
            "genipin",
            source="crosslinkers",
            temperature_K=310.15,
            time_s=86400,
            concentration_mM=2.0,
            pH=7.4,
        )
        assert len(doc.reagent_table) >= 2

    def test_generate_crosslinker_genipin_procedure_steps_min_length(self):
        """Genipin protocol procedure_steps must have at least 5 steps."""
        doc = ProtocolGenerator().generate(
            "genipin",
            source="crosslinkers",
            temperature_K=310.15,
            time_s=86400,
            concentration_mM=2.0,
            pH=7.4,
        )
        assert len(doc.procedure_steps) >= 5

    def test_generate_crosslinker_genipin_safety_warnings_min_length(self):
        """Genipin protocol must include at least 1 safety warning."""
        doc = ProtocolGenerator().generate(
            "genipin",
            source="crosslinkers",
            temperature_K=310.15,
            time_s=86400,
            concentration_mM=2.0,
            pH=7.4,
        )
        assert len(doc.safety_warnings) >= 1

    def test_generate_crosslinker_genipin_qc_targets_min_length(self):
        """Genipin protocol must include at least 1 QC target."""
        doc = ProtocolGenerator().generate(
            "genipin",
            source="crosslinkers",
            temperature_K=310.15,
            time_s=86400,
            concentration_mM=2.0,
            pH=7.4,
        )
        assert len(doc.qc_targets) >= 1

    def test_generate_uv_crosslinker_has_photoinitiator(self):
        """PEGDA/UV protocol reagent_table must mention 'photoinitiator' or 'Irgacure'."""
        doc = ProtocolGenerator().generate(
            "pegda_uv",
            source="crosslinkers",
            temperature_K=298.15,
            time_s=300,
            concentration_mM=10.0,
            pH=7.4,
        )
        names_and_roles = " ".join(
            r.name + " " + r.role for r in doc.reagent_table
        )
        assert "photoinitiator" in names_and_roles.lower() or "irgacure" in names_and_roles.lower()

    def test_generate_edc_nhs_has_nhs_coreagent(self):
        """EDC/NHS protocol reagent_table must include an NHS entry."""
        doc = ProtocolGenerator().generate(
            "edc_nhs",
            source="crosslinkers",
            temperature_K=298.15,
            time_s=3600,
            concentration_mM=5.0,
            pH=7.4,
        )
        reagent_names = [r.name for r in doc.reagent_table]
        assert any("NHS" in name for name in reagent_names), (
            f"No NHS entry found in reagent_table: {reagent_names}"
        )

    def test_generate_invalid_key_raises_value_error(self):
        """generate() with an unknown key must raise ValueError."""
        with pytest.raises(ValueError):
            ProtocolGenerator().generate(
                "bad_key_xyz",
                source="crosslinkers",
                temperature_K=298.15,
                time_s=3600,
                concentration_mM=1.0,
                pH=7.4,
            )


# ═══════════════════════════════════════════════════════════════════════════
#  TestProtocolGeneratorReagentProfiles
# ═══════════════════════════════════════════════════════════════════════════


class TestProtocolGeneratorReagentProfiles:
    """Tests covering generate() with source='reagent_profiles'."""

    def test_generate_reagent_profile_ech_safety_warnings_carcinogen_and_alkaline(self):
        """ECH activation at pH 12 must produce at least 2 safety warnings."""
        doc = ProtocolGenerator().generate(
            "ech_activation",
            source="reagent_profiles",
            temperature_K=298.15,
            time_s=7200,
            concentration_mM=10.0,
            pH=12.0,
        )
        assert len(doc.safety_warnings) >= 2


# ═══════════════════════════════════════════════════════════════════════════
#  TestReagentMassCalculation
# ═══════════════════════════════════════════════════════════════════════════


class TestReagentMassCalculation:
    """Tests covering mass calculation correctness."""

    def test_reagent_mass_genipin(self):
        """Genipin at 2.0 mM, 1.0 mL bead volume: mass = 2.0 * 226.23 * 10 / 1000 = 4.5246 mg.

        The amount string in the reagent_table entry must contain '4.5'.
        """
        doc = ProtocolGenerator().generate(
            "genipin",
            source="crosslinkers",
            temperature_K=310.15,
            time_s=86400,
            concentration_mM=2.0,
            pH=7.4,
            bead_volume_mL=1.0,
        )
        genipin_entry = next(
            (r for r in doc.reagent_table if "enipin" in r.name and "secondary" not in r.name.lower()),
            None,
        )
        assert genipin_entry is not None, "No genipin entry found in reagent_table"
        assert "4.5" in genipin_entry.amount, (
            f"Expected '4.5' in amount string, got: {genipin_entry.amount!r}"
        )


# ═══════════════════════════════════════════════════════════════════════════
#  TestUserParamsPropagation
# ═══════════════════════════════════════════════════════════════════════════


class TestUserParamsPropagation:
    """Tests that user-supplied parameters are stored faithfully in the document."""

    def _make_doc(self) -> ProtocolDocument:
        return ProtocolGenerator().generate(
            "genipin",
            source="crosslinkers",
            temperature_K=310.15,
            time_s=86400,
            concentration_mM=2.0,
            pH=7.4,
            bead_volume_mL=1.0,
        )

    def test_user_temperature_propagates(self):
        assert self._make_doc().user_temperature_K == 310.15

    def test_user_time_propagates(self):
        assert self._make_doc().user_time_s == 86400

    def test_user_concentration_propagates(self):
        assert self._make_doc().user_concentration_mM == 2.0

    def test_user_ph_propagates(self):
        assert self._make_doc().user_pH == 7.4


# ═══════════════════════════════════════════════════════════════════════════
#  TestBroadCoverage
# ═══════════════════════════════════════════════════════════════════════════


class TestBroadCoverage:
    """Smoke-test: generate() succeeds for every registered key."""

    def test_all_crosslinkers_generate(self):
        """Every key in CROSSLINKERS must generate a ProtocolDocument without error."""
        gen = ProtocolGenerator()
        for key in CROSSLINKERS:
            doc = gen.generate(
                key,
                source="crosslinkers",
                temperature_K=298.15,
                time_s=3600,
                concentration_mM=1.0,
                pH=7.4,
            )
            assert isinstance(doc, ProtocolDocument), (
                f"generate('{key}', source='crosslinkers') did not return ProtocolDocument"
            )

    def test_all_reagent_profiles_generate(self):
        """Every key in REAGENT_PROFILES must generate a ProtocolDocument without error."""
        gen = ProtocolGenerator()
        for key in REAGENT_PROFILES:
            doc = gen.generate(
                key,
                source="reagent_profiles",
                temperature_K=298.15,
                time_s=3600,
                concentration_mM=1.0,
                pH=7.4,
            )
            assert isinstance(doc, ProtocolDocument), (
                f"generate('{key}', source='reagent_profiles') did not return ProtocolDocument"
            )
