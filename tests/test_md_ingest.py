"""Tests for Node F5 Phase 1: MARTINI MD parameter ingest.

10 tests per docs/f5_md_ingest_protocol.md §6.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from emulsim.datatypes import MaterialProperties
from emulsim.md_ingest import (
    MartiniRecord,
    apply_chi_to_props,
    load_martini_record,
    save_martini_record,
)

_FIXTURE = (
    Path(__file__).resolve().parent.parent
    / "data" / "validation" / "md" / "example_martini_cellulose.json"
)


def _minimal_dict() -> dict:
    return {
        "source": "martini_3.0",
        "system_description": "test system",
        "beads": {"polymer": "P3", "solvent": "P4", "nonsolvent": "P4"},
        "chi": {
            "polymer_solvent": 0.4,
            "polymer_nonsolvent": 0.9,
            "solvent_nonsolvent": 0.3,
        },
        "diagnostics": {"temperature_K": 298.15},
    }


class TestLoadValid:
    def test_loads_reference_fixture(self):
        record = load_martini_record(_FIXTURE)
        assert record.source == "martini_3.0"
        assert "NaOH-urea" in record.system_description
        assert record.chi["polymer_solvent"] == pytest.approx(0.41)
        assert record.chi["polymer_nonsolvent"] == pytest.approx(0.88)
        assert record.chi["solvent_nonsolvent"] == pytest.approx(0.28)
        assert record.diagnostics["temperature_K"] == pytest.approx(298.15)
        assert record.paper_doi.startswith("10.")


class TestMissingKeys:
    def test_missing_required_key_raises(self, tmp_path):
        bad = _minimal_dict()
        del bad["chi"]
        p = tmp_path / "bad.json"
        p.write_text(json.dumps(bad))
        with pytest.raises(ValueError, match="missing required keys"):
            load_martini_record(p)

    def test_missing_source_raises(self, tmp_path):
        bad = _minimal_dict()
        del bad["source"]
        p = tmp_path / "bad.json"
        p.write_text(json.dumps(bad))
        with pytest.raises(ValueError):
            load_martini_record(p)


class TestExtraKeysPreserved:
    def test_forward_compatible_keys_go_into_extra(self, tmp_path):
        payload = _minimal_dict()
        payload["future_schema_v2_field"] = {"foo": "bar"}
        p = tmp_path / "extra.json"
        p.write_text(json.dumps(payload))
        record = load_martini_record(p)
        assert "future_schema_v2_field" in record.extra
        assert record.extra["future_schema_v2_field"] == {"foo": "bar"}


class TestPartialChi:
    def test_missing_chi_subkey_allowed(self, tmp_path):
        """Not every MD run measures all three χ. Missing sub-keys
        should parse OK; downstream apply leaves those fields alone.
        """
        payload = _minimal_dict()
        payload["chi"] = {"polymer_solvent": 0.4}   # only one sub-key
        p = tmp_path / "partial.json"
        p.write_text(json.dumps(payload))
        record = load_martini_record(p)
        assert "polymer_solvent" in record.chi
        assert "polymer_nonsolvent" not in record.chi


class TestApplyToProps:
    def test_applies_all_three_chi_fields(self):
        props = MaterialProperties()
        record = load_martini_record(_FIXTURE)
        written = apply_chi_to_props(props, record)
        assert set(written) == {
            "chi_PS_cellulose", "chi_PN_cellulose", "chi_SN_cellulose",
        }
        assert props.chi_PS_cellulose == pytest.approx(0.41)
        assert props.chi_PN_cellulose == pytest.approx(0.88)
        assert props.chi_SN_cellulose == pytest.approx(0.28)

    def test_does_not_modify_non_cellulose_fields(self):
        """Apply should not touch alginate / PLGA / agarose-chitosan
        fields on MaterialProperties.
        """
        props = MaterialProperties()
        pristine_K_alg = props.K_alg_modulus
        pristine_D_DCM = props.D_DCM_plga
        pristine_f_G = props.f_guluronate

        record = load_martini_record(_FIXTURE)
        apply_chi_to_props(props, record)

        assert props.K_alg_modulus == pristine_K_alg
        assert props.D_DCM_plga == pristine_D_DCM
        assert props.f_guluronate == pristine_f_G

    def test_non_finite_chi_raises_on_apply(self):
        # Construct in-memory (bypass load-time validation)
        record = MartiniRecord(
            source="x", system_description="x",
            beads={}, diagnostics={},
            chi={"polymer_solvent": float("inf")},
        )
        props = MaterialProperties()
        with pytest.raises(ValueError):
            apply_chi_to_props(props, record)


class TestChiValueValidation:
    def test_nan_chi_at_load_raises(self, tmp_path):
        payload = _minimal_dict()
        payload["chi"]["polymer_solvent"] = float("nan")
        p = tmp_path / "nan.json"
        # json.dumps turns NaN into "NaN" literal; write directly.
        p.write_text(
            '{"source":"x","system_description":"x","beads":{},'
            '"chi":{"polymer_solvent":NaN},'
            '"diagnostics":{}}'
        )
        with pytest.raises(ValueError, match="not finite"):
            load_martini_record(p)

    def test_negative_chi_allowed(self, tmp_path):
        """χ < 0 is physically valid (attractive mixing; better than
        athermal). Should parse and apply without error.
        """
        payload = _minimal_dict()
        payload["chi"]["polymer_solvent"] = -0.1
        p = tmp_path / "neg.json"
        p.write_text(json.dumps(payload))
        record = load_martini_record(p)
        assert record.chi["polymer_solvent"] == pytest.approx(-0.1)
        props = MaterialProperties()
        apply_chi_to_props(props, record)
        assert props.chi_PS_cellulose == pytest.approx(-0.1)


class TestRoundTrip:
    def test_save_then_load_is_identity(self, tmp_path):
        record = load_martini_record(_FIXTURE)
        out = tmp_path / "roundtrip.json"
        save_martini_record(record, out)
        loaded = load_martini_record(out)
        assert loaded.source == record.source
        assert loaded.system_description == record.system_description
        assert loaded.beads == record.beads
        assert loaded.chi == record.chi
        assert loaded.diagnostics == record.diagnostics
        assert loaded.paper_doi == record.paper_doi
        assert loaded.notes == record.notes
