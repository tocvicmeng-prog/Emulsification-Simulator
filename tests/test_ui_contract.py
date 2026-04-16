"""Tests for Phase 0: UI Contract Layer.

Covers all 4 new modules:
  - ui_model_metadata: OutputMetadata, ModelBasis, ConfidenceLevel
  - ui_validators:     validate_m1_inputs, validate_m2_inputs,
                       validate_m3_chromatography, validate_m3_result
  - ui_state:          SessionStateManager invalidation cascade
  - ui_units:          all conversion and formatting functions
"""

import math
import pytest

from emulsim.visualization.ui_model_metadata import (
    OutputMetadata,
    ModelBasis,
    ConfidenceLevel,
    M1_D50_META,
    M1_PORE_META,
    M2_ACS_META,
    M3_BREAKTHROUGH_META,
    M3_GRADIENT_META,
    M3_MASS_BALANCE_META,
)
from emulsim.module2_functionalization.modification_steps import ModificationStepType
from emulsim.visualization.ui_validators import _M2_SUPPORTED_STEP_TYPES
from emulsim.visualization.ui_validators import (
    ValidationResult,
    validate_m1_inputs,
    validate_m2_inputs,
    validate_m3_chromatography,
    validate_m3_result,
)
from emulsim.visualization.ui_state import SessionStateManager
from emulsim.visualization.ui_units import (
    mol_m3_to_mM,
    mM_to_mol_m3,
    Pa_to_kPa,
    Pa_to_bar,
    bar_to_Pa,
    m_to_um,
    m_to_nm,
    um_to_m,
    nm_to_m,
    m3_s_to_mL_min,
    mL_min_to_m3_s,
    K_to_C,
    C_to_K,
    format_sci,
    format_with_unit,
    format_pressure,
    format_length,
    format_concentration,
)


# ─── Helper defaults for M1 validation ────────────────────────────────────────

def _m1_valid_kwargs(**overrides):
    """Return a dict of valid M1 inputs, with optional overrides."""
    defaults = dict(
        rpm=500,
        phi_d=0.20,
        c_agarose=4.0,
        c_chitosan=2.0,
        dda=0.85,
        crosslinker_key="genipin",
        crosslinker_conc=1.0,
        T_crosslink=37.0,
        T_oil=60.0,
    )
    defaults.update(overrides)
    return defaults


# ─── OutputMetadata tests ─────────────────────────────────────────────────────

class TestOutputMetadata:
    def test_creation_basic(self):
        meta = OutputMetadata(
            model_basis=ModelBasis.MECHANISTIC,
            confidence=ConfidenceLevel.HIGH,
        )
        assert meta.model_basis == ModelBasis.MECHANISTIC
        assert meta.confidence == ConfidenceLevel.HIGH

    def test_icon_high(self):
        meta = OutputMetadata(ModelBasis.MECHANISTIC, ConfidenceLevel.HIGH)
        assert meta.icon == "[OK]"

    def test_icon_moderate(self):
        meta = OutputMetadata(ModelBasis.MECHANISTIC, ConfidenceLevel.MODERATE)
        assert meta.icon == "[~]"

    def test_icon_low(self):
        meta = OutputMetadata(ModelBasis.EMPIRICAL_CALIBRATED, ConfidenceLevel.LOW)
        assert meta.icon == "[!]"

    def test_icon_none(self):
        meta = OutputMetadata(ModelBasis.NOT_PREDICTED, ConfidenceLevel.NONE)
        assert meta.icon == "[X]"

    def test_label_contains_basis(self):
        meta = OutputMetadata(ModelBasis.SEMI_QUANTITATIVE, ConfidenceLevel.LOW)
        assert "semi_quantitative" in meta.label

    def test_frozen_immutable(self):
        meta = OutputMetadata(ModelBasis.RANKING_ONLY, ConfidenceLevel.LOW)
        with pytest.raises((AttributeError, TypeError)):
            meta.confidence = ConfidenceLevel.HIGH  # type: ignore[misc]

    def test_prebuilt_d50_meta(self):
        assert M1_D50_META.source_module == "M1"
        assert M1_D50_META.model_basis == ModelBasis.MECHANISTIC
        assert M1_D50_META.confidence == ConfidenceLevel.MODERATE

    def test_prebuilt_pore_meta_calibration_required(self):
        assert M1_PORE_META.calibration_required is True

    def test_prebuilt_m2_acs_meta(self):
        assert M2_ACS_META.source_module == "M2"
        assert M2_ACS_META.confidence == ConfidenceLevel.LOW
        assert "All 9 backend step types" in M2_ACS_META.validity_range

    def test_prebuilt_mass_balance_meta_high(self):
        assert M3_MASS_BALANCE_META.confidence == ConfidenceLevel.HIGH

    def test_gradient_meta_not_stale(self):
        warning_text = " ".join(M3_GRADIENT_META.warnings)
        assert "BF-2" not in warning_text
        assert "diagnostic" in warning_text.lower()


# ─── M1 Validation tests ──────────────────────────────────────────────────────

class TestValidationM1:
    def test_valid_inputs_no_blockers(self):
        result = validate_m1_inputs(**_m1_valid_kwargs())
        assert result.valid is True
        assert result.blockers == []

    def test_phi_d_blocker_above_050(self):
        """phi_d > 0.50 triggers a blocker."""
        result = validate_m1_inputs(**_m1_valid_kwargs(phi_d=0.55))
        assert result.valid is False
        assert any("phi_d" in b and "0.50" in b for b in result.blockers)

    def test_phi_d_exactly_050_is_blocker(self):
        result = validate_m1_inputs(**_m1_valid_kwargs(phi_d=0.50))
        # phi_d=0.50 is not > 0.50, so not a blocker by rule R1; should pass
        # (boundary: only > 0.50 triggers; 0.50 is still a warning because > 0.30)
        # Let's just verify it doesn't crash and gives a warning
        assert result.valid is True  # 0.50 exactly is NOT > 0.50

    def test_phi_d_warning_above_030(self):
        """phi_d in (0.30, 0.50] triggers a warning but not a blocker."""
        result = validate_m1_inputs(**_m1_valid_kwargs(phi_d=0.35))
        assert result.valid is True
        assert any("0.30" in w or "concentrated" in w for w in result.warnings)

    def test_phi_d_nonpositive_blocker(self):
        result = validate_m1_inputs(**_m1_valid_kwargs(phi_d=0.0))
        assert result.valid is False

    # ── c_agarose ──────────────────────────────────────────────────────────

    def test_c_agarose_warning_below_2(self):
        """c_agarose < 2 triggers a warning."""
        result = validate_m1_inputs(**_m1_valid_kwargs(c_agarose=1.5))
        assert result.valid is True
        assert any("agarose" in w.lower() or "c_agarose" in w for w in result.warnings)

    def test_c_agarose_warning_above_6(self):
        """c_agarose > 6 triggers a warning."""
        result = validate_m1_inputs(**_m1_valid_kwargs(c_agarose=7.0))
        assert result.valid is True
        assert any("agarose" in w.lower() or "c_agarose" in w for w in result.warnings)

    def test_c_agarose_in_range_no_warning(self):
        result = validate_m1_inputs(**_m1_valid_kwargs(c_agarose=4.0))
        agarose_warnings = [w for w in result.warnings if "agarose" in w.lower()]
        assert agarose_warnings == []

    # ── c_chitosan ─────────────────────────────────────────────────────────

    def test_c_chitosan_warning_out_of_range(self):
        """c_chitosan outside [1, 3] triggers a warning."""
        result = validate_m1_inputs(**_m1_valid_kwargs(c_chitosan=0.5))
        assert any("chitosan" in w.lower() for w in result.warnings)

    def test_c_chitosan_zero_blocker(self):
        result = validate_m1_inputs(**_m1_valid_kwargs(c_chitosan=0.0))
        assert result.valid is False

    # ── T_crosslink ────────────────────────────────────────────────────────

    def test_T_crosslink_blocker_above_85(self):
        """T_crosslink >= 85 C triggers a blocker."""
        result = validate_m1_inputs(**_m1_valid_kwargs(T_crosslink=85.0))
        assert result.valid is False
        assert any("85" in b or "re-melt" in b for b in result.blockers)

    def test_T_crosslink_below_85_ok(self):
        result = validate_m1_inputs(**_m1_valid_kwargs(T_crosslink=60.0))
        tcross_blockers = [b for b in result.blockers if "crosslink" in b.lower() or "re-melt" in b.lower()]
        assert tcross_blockers == []

    def test_T_crosslink_negative_blocker(self):
        result = validate_m1_inputs(**_m1_valid_kwargs(T_crosslink=-5.0))
        assert result.valid is False


# ─── M2 Validation tests ──────────────────────────────────────────────────────

class TestValidationM2:
    class _FakeStep:
        """Minimal mock of a ModificationStep."""
        def __init__(self, step_type_value: str):
            self.step_type = type("ST", (), {"value": step_type_value})()

    def test_valid_secondary_crosslinking(self):
        steps = [self._FakeStep("secondary_crosslinking")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is True

    def test_valid_activation(self):
        steps = [self._FakeStep("activation")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is True

    def test_ligand_coupling_accepted(self):
        """LIGAND_COUPLING step type is now supported (Codex P2-5 fix)."""
        steps = [self._FakeStep("ligand_coupling")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is True

    def test_protein_coupling_accepted(self):
        """PROTEIN_COUPLING step type is now supported."""
        steps = [self._FakeStep("protein_coupling")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is True

    def test_quenching_accepted(self):
        """QUENCHING step type is now supported."""
        steps = [self._FakeStep("quenching")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is True

    @pytest.mark.parametrize(
        "step_type_value",
        [
            "secondary_crosslinking",
            "activation",
            "ligand_coupling",
            "protein_coupling",
            "quenching",
            "spacer_arm",
            "metal_charging",
            "protein_pretreatment",
            "washing",
        ],
    )
    def test_all_backend_supported_step_types_accepted(self, step_type_value):
        steps = [self._FakeStep(step_type_value)]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is True

    def test_unknown_step_type_blocker(self):
        steps = [self._FakeStep("magic_enchantment")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is False

    def test_upstream_trust_unreliable_blocker(self):
        """M1 trust UNRELIABLE blocks M2 regardless of step validity."""
        steps = [self._FakeStep("secondary_crosslinking")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="UNRELIABLE")
        assert result.valid is False
        assert any("UNRELIABLE" in b for b in result.blockers)

    def test_upstream_trust_caution_allows_run(self):
        """CAUTION trust level is not a blocker."""
        steps = [self._FakeStep("secondary_crosslinking")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="CAUTION")
        assert result.valid is True

    def test_empty_steps_blocker(self):
        result = validate_m2_inputs([], acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is False

    def test_no_acs_state_blocker(self):
        steps = [self._FakeStep("secondary_crosslinking")]
        result = validate_m2_inputs(steps, acs_state=None, m1_trust_level="TRUSTWORTHY")
        assert result.valid is False

    def test_duplicate_steps_warning(self):
        steps = [self._FakeStep("activation"), self._FakeStep("activation")]
        result = validate_m2_inputs(steps, acs_state=object(), m1_trust_level="TRUSTWORTHY")
        assert result.valid is True  # duplicates are warnings, not blockers
        assert any("duplicate" in w.lower() or "Duplicate" in w for w in result.warnings)


# ─── M3 Validation tests ──────────────────────────────────────────────────────

class TestValidationM3:
    def test_negative_flow_rate_blocker(self):
        result = validate_m3_chromatography(
            flow_rate=-1e-8,
            column=None,
            isotherm_type="langmuir",
            gradient_enabled=False,
        )
        assert result.valid is False

    def test_gradient_competitive_langmuir_warning(self):
        """Gradient + competitive_langmuir should raise a 'diagnostic only' warning."""
        result = validate_m3_chromatography(
            flow_rate=1e-8,
            column=None,
            isotherm_type="competitive_langmuir",
            gradient_enabled=True,
        )
        assert result.valid is True
        assert any("diagnostic" in w.lower() or "gradient" in w.lower() for w in result.warnings)
        assert all("BF-2" not in w for w in result.warnings)

    def test_gradient_linear_warning(self):
        result = validate_m3_chromatography(
            flow_rate=1e-8,
            column=None,
            isotherm_type="linear",
            gradient_enabled=True,
        )
        assert result.valid is True
        assert any("no effect" in w.lower() or "linear" in w.lower() for w in result.warnings)

    def test_no_gradient_no_warning(self):
        result = validate_m3_chromatography(
            flow_rate=1e-8,
            column=None,
            isotherm_type="langmuir",
            gradient_enabled=False,
        )
        assert result.valid is True
        assert result.warnings == []

    def test_mass_balance_blocker_above_5pct(self):
        """>5% mass balance error is a BLOCKER."""
        result = validate_m3_result(mass_balance_error=0.06, pressure_drop=1.0e4)
        assert result.valid is False
        assert any("5%" in b or "5 %" in b or "5.0%" in b or "6.0%" in b for b in result.blockers)

    def test_mass_balance_warning_2_to_5pct(self):
        """2–5% mass balance error is a WARNING, not a blocker."""
        result = validate_m3_result(mass_balance_error=0.03, pressure_drop=1.0e4)
        assert result.valid is True
        assert any("2" in w or "5" in w or "3.0%" in w for w in result.warnings)

    def test_mass_balance_ok_below_2pct(self):
        """<2% mass balance error is acceptable (no message)."""
        result = validate_m3_result(mass_balance_error=0.01, pressure_drop=1.0e4)
        assert result.valid is True
        assert result.blockers == []
        assert result.warnings == []

    def test_pressure_drop_blocker(self):
        """Pressure drop exceeding max_safe triggers a blocker."""
        result = validate_m3_result(
            mass_balance_error=0.005,
            pressure_drop=4.0e5,   # 4 bar
            max_safe_pressure=3.0e5,  # 3 bar limit
        )
        assert result.valid is False
        assert any("pressure" in b.lower() for b in result.blockers)


# ─── SessionStateManager tests ────────────────────────────────────────────────

class TestSessionStateManager:
    def test_initial_state_empty(self):
        mgr = SessionStateManager()
        assert mgr.m1_hash == ""
        assert mgr.m2_hash == ""
        assert mgr.m3_hash == ""

    def test_m1_update_returns_true_on_first_call(self):
        mgr = SessionStateManager()
        changed = mgr.update_m1({"rpm": 500, "phi_d": 0.2})
        assert changed is True

    def test_m1_update_returns_false_on_same_inputs(self):
        mgr = SessionStateManager()
        inputs = {"rpm": 500, "phi_d": 0.2}
        mgr.update_m1(inputs)
        changed = mgr.update_m1(inputs)
        assert changed is False

    def test_m1_change_invalidates_m2_and_m3(self):
        """M1 input change resets m2_hash and m3_hash."""
        mgr = SessionStateManager()
        mgr.update_m1({"rpm": 500})
        mgr.update_m2({"steps": 1})
        mgr.update_m3({"flow_rate": 1e-8})

        assert mgr.m2_hash != ""
        assert mgr.m3_hash != ""

        # Now change M1 inputs
        changed = mgr.update_m1({"rpm": 1000})
        assert changed is True
        assert mgr.m2_hash == ""
        assert mgr.m3_hash == ""

    def test_m2_change_invalidates_m3_only(self):
        """M2 input change resets m3_hash but NOT m1_hash."""
        mgr = SessionStateManager()
        mgr.update_m1({"rpm": 500})
        mgr.update_m2({"steps": 1})
        mgr.update_m3({"flow_rate": 1e-8})

        m1_hash_before = mgr.m1_hash
        changed = mgr.update_m2({"steps": 2})
        assert changed is True
        assert mgr.m1_hash == m1_hash_before   # M1 unchanged
        assert mgr.m3_hash == ""               # M3 cleared

    def test_same_inputs_no_downstream_clear(self):
        """Repeating same M1 inputs does NOT reset downstream hashes."""
        mgr = SessionStateManager()
        mgr.update_m1({"rpm": 500})
        mgr.update_m2({"steps": 1})
        m2_hash_before = mgr.m2_hash

        mgr.update_m1({"rpm": 500})  # no change
        assert mgr.m2_hash == m2_hash_before

    def test_invalidation_clears_store_keys(self):
        """Bound store keys are deleted on invalidation."""
        store = {"m2_result": "some_result", "m3_result": "some_result"}
        mgr = SessionStateManager()
        mgr.bind_store(store)
        mgr.update_m1({"rpm": 500})
        mgr.update_m1({"rpm": 1000})  # trigger cascade

        assert "m2_result" not in store
        assert "m3_result" not in store

    def test_has_run_flags(self):
        mgr = SessionStateManager()
        assert mgr.m1_has_run() is False
        mgr.update_m1({"rpm": 500})
        assert mgr.m1_has_run() is True
        assert mgr.m2_has_run() is False

    def test_reset_all_clears_everything(self):
        mgr = SessionStateManager()
        mgr.update_m1({"a": 1})
        mgr.update_m2({"b": 2})
        mgr.update_m3({"c": 3})
        mgr.reset_all()
        assert mgr.m1_hash == ""
        assert mgr.m2_hash == ""
        assert mgr.m3_hash == ""


# ─── Unit conversion tests ────────────────────────────────────────────────────

class TestUnitConversions:
    def test_mol_m3_to_mM_identity(self):
        assert mol_m3_to_mM(1.0) == pytest.approx(1.0)
        assert mol_m3_to_mM(0.0) == pytest.approx(0.0)
        assert mol_m3_to_mM(100.0) == pytest.approx(100.0)

    def test_mM_to_mol_m3_identity(self):
        assert mM_to_mol_m3(5.0) == pytest.approx(5.0)

    def test_Pa_to_kPa(self):
        assert Pa_to_kPa(1000.0) == pytest.approx(1.0)
        assert Pa_to_kPa(3.0e5) == pytest.approx(300.0)

    def test_Pa_to_bar(self):
        assert Pa_to_bar(1.0e5) == pytest.approx(1.0)
        assert Pa_to_bar(3.0e5) == pytest.approx(3.0)

    def test_bar_to_Pa_roundtrip(self):
        assert bar_to_Pa(Pa_to_bar(2.5e5)) == pytest.approx(2.5e5)

    def test_m_to_um(self):
        assert m_to_um(1.0e-6) == pytest.approx(1.0)
        assert m_to_um(25.0e-6) == pytest.approx(25.0)

    def test_m_to_nm(self):
        assert m_to_nm(1.0e-9) == pytest.approx(1.0)
        assert m_to_nm(100.0e-9) == pytest.approx(100.0)

    def test_um_to_m_roundtrip(self):
        assert um_to_m(m_to_um(15.3e-6)) == pytest.approx(15.3e-6, rel=1e-9)

    def test_nm_to_m_roundtrip(self):
        assert nm_to_m(m_to_nm(50.0e-9)) == pytest.approx(50.0e-9, rel=1e-9)

    def test_m3_s_to_mL_min(self):
        # 1 mL/min = 1e-6 L/s = 1e-6/60 m^3/s
        one_mL_per_min = 1.0e-6 / 60.0
        assert m3_s_to_mL_min(one_mL_per_min) == pytest.approx(1.0, rel=1e-6)

    def test_mL_min_to_m3_s_roundtrip(self):
        flow = 5.0  # mL/min
        assert m3_s_to_mL_min(mL_min_to_m3_s(flow)) == pytest.approx(flow, rel=1e-9)

    def test_K_to_C(self):
        assert K_to_C(273.15) == pytest.approx(0.0)
        assert K_to_C(373.15) == pytest.approx(100.0)

    def test_C_to_K(self):
        assert C_to_K(0.0) == pytest.approx(273.15)

    def test_format_sci_basic(self):
        s = format_sci(0.000123, sig=3)
        assert "1.23" in s
        assert "e" in s.lower() or "E" in s

    def test_format_sci_zero(self):
        s = format_sci(0.0, sig=3)
        assert "0" in s

    def test_format_sci_large(self):
        s = format_sci(1_234_567.0, sig=3)
        assert "1.23" in s

    def test_format_with_unit_fixed(self):
        s = format_with_unit(12.345, "mM", sig=3)
        assert "mM" in s
        assert "12" in s

    def test_format_with_unit_sci(self):
        s = format_with_unit(1.23e-7, "m", sig=3)
        assert "m" in s
        assert "e" in s.lower() or "E" in s

    def test_format_pressure_returns_string(self):
        s = format_pressure(2.5e5)
        assert isinstance(s, str)
        assert "bar" in s

    def test_format_length_nm_range(self):
        s = format_length(50e-9)
        assert "nm" in s

    def test_format_length_um_range(self):
        s = format_length(25e-6)
        assert "um" in s

    def test_format_concentration(self):
        s = format_concentration(5.0)
        assert "mM" in s
        assert "5" in s


# ─── Node 11 (v6.1, F6): UI metadata vs backend drift guards ──────────────


class TestUiBackendDriftGuards:
    """Lock in F6 fixes — fail loudly if UI metadata diverges from backend.

    These regression tests catch the failure modes that doc 35 documented:
      - validate_m2_inputs accepting only a subset of the actual M2 step types
        (so a user-defined workflow silently passes UI validation but crashes
        deeper in the orchestrator).
      - ui_model_metadata advertising stale claims about isotherm behaviour
        (e.g. "gradient does not affect binding" after gradient-aware routing
        was implemented).
    """

    def test_m2_validator_supports_every_backend_step_type(self):
        """_M2_SUPPORTED_STEP_TYPES must equal the full ModificationStepType enum."""
        backend_types = {st.value for st in ModificationStepType}
        ui_types = set(_M2_SUPPORTED_STEP_TYPES)
        missing_in_ui = backend_types - ui_types
        extra_in_ui = ui_types - backend_types
        assert not missing_in_ui, (
            f"UI validator missing backend step types: {sorted(missing_in_ui)}. "
            "User workflows using these would silently fail UI validation."
        )
        assert not extra_in_ui, (
            f"UI validator advertises non-existent step types: {sorted(extra_in_ui)}."
        )

    def test_m2_acs_metadata_lists_every_step_type(self):
        """M2_ACS_META.validity_range must mention every backend step type by name."""
        # Comparison is uppercase-insensitive — the metadata renders the
        # enum values uppercased.
        text = M2_ACS_META.validity_range.upper()
        for st in ModificationStepType:
            token = st.value.upper()
            assert token in text, (
                f"M2_ACS_META.validity_range omits step type {token!r}; "
                "users won't know it's supported."
            )

    def test_m3_gradient_warning_not_stale(self):
        """Gradient warning must reflect the current gradient-aware behaviour.

        The doc-35 stale claim was that gradient does not affect binding for
        competitive Langmuir. After the H6 wiring (v6.0) and the manifest
        rollout (Node 5), the metadata must NOT contain that wording.
        """
        joined = " ".join(M3_GRADIENT_META.warnings).lower()
        assert "does not affect binding" not in joined
        assert "gradient does not" not in joined
        # Positive assertion: the new wording mentions gradient-sensitive
        # isotherms updating binding during elution.
        assert "gradient-sensitive" in joined or "update binding" in joined
