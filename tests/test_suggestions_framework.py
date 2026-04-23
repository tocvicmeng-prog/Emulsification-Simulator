"""Tests for the suggestions package foundation (types, registry, URL codec)."""

from __future__ import annotations

import pytest

from emulsim.suggestions import REGISTRY_KEYS, generate_all, get_module
from emulsim.suggestions.serialization import (
    ctx_from_query_params,
    ctx_to_query_dict,
    extras_from_query_params,
    suggestion_to_url,
    url_length_ok,
)
from emulsim.suggestions.types import Suggestion, SuggestionContext, TargetRange


def _baseline_ctx(**overrides: object) -> SuggestionContext:
    """Default SuggestionContext used across tests."""
    defaults: dict[str, object] = dict(
        family="agarose_chitosan",
        d32_actual=2e-6,
        d50_actual=3e-6,
        pore_actual=100e-9,
        l2_mode="ch_2d",
        cooling_rate_effective=0.17,
        p_final=0.25,
        G_DN_actual=10_000.0,
        target_d32=2e-6,
        target_pore=100e-9,
        target_G=10_000.0,
        rpm=8000.0,
        T_oil=363.15,
        cooling_rate_input=0.17,
        c_agarose=42.0,
        c_chitosan=18.0,
        c_crosslinker_mM=2.0,
        crosslinker_key="genipin",
        rho_oil=850.0,
        mu_oil=0.005,
        rho_d=1000.0,
        cp_d=4180.0,
        k_oil=0.15,
        h_coeff=500.0,
        T_bath=293.15,
        T_gel=311.15,
        DDA=0.85,
        M_GlcN=161.16,
        f_bridge=0.4,
        impeller_D=0.05,
        phi_d=0.05,
        run_id="test-run-001",
    )
    defaults.update(overrides)
    return SuggestionContext(**defaults)  # type: ignore[arg-type]


class TestRegistry:
    def test_registry_keys_are_stable(self):
        expected = {
            "adjust_cooling_rate",
            "increase_rpm",
            "decrease_rpm",
            "increase_crosslinker",
            "reduce_polymer",
        }
        assert set(REGISTRY_KEYS) == expected

    def test_every_registry_key_resolves_to_a_module(self):
        for key in REGISTRY_KEYS:
            module = get_module(key)
            assert hasattr(module, "generate")
            assert hasattr(module, "derive_target")
            assert hasattr(module, "render_derivation")

    def test_unknown_key_raises(self):
        with pytest.raises(KeyError):
            get_module("this_key_does_not_exist")


class TestURLCodec:
    def test_ctx_round_trip_is_identity(self):
        ctx = _baseline_ctx()
        qd = ctx_to_query_dict(ctx)
        ctx2 = ctx_from_query_params(qd)
        assert ctx == ctx2

    def test_round_trip_preserves_floats_exactly(self):
        # Use a value where a lossy string conversion would bite.
        ctx = _baseline_ctx(rho_oil=853.1234567890)
        qd = ctx_to_query_dict(ctx)
        ctx2 = ctx_from_query_params(qd)
        assert ctx.rho_oil == ctx2.rho_oil

    def test_missing_key_raises_key_error(self):
        ctx = _baseline_ctx()
        qd = ctx_to_query_dict(ctx)
        del qd["rho_oil"]
        with pytest.raises(KeyError):
            ctx_from_query_params(qd)

    def test_suggestion_to_url_includes_key_and_extras(self):
        ctx = _baseline_ctx()
        s = Suggestion(
            key="increase_rpm",
            display_text="x",
            severity="warning",
            context=ctx,
            extras={"direction": "increase"},
        )
        url = suggestion_to_url(s)
        assert url.startswith("/suggestion_detail?")
        assert "key=increase_rpm" in url
        assert "x_direction=increase" in url

    def test_extras_round_trip(self):
        ctx = _baseline_ctx()
        s = Suggestion(
            key="decrease_rpm",
            display_text="x",
            severity="warning",
            context=ctx,
            extras={"direction": "decrease", "note": "foo"},
        )
        url = suggestion_to_url(s)
        # Emulate st.query_params parsing.
        from urllib.parse import parse_qs, urlparse
        qs = parse_qs(urlparse(url).query)
        flat = {k: v[0] for k, v in qs.items()}
        extras = extras_from_query_params(flat)
        assert extras == {"direction": "decrease", "note": "foo"}

    def test_baseline_url_fits_streamlit_limit(self):
        ctx = _baseline_ctx()
        s = Suggestion(key="adjust_cooling_rate", display_text="x", severity="warning", context=ctx)
        url = suggestion_to_url(s)
        assert url_length_ok(url), f"URL length {len(url)} exceeds the Streamlit-safe cap"


class TestGenerateAll:
    def test_no_deviation_returns_empty_list(self):
        ctx = _baseline_ctx()  # all actuals match targets
        assert generate_all(ctx) == []

    def test_large_pore_triggers_cooling_rate(self):
        ctx = _baseline_ctx(pore_actual=300e-9)  # 3x target
        suggestions = generate_all(ctx)
        keys = [s.key for s in suggestions]
        assert "adjust_cooling_rate" in keys

    def test_large_d32_triggers_increase_rpm(self):
        ctx = _baseline_ctx(d32_actual=6e-6)  # 3x target
        suggestions = generate_all(ctx)
        assert any(s.key == "increase_rpm" for s in suggestions)

    def test_small_d32_triggers_decrease_rpm(self):
        ctx = _baseline_ctx(d32_actual=0.5e-6)  # quarter of target
        suggestions = generate_all(ctx)
        assert any(s.key == "decrease_rpm" for s in suggestions)

    def test_low_G_triggers_increase_crosslinker_only(self):
        ctx = _baseline_ctx(G_DN_actual=1_000.0)  # 10x below target
        suggestions = generate_all(ctx)
        keys = [s.key for s in suggestions]
        assert "increase_crosslinker" in keys
        assert "reduce_polymer" not in keys

    def test_high_G_triggers_reduce_polymer_only(self):
        ctx = _baseline_ctx(G_DN_actual=100_000.0)  # 10x above target
        suggestions = generate_all(ctx)
        keys = [s.key for s in suggestions]
        assert "reduce_polymer" in keys
        assert "increase_crosslinker" not in keys


class TestTargetRange:
    def test_is_qualitative_only_default_false(self):
        t = TargetRange(
            nominal=1.0, min=0.9, max=1.1,
            unit="K/s", limited_by="ok",
            confidence_tier="SEMI_QUANTITATIVE",
        )
        assert t.is_qualitative_only is False
        assert t.qualitative_reason == ""

    def test_qualitative_path(self):
        t = TargetRange(
            nominal=0.0, min=0.0, max=0.0,
            unit="K/s", limited_by="ok",
            confidence_tier="QUALITATIVE_TREND",
            is_qualitative_only=True,
            qualitative_reason="empirical L2 model",
        )
        assert t.is_qualitative_only is True
        assert "empirical" in t.qualitative_reason
