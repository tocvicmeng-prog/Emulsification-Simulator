"""Suggestion derivation page.

Reads the SuggestionContext + key from the URL query string, dispatches to
the per-key suggestion module's render_derivation(), and renders a uniform
three-section layout:

    1. Derivation logic and formula pathway
    2. Target value or range
    3. Assumptions and confidence

See docs/user_manual §8b for the user-facing explanation of these sections.
"""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

_root = Path(__file__).resolve().parents[4]
if str(_root / "src") not in sys.path:
    sys.path.insert(0, str(_root / "src"))

from emulsim.suggestions import REGISTRY_KEYS, get_module
from emulsim.suggestions.serialization import (
    ctx_from_query_params,
    extras_from_query_params,
)

st.set_page_config(
    page_title="Suggestion Derivation",
    page_icon="📊",
    layout="wide",
)

# ─── Query-param parsing ───────────────────────────────────────────────

params = dict(st.query_params)
key = params.get("key", "")

if key not in REGISTRY_KEYS:
    st.error(
        f"Unknown suggestion key '{key}'. Known keys: {sorted(REGISTRY_KEYS)}"
    )
    st.stop()

try:
    ctx = ctx_from_query_params(params)
except (KeyError, ValueError) as exc:
    st.error(
        "Could not parse the suggestion context from the URL. This page "
        "is meant to be opened via the [📊] link from the M1 tab after a run.\n\n"
        f"Parse error: {exc}"
    )
    st.stop()

extras = extras_from_query_params(params)

# ─── Header + back link ────────────────────────────────────────────────

module = get_module(key)
target = module.derive_target(ctx)

_title_map = {
    "adjust_cooling_rate": "Why: adjust the cooling rate",
    "increase_rpm":        "Why: increase RPM",
    "decrease_rpm":        "Why: decrease RPM",
    "increase_crosslinker": "Why: increase crosslinker concentration",
    "reduce_polymer":      "Why: reduce polymer concentration",
}

st.title(_title_map.get(key, f"Why: {key}"))
st.markdown("[← Back to Simulator](/)")
if ctx.run_id:
    st.caption(f"Derivation tied to run_id: `{ctx.run_id}`")

# ─── Section 1 — derivation ────────────────────────────────────────────

st.header("1. Derivation logic and formula pathway")
module.render_derivation(ctx, target)

# ─── Section 2 — target value / range ──────────────────────────────────

st.header("2. Target value or range")

if target.is_qualitative_only:
    st.warning("**Numeric target withheld — qualitative-only suggestion.**")
    st.markdown(target.qualitative_reason)
else:
    c1, c2, c3 = st.columns(3)
    c1.metric("Nominal", f"{target.nominal:.3f} {target.unit}")
    c2.metric("Lower bound", f"{target.min:.3f} {target.unit}")
    c3.metric("Upper bound", f"{target.max:.3f} {target.unit}")
    st.caption(
        f"Feasibility constraint: **{target.limited_by}**. "
        f"The band widens or narrows based on the tolerance allowed on the "
        f"target-of-interest — see Section 3 for the exact tolerances used."
    )

# ─── Section 3 — assumptions + confidence ──────────────────────────────

st.header("3. Assumptions and confidence")

_tier_color = {
    "VALIDATED":          "✅",
    "SEMI_QUANTITATIVE":  "🟡",
    "QUALITATIVE_TREND":  "🟠",
    "UNSUPPORTED":        "❌",
}
tier_icon = _tier_color.get(target.confidence_tier, "❔")
st.markdown(f"**Confidence tier:** {tier_icon} `{target.confidence_tier}`")

if target.assumptions:
    st.markdown("**Assumptions that went into this derivation:**")
    for a in target.assumptions:
        st.markdown(f"- {a}")
else:
    st.caption("No additional assumptions beyond those stated in Section 1.")
