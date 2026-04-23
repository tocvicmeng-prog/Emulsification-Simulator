"""Pure-function generators for the M1 suggestion list.

Relocated from visualization/tabs/tab_m1.py:772-792. Each per-key module's
`generate(ctx)` calls into helpers here to compute deviations + build the
display text. Keeping the text generation separate from derivation lets the
M1 tab render the list without importing streamlit-page dependencies.
"""

from __future__ import annotations

from .types import SuggestionContext


# ── Deviation helpers ────────────────────────────────────────────────────

def d_deviation(ctx: SuggestionContext) -> float:
    """Normalized |d_actual - d_target| / d_target. Uses d32."""
    if ctx.target_d32 <= 0:
        return 0.0
    return abs(ctx.d32_actual - ctx.target_d32) / ctx.target_d32


def pore_deviation(ctx: SuggestionContext) -> float:
    if ctx.target_pore <= 0:
        return 0.0
    return abs(ctx.pore_actual - ctx.target_pore) / ctx.target_pore


def G_deviation(ctx: SuggestionContext) -> float:
    """Log-distance between actual and target G; lets G span orders of magnitude."""
    import math
    if ctx.G_DN_actual <= 0 or ctx.target_G <= 0:
        return 0.0
    return abs(math.log10(ctx.G_DN_actual) - math.log10(ctx.target_G))


# ── Deviation thresholds (unchanged from tab_m1.py:773,779,782) ────────

DEVIATION_TRIGGER: float = 0.5  # match the existing UX threshold


# ── Display-text helpers (match the voice of the old `recs` list) ─────────

def rpm_up_text(ctx: SuggestionContext) -> str:
    return (
        f"**Increase RPM** (currently {ctx.rpm:.0f}) — droplet d32 "
        f"{ctx.d32_actual*1e6:.1f} µm > target {ctx.target_d32*1e6:.0f} µm. "
        f"Try higher RPM or increase Span-80 concentration."
    )


def rpm_down_text(ctx: SuggestionContext) -> str:
    return (
        f"**Decrease RPM** (currently {ctx.rpm:.0f}) — droplet d32 "
        f"{ctx.d32_actual*1e6:.1f} µm < target {ctx.target_d32*1e6:.0f} µm."
    )


def cooling_rate_text(ctx: SuggestionContext) -> str:
    direction = "Slower cooling → larger pores; faster → finer"
    return (
        f"**Adjust cooling rate** — pore size "
        f"{ctx.pore_actual*1e9:.0f} nm deviates from target "
        f"{ctx.target_pore*1e9:.0f} nm. {direction}."
    )


def crosslinker_up_text(ctx: SuggestionContext) -> str:
    return (
        f"**Increase crosslinker concentration** or crosslinking time — "
        f"G_DN ({ctx.G_DN_actual/1000:.1f} kPa) below target "
        f"({ctx.target_G/1000:.1f} kPa). Current is stoichiometry-limited "
        f"at p={ctx.p_final:.1%}."
    )


def polymer_down_text(ctx: SuggestionContext) -> str:
    return (
        f"**Reduce polymer concentration** — G_DN "
        f"({ctx.G_DN_actual/1000:.1f} kPa) exceeds target "
        f"({ctx.target_G/1000:.1f} kPa)."
    )


__all__ = [
    "d_deviation",
    "pore_deviation",
    "G_deviation",
    "DEVIATION_TRIGGER",
    "rpm_up_text",
    "rpm_down_text",
    "cooling_rate_text",
    "crosslinker_up_text",
    "polymer_down_text",
]
