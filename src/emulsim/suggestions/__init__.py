"""Structured parameter-optimization suggestions with hyperlinkable derivation pages.

Each suggestion (e.g. "increase RPM", "adjust cooling rate") is produced by a
per-key module in this package. A module exports three callables:

    generate(ctx)         -> Suggestion | None
    derive_target(ctx)    -> TargetRange
    render_derivation(ctx, target) -> None   (Streamlit page body)

The M1 tab calls generate_all(ctx) after a run to obtain all applicable
Suggestions, renders each with an [📊] link, and the /suggestion_detail page
dispatches back to the right module's render_derivation().

Scope of suggestion keys (v9.2.0 big-bang scope per JCP-EMULSIM-DERIV-001):
    increase_rpm, decrease_rpm, adjust_cooling_rate,
    increase_crosslinker, reduce_polymer
"""

from __future__ import annotations

from typing import Protocol

from .types import Suggestion, SuggestionContext, TargetRange


class SuggestionModule(Protocol):
    """The three-callable contract every per-key module must implement."""

    def generate(self, ctx: SuggestionContext) -> Suggestion | None: ...
    def derive_target(self, ctx: SuggestionContext) -> TargetRange: ...
    def render_derivation(self, ctx: SuggestionContext, target: TargetRange) -> None: ...


def _lazy_import(name: str) -> SuggestionModule:
    # Lazy: keeps streamlit + matplotlib off the hot import path of the M1 tab.
    import importlib
    return importlib.import_module(f".{name}", package=__name__)  # type: ignore[return-value]


# Registry: URL key → module. Multi-key entries (increase_rpm vs decrease_rpm)
# share a module but branch inside generate().
REGISTRY_KEYS: dict[str, str] = {
    "adjust_cooling_rate": "cooling_rate",
    "increase_rpm": "rpm",
    "decrease_rpm": "rpm",
    "increase_crosslinker": "crosslinker",
    "reduce_polymer": "polymer",
}


def get_module(key: str) -> SuggestionModule:
    """Fetch the module implementing a given suggestion key. Raises KeyError."""
    module_name = REGISTRY_KEYS[key]
    return _lazy_import(module_name)


def generate_all(ctx: SuggestionContext) -> list[Suggestion]:
    """Run every registered suggestion generator against the context.

    Order is stable: matches REGISTRY_KEYS insertion order. Modules that
    return None are filtered out (suggestion not applicable).
    """
    suggestions: list[Suggestion] = []
    seen_modules: set[str] = set()
    # RPM module is shared between increase/decrease keys — deduplicate so the
    # module runs once and internally decides which branch (if any) to emit.
    for key, mod_name in REGISTRY_KEYS.items():
        if mod_name in seen_modules:
            continue
        seen_modules.add(mod_name)
        module = _lazy_import(mod_name)
        result = module.generate(ctx)
        if result is not None:
            suggestions.append(result)
    return suggestions


__all__ = [
    "Suggestion",
    "SuggestionContext",
    "SuggestionModule",
    "TargetRange",
    "REGISTRY_KEYS",
    "generate_all",
    "get_module",
]
