"""URL round-trip codec for SuggestionContext.

Encodes a frozen dataclass to a query-string dict and back. Chosen over
session-state transport so derivation pages are bookmarkable and shareable.
Emits a warning in the caller if the encoded length would exceed Streamlit's
~2 KB URL limit, but proceeds.

Format: each field encoded as its string repr. Floats use repr() for
round-trip precision; strings are URL-safe quoted by urllib on the consumer.
"""

from __future__ import annotations

import urllib.parse
from dataclasses import fields

from .types import Suggestion, SuggestionContext


_STREAMLIT_SAFE_URL_LEN = 1800  # conservative under the ~2 KB ceiling


def ctx_to_query_dict(ctx: SuggestionContext) -> dict[str, str]:
    """Serialise every field to a flat {str: str} dict."""
    out: dict[str, str] = {}
    for f in fields(ctx):
        value = getattr(ctx, f.name)
        if isinstance(value, float):
            out[f.name] = repr(value)
        else:
            out[f.name] = str(value)
    return out


def ctx_from_query_params(params: dict[str, str]) -> SuggestionContext:
    """Rebuild a SuggestionContext from a Streamlit st.query_params-like dict.

    Missing keys raise KeyError, malformed floats raise ValueError — both
    surface to the page as an error message rather than silent defaults.
    """
    kwargs: dict[str, object] = {}
    for f in fields(SuggestionContext):
        raw = params[f.name]
        if f.type == "float" or f.type is float:
            kwargs[f.name] = float(raw)
        elif f.type == "int" or f.type is int:
            kwargs[f.name] = int(raw)
        elif f.type == "bool" or f.type is bool:
            kwargs[f.name] = raw in ("True", "true", "1")
        else:
            kwargs[f.name] = raw
    return SuggestionContext(**kwargs)  # type: ignore[arg-type]


def suggestion_to_url(suggestion: Suggestion, page_path: str = "/suggestion_detail") -> str:
    """Build the full href for a Suggestion's derivation page."""
    qd = ctx_to_query_dict(suggestion.context)
    qd["key"] = suggestion.key
    for k, v in suggestion.extras.items():
        qd[f"x_{k}"] = v
    query = urllib.parse.urlencode(qd, safe=",.-")
    return f"{page_path}?{query}"


def extras_from_query_params(params: dict[str, str]) -> dict[str, str]:
    """Pull out the x_-prefixed extras a Suggestion can stash in its URL."""
    return {k[2:]: v for k, v in params.items() if k.startswith("x_")}


def url_length_ok(url: str) -> bool:
    return len(url) <= _STREAMLIT_SAFE_URL_LEN


__all__ = [
    "ctx_to_query_dict",
    "ctx_from_query_params",
    "suggestion_to_url",
    "extras_from_query_params",
    "url_length_ok",
]
