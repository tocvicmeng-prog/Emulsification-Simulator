"""Contextual reagent-detail link builder (v9.0, milestone M8).

Single source of truth for building the /reagent_detail?... query-param
URL used by per-family formulation modules for inline "View mechanism
& protocol" links.
"""

from __future__ import annotations

from urllib.parse import urlencode


def build_reagent_link(
    *,
    key: str,
    source: str,
    T_K: float | None = None,
    t_s: float | None = None,
    c_mM: float | None = None,
    pH: float | None = None,
) -> str:
    """Build a relative URL to the reagent_detail page with contextual params."""
    qp: dict = {"key": key, "source": source}
    if T_K is not None:
        qp["T"] = f"{T_K:.2f}"
    if t_s is not None:
        qp["t"] = f"{t_s:.0f}"
    if c_mM is not None:
        qp["c"] = f"{c_mM:.3f}"
    if pH is not None:
        qp["pH"] = f"{pH:.2f}"
    return f"/reagent_detail?{urlencode(qp)}"
