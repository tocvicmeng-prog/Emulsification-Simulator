"""Uncertainty configuration panel for EmulSim UI.

Node 30b (v7.1): builds a :class:`UnifiedUncertaintySpec` from streamlit
inputs and persists it at ``st.session_state["_unc_spec"]``. Also
exposes the parallel-workers setting at ``st.session_state["_unc_n_jobs"]``.

The M1-L4 Monte Carlo (``UnifiedUncertaintyEngine.run_m1l4``) always
applies its built-in MaterialProperties perturbation set (IFT, viscosity,
kinetic prefactors, L2 pore coefficients, model-form bridge + IPN
coupling); the panel below configures the *additional* spec-driven
sources (custom user CVs) plus the sampling controls (n_samples, seed,
n_jobs) and surfaces a count of calibration-posterior sources that will
be absorbed from ``st.session_state["_cal_store"]`` at engine-construction
time.

The core helper :func:`build_uncertainty_spec` is a pure function of
its inputs so unit tests can verify spec construction without starting
streamlit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import streamlit as st

from emulsim.uncertainty_unified import (
    UncertaintyKind,
    UncertaintySource,
    UnifiedUncertaintySpec,
)


# ─── Pure helper (testable without streamlit) ─────────────────────────────


@dataclass(frozen=True)
class CustomSourceInput:
    """One user-defined uncertainty source from the panel's advanced expander."""
    name: str
    kind: UncertaintyKind
    value: float
    std: float
    distribution: str = "normal"  # "normal" | "lognormal"


def build_uncertainty_spec(
    *,
    n_samples: int,
    seed: int,
    custom_sources: list[CustomSourceInput] | None = None,
) -> UnifiedUncertaintySpec:
    """Build a UnifiedUncertaintySpec from panel inputs.

    Pure function — no streamlit calls, no session state. The engine
    absorbs calibration-store posteriors separately at construction, so
    the panel does not inject them into the spec here.

    Parameters
    ----------
    n_samples : int
        MC sample count. Must be >= 1.
    seed : int
        Master RNG seed.
    custom_sources : list of CustomSourceInput, optional
        Additional spec sources beyond the engine's built-in defaults.
        Each valid entry (non-empty name, std > 0) is appended.
    """
    if n_samples < 1:
        raise ValueError(f"n_samples must be >= 1, got {n_samples}")
    spec = UnifiedUncertaintySpec(n_samples=int(n_samples), seed=int(seed))
    for src in custom_sources or []:
        if not src.name.strip() or src.std <= 0:
            continue
        spec.add(UncertaintySource(
            name=src.name,
            kind=src.kind,
            distribution=src.distribution,
            value=float(src.value),
            std=float(src.std),
        ))
    return spec


def count_store_posteriors(store) -> int:
    """Count calibration entries with ``posterior_uncertainty > 0``.

    The UnifiedUncertaintyEngine absorbs one UncertaintySource per such
    entry at construction; the panel shows this count as user feedback.
    """
    if store is None:
        return 0
    total = 0
    for entry in getattr(store, "entries", []):
        sigma = float(getattr(entry, "posterior_uncertainty", 0.0) or 0.0)
        if sigma > 0.0:
            total += 1
    return total


# ─── Streamlit panel ──────────────────────────────────────────────────────


def render_uncertainty_panel() -> Optional[UnifiedUncertaintySpec]:
    """Render the uncertainty panel and return the configured spec.

    Side effects:
      - writes ``st.session_state["_unc_spec"]`` with the built spec
      - writes ``st.session_state["_unc_n_jobs"]`` with the selected
        parallel-workers value (integer; -1 means "all cores")
    """
    with st.expander("\U0001f4ca Uncertainty (MC sampling)", expanded=False):
        st.caption(
            "Configure Monte Carlo propagation for the M1\u2192L4 pipeline. "
            "The engine's built-in MaterialProperties perturbations are "
            "always applied; sources configured below are added on top."
        )

        # ── Sampling controls ──
        col1, col2, col3 = st.columns(3)
        with col1:
            n_samples = st.number_input(
                "Samples", min_value=10, max_value=5000, value=100, step=10,
                key="unc_n_samples",
                help="Monte Carlo sample count. Higher = tighter bounds, slower.",
            )
        with col2:
            seed = st.number_input(
                "Seed", min_value=0, max_value=2**31 - 1, value=42, step=1,
                key="unc_seed",
                help="RNG seed for reproducibility.",
            )
        with col3:
            n_jobs_choice = st.selectbox(
                "Parallelism",
                ["1 (serial)", "2", "4", "-1 (all cores)"],
                index=0,
                key="unc_n_jobs_choice",
                help=(
                    "Joblib workers. Auto-falls back to serial when "
                    "n_samples < 4 * n_jobs (Audit N7)."
                ),
            )
        _n_jobs_map = {"1 (serial)": 1, "2": 2, "4": 4, "-1 (all cores)": -1}
        n_jobs = _n_jobs_map[n_jobs_choice]
        st.session_state["_unc_n_jobs"] = n_jobs

        # ── Calibration posterior summary ──
        store = st.session_state.get("_cal_store")
        posterior_count = count_store_posteriors(store)
        if posterior_count > 0:
            st.info(
                f"\U0001f4cc {posterior_count} calibration posterior(s) "
                "will be absorbed from the calibration store at "
                "engine construction."
            )
        else:
            st.caption(
                "No calibration posteriors configured. Add entries with "
                "`posterior_uncertainty > 0` to the calibration store to "
                "propagate them through the MC."
            )

        # ── Advanced: custom uncertainty sources ──
        custom_sources: list[CustomSourceInput] = []
        with st.expander(
            "Advanced: custom uncertainty sources",
            expanded=False,
        ):
            st.caption(
                "Additional UncertaintySource entries appended to the "
                "spec. Use for measurement CVs or literature-spread "
                "material properties that the defaults do not cover."
            )
            n_custom = st.number_input(
                "Number of custom sources", 0, 10, 0, step=1,
                key="unc_n_custom",
            )
            _kind_values = [k.value for k in UncertaintyKind]
            for i in range(int(n_custom)):
                st.markdown(f"**Source {i + 1}**")
                row1 = st.columns([2, 1])
                with row1[0]:
                    name = st.text_input(
                        "Name",
                        value=f"custom_{i + 1}",
                        key=f"unc_custom_{i}_name",
                        help="Identifier (e.g. 'props.sigma').",
                    )
                with row1[1]:
                    kind_str = st.selectbox(
                        "Kind",
                        _kind_values,
                        index=_kind_values.index(
                            UncertaintyKind.MATERIAL_PROPERTY.value
                        ),
                        key=f"unc_custom_{i}_kind",
                    )
                row2 = st.columns([1, 1, 1])
                with row2[0]:
                    dist = st.selectbox(
                        "Distribution",
                        ["normal", "lognormal"],
                        index=0,
                        key=f"unc_custom_{i}_dist",
                    )
                with row2[1]:
                    value = st.number_input(
                        "Value",
                        value=1.0,
                        format="%.4g",
                        key=f"unc_custom_{i}_value",
                        help="Mean (normal) or median (lognormal).",
                    )
                with row2[2]:
                    std = st.number_input(
                        "Std",
                        min_value=0.0,
                        value=0.1,
                        format="%.4g",
                        key=f"unc_custom_{i}_std",
                        help="Std (normal) or sigma_log (lognormal).",
                    )
                custom_sources.append(CustomSourceInput(
                    name=name,
                    kind=UncertaintyKind(kind_str),
                    value=float(value),
                    std=float(std),
                    distribution=dist,
                ))

        # ── Build spec and persist ──
        spec = build_uncertainty_spec(
            n_samples=int(n_samples),
            seed=int(seed),
            custom_sources=custom_sources,
        )
        st.session_state["_unc_spec"] = spec

        # ── Summary footer ──
        extra = len(spec.sources)
        if extra:
            st.caption(
                f"Spec: n_samples={spec.n_samples}, seed={spec.seed}, "
                f"n_jobs={n_jobs}, +{extra} custom source(s) "
                f"+{posterior_count} posterior(s) absorbed at engine init."
            )
        else:
            st.caption(
                f"Spec: n_samples={spec.n_samples}, seed={spec.seed}, "
                f"n_jobs={n_jobs}, default perturbations only "
                f"+{posterior_count} posterior(s) absorbed at engine init."
            )

    return spec
