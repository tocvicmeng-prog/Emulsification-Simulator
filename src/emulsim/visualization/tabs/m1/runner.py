"""M1 pipeline runner (v9.0, milestones M3-M7).

Assembles SimulationParameters from the per-section contexts and calls
PipelineOrchestrator.run_single. Single writer of the
FullResult into session state.

Contract:
    family_context + hardware_context + formulation_context
    [+ crosslinking_context if AGAROSE_CHITOSAN]
    [+ targets_context, material_overrides]
      → SimulationParameters
      → orch.run_single(...)
      → FullResult in st.session_state["result"]

Status: SKELETON (populated progressively M3 → M7).
"""

from __future__ import annotations

from typing import Any


def run_m1_pipeline(
    *,
    family_context: Any,
    hardware_context: Any,
    formulation_context: Any,
    crosslinking_context: Any | None,
    targets_context: Any,
    material_overrides: dict,
) -> Any:
    """Execute the L1→L2→L3→L4 pipeline with the assembled contexts.

    Returns the FullResult. Raises any exception from the orchestrator
    so the caller can surface it via st.exception (see tab_m1.py:445).
    """
    raise NotImplementedError("Populated in M3 (skeleton) and fleshed out in M4-M7")
