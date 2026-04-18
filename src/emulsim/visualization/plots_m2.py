"""Module 2 visualization helpers for the EmulSim UI.

Provides:
  - plot_acs_waterfall:            Stacked bar chart of ACS site states per step.
  - plot_surface_area_comparison:  Bar chart of surface area breakdown.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    _PLOTLY_AVAILABLE = True
except ImportError:
    _PLOTLY_AVAILABLE = False

if TYPE_CHECKING:
    pass


def plot_acs_waterfall(
    history: list,
) -> Optional[object]:
    """Stacked bar chart showing ACS site states (remaining vs consumed) per step.

    Each bar represents one modification step. The bar is split into:
      - Remaining sites (after the step)
      - Sites consumed by this step

    Args:
        history: List of ModificationResult objects from ModificationOrchestrator.

    Returns:
        Plotly Figure, or None if plotly is unavailable or history is empty.
    """
    if not _PLOTLY_AVAILABLE or not history:
        return None

    step_labels = []
    remaining_vals = []
    consumed_vals = []

    for i, mr in enumerate(history):
        reagent = mr.step.reagent_key
        tgt = mr.step.target_acs
        label = f"Step {i + 1}: {reagent}"
        step_labels.append(label)

        before_profile = mr.acs_before.get(tgt)
        after_profile = mr.acs_after.get(tgt)

        if before_profile is not None and after_profile is not None:
            consumed = before_profile.remaining_sites - after_profile.remaining_sites
            remaining = after_profile.remaining_sites
        elif before_profile is not None:
            consumed = mr.conversion * before_profile.remaining_sites
            remaining = before_profile.remaining_sites - consumed
        else:
            consumed = 0.0
            remaining = 0.0

        # Convert to nmol/particle for readability
        consumed_vals.append(consumed * 1e9)
        remaining_vals.append(remaining * 1e9)

    fig = go.Figure()

    fig.add_trace(go.Bar(
        name="Sites Consumed",
        x=step_labels,
        y=consumed_vals,
        marker_color="tomato",
    ))

    fig.add_trace(go.Bar(
        name="Sites Remaining",
        x=step_labels,
        y=remaining_vals,
        marker_color="steelblue",
    ))

    fig.update_layout(
        barmode="stack",
        title="ACS Site States per Modification Step",
        xaxis_title="Modification Step",
        yaxis_title="Sites (nmol/particle)",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        height=380,
        template="plotly_white",
    )

    return fig


def plot_surface_area_comparison(
    surface_model,
) -> Optional[object]:
    """Bar chart comparing external, internal geometric, and accessible surface areas.

    Args:
        surface_model: AccessibleSurfaceModel instance from FunctionalMicrosphere.

    Returns:
        Plotly Figure, or None if plotly is unavailable or surface_model is None.
    """
    if not _PLOTLY_AVAILABLE or surface_model is None:
        return None

    area_labels = [
        "External",
        "Internal (geometric)",
        "Reagent-accessible",
        "Ligand-accessible",
    ]
    # Convert m^2 to um^2 for readability (1 m^2 = 1e12 um^2)
    area_values_um2 = [
        surface_model.external_area * 1e12,
        surface_model.internal_geometric_area * 1e12,
        surface_model.reagent_accessible_area * 1e12,
        surface_model.ligand_accessible_area * 1e12,
    ]
    bar_colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=area_labels,
        y=area_values_um2,
        marker_color=bar_colors,
        text=[f"{v:.2e}" for v in area_values_um2],
        textposition="outside",
    ))

    fig.update_layout(
        title=f"Surface Area Breakdown (tier: {surface_model.tier.value})",
        xaxis_title="Area Type",
        yaxis_title="Area (um^2/particle)",
        height=380,
        template="plotly_white",
        showlegend=False,
    )

    _tier_trust = getattr(surface_model, "trust_level", "CAUTION")
    fig.add_annotation(
        text=f"Trust: {_tier_trust}",
        xref="paper", yref="paper",
        x=1.0, y=1.05,
        showarrow=False,
        font=dict(size=11, color="gray"),
    )

    return fig
