"""Plotly-based visualization functions for simulation results."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..datatypes import EmulsificationResult, GelationResult, CrosslinkingResult, MechanicalResult, FullResult


def plot_droplet_size_distribution(result: EmulsificationResult) -> go.Figure:
    """Volume-weighted droplet size distribution from Level 1."""
    d_bins_safe = np.maximum(result.d_bins * 1e6, 1e-6)
    vol = result.n_d * (np.pi / 6.0 * result.d_bins**3)
    vol_pct = vol / np.sum(vol) * 100 if np.sum(vol) > 0 else vol

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=d_bins_safe,
        y=vol_pct,
        marker_color="steelblue",
        name="Volume %",
    ))
    fig.add_vline(x=result.d32 * 1e6, line_dash="dash", line_color="red",
                  annotation_text=f"d32={result.d32*1e6:.1f} µm")
    fig.add_vline(x=result.d50 * 1e6, line_dash="dot", line_color="green",
                  annotation_text=f"d50={result.d50*1e6:.1f} µm")
    fig.update_layout(
        title="Droplet Size Distribution (L1)",
        xaxis_title="Diameter (µm)", xaxis_type="log",
        yaxis_title="Volume %",
        height=400,
    )
    return fig


def plot_phase_field(result: GelationResult) -> go.Figure:
    """2D phase-field snapshot from Level 2."""
    phi = result.phi_field
    if phi.ndim == 1:
        fig = go.Figure(go.Scatter(
            x=result.r_grid * 1e9, y=phi, mode="lines", line_color="darkgreen",
        ))
        fig.update_layout(
            title="Radial Composition Profile (L2, 1D)",
            xaxis_title="r (nm)", yaxis_title="φ (polymer vol. frac.)",
            height=400,
        )
    else:
        r = result.r_grid
        fig = go.Figure(go.Heatmap(
            z=phi, x=r * 1e9, y=r * 1e9,
            colorscale="Viridis", colorbar_title="φ",
        ))
        fig.update_layout(
            title=f"Polymer Volume Fraction (L2, {phi.shape[0]}×{phi.shape[1]})",
            xaxis_title="x (nm)", yaxis_title="y (nm)",
            height=500, width=550,
        )
    return fig


def plot_crosslinking_kinetics(result: CrosslinkingResult) -> go.Figure:
    """Crosslinking evolution over time from Level 3."""
    t_h = result.t_array / 3600.0

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                        subplot_titles=("Crosslink Density & Modulus", "Mesh Size"))

    fig.add_trace(go.Scatter(
        x=t_h, y=result.G_chitosan_array,
        name="G_chitosan (Pa)", line_color="crimson",
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=t_h, y=result.xi_array * 1e9,
        name="ξ (nm)", line_color="teal",
    ), row=2, col=1)

    fig.update_xaxes(title_text="Time (hours)", row=2, col=1)
    fig.update_yaxes(title_text="G (Pa)", row=1, col=1)
    fig.update_yaxes(title_text="Mesh size (nm)", row=2, col=1)
    fig.update_layout(title="Crosslinking Kinetics (L3)", height=500)
    return fig


def plot_hertz_contact(result: MechanicalResult) -> go.Figure:
    """Force-displacement curve from Level 4."""
    fig = go.Figure(go.Scatter(
        x=result.delta_array * 1e9, y=result.F_array * 1e9,
        mode="lines", line_color="darkorange", name="Hertz",
    ))
    fig.update_layout(
        title=f"Microsphere Compression (G_DN={result.G_DN:.0f} Pa)",
        xaxis_title="Indentation (nm)", yaxis_title="Force (nN)",
        height=400,
    )
    return fig


def plot_kav_curve(result: MechanicalResult) -> go.Figure:
    """Partition coefficient vs solute size from Level 4."""
    fig = go.Figure(go.Scatter(
        x=result.rh_array * 1e9, y=result.Kav_array,
        mode="lines", line_color="purple",
    ))
    fig.update_layout(
        title="Chromatographic Partitioning (Ogston Model)",
        xaxis_title="Solute hydrodynamic radius (nm)",
        yaxis_title="K_av",
        xaxis_type="log",
        height=400,
    )
    return fig


def plot_results_dashboard(result: FullResult) -> go.Figure:
    """Combined 4-panel dashboard of all levels."""
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            f"L1: Droplet Size (d32={result.emulsification.d32*1e6:.1f} µm)",
            f"L2: Phase Field (pore={result.gelation.pore_size_mean*1e9:.0f} nm)",
            f"L3: Crosslinking (p={result.crosslinking.p_final:.3f})",
            f"L4: Compression (G_DN={result.mechanical.G_DN:.0f} Pa)",
        ),
    )

    # L1: Size distribution
    e = result.emulsification
    vol = e.n_d * (np.pi / 6.0 * e.d_bins**3)
    vol_norm = vol / np.max(vol) if np.max(vol) > 0 else vol
    d_bins_safe_dash = np.maximum(e.d_bins * 1e6, 1e-6)
    fig.add_trace(go.Bar(x=d_bins_safe_dash, y=vol_norm, marker_color="steelblue",
                         showlegend=False), row=1, col=1)

    # L2: Phase field (1D or diagonal slice of 2D)
    g = result.gelation
    phi = g.phi_field
    if phi.ndim == 2:
        diag = np.diag(phi)
        r_diag = g.r_grid[:len(diag)]
        fig.add_trace(go.Scatter(x=r_diag * 1e9, y=diag, mode="lines",
                                 line_color="darkgreen", showlegend=False), row=1, col=2)
    else:
        fig.add_trace(go.Scatter(x=g.r_grid * 1e9, y=phi, mode="lines",
                                 line_color="darkgreen", showlegend=False), row=1, col=2)

    # L3: Modulus evolution
    x = result.crosslinking
    fig.add_trace(go.Scatter(x=x.t_array / 3600, y=x.G_chitosan_array,
                             mode="lines", line_color="crimson", showlegend=False), row=2, col=1)

    # L4: Force-displacement
    m = result.mechanical
    fig.add_trace(go.Scatter(x=m.delta_array * 1e9, y=m.F_array * 1e9,
                             mode="lines", line_color="darkorange", showlegend=False), row=2, col=2)

    fig.update_xaxes(title_text="d (µm)", type="log", row=1, col=1)
    fig.update_xaxes(title_text="position (nm)", row=1, col=2)
    fig.update_xaxes(title_text="Time (h)", row=2, col=1)
    fig.update_xaxes(title_text="δ (nm)", row=2, col=2)
    fig.update_layout(height=700, title_text="Simulation Results Dashboard")
    return fig
