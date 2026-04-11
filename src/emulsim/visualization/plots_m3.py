"""Plotly-based visualization functions for Module 3 (chromatography & catalysis) results."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..module3_performance.gradient import GradientProgram
from ..module3_performance.hydrodynamics import ColumnGeometry
from ..module3_performance.catalysis.kinetics import effectiveness_factor


# ─── 1. Chromatogram ─────────────────────────────────────────────────────────

def plot_chromatogram(
    time: np.ndarray,
    uv_signal: np.ndarray,
    gradient_values: np.ndarray | None = None,
    gradient_affects_binding: bool = False,
) -> go.Figure:
    """UV chromatogram with optional gradient overlay.

    Args:
        time: Time array [s].
        uv_signal: UV absorbance signal [mAU].
        gradient_values: Gradient profile (salt [M] or pH) vs time.
            If provided, plotted on a secondary y-axis as a dashed line.
        gradient_affects_binding: If True, annotates the gradient as
            mechanistically coupled to binding; otherwise marks it as
            diagnostic/display-only.

    Returns:
        go.Figure with primary UV trace and optional gradient overlay.
    """
    time_min = np.asarray(time, dtype=float) / 60.0
    uv = np.asarray(uv_signal, dtype=float)

    if gradient_values is not None:
        fig = make_subplots(specs=[[{"secondary_y": True}]])
    else:
        fig = go.Figure()

    # Primary: UV trace
    uv_trace = go.Scatter(
        x=time_min,
        y=uv,
        mode="lines",
        name="UV 280 nm",
        line=dict(color="steelblue", width=2),
    )

    if gradient_values is not None:
        fig.add_trace(uv_trace, secondary_y=False)

        grad = np.asarray(gradient_values, dtype=float)
        grad_trace = go.Scatter(
            x=time_min,
            y=grad,
            mode="lines",
            name="Gradient",
            line=dict(color="tomato", width=1.5, dash="dash"),
        )
        fig.add_trace(grad_trace, secondary_y=True)

        # Badge annotation
        if gradient_affects_binding:
            badge_text = "Gradient drives binding"
            badge_color = "rgba(255, 80, 80, 0.15)"
        else:
            badge_text = "Gradient: diagnostic only"
            badge_color = "rgba(100, 100, 100, 0.12)"

        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.01, y=0.97,
            text=badge_text,
            showarrow=False,
            bgcolor=badge_color,
            bordercolor="gray",
            borderwidth=1,
            font=dict(size=11),
            xanchor="left",
            yanchor="top",
        )

        fig.update_yaxes(title_text="UV Absorbance (mAU)", secondary_y=False)
        fig.update_yaxes(title_text="Gradient [M]", secondary_y=True,
                         showgrid=False)
    else:
        fig.add_trace(uv_trace)
        fig.update_yaxes(title_text="UV Absorbance (mAU)")

    fig.update_xaxes(title_text="Time (min)")
    fig.update_layout(
        title="Chromatogram",
        height=400,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


# ─── 2. Breakthrough Curve ───────────────────────────────────────────────────

def plot_breakthrough_curve(
    time: np.ndarray,
    C_outlet: np.ndarray,
    C_feed: float,
    dbc_5: float | None = None,
    dbc_10: float | None = None,
    dbc_50: float | None = None,
) -> go.Figure:
    """C/C0 breakthrough curve with threshold lines and DBC annotations.

    Args:
        time: Time array [s].
        C_outlet: Outlet concentration [mol/m^3].
        C_feed: Feed concentration [mol/m^3].
        dbc_5: Dynamic binding capacity at 5% breakthrough [mol/m^3 column].
        dbc_10: Dynamic binding capacity at 10% breakthrough [mol/m^3 column].
        dbc_50: Dynamic binding capacity at 50% breakthrough [mol/m^3 column].

    Returns:
        go.Figure showing normalised breakthrough profile.
    """
    time_min = np.asarray(time, dtype=float) / 60.0
    C_out = np.asarray(C_outlet, dtype=float)
    C0 = float(C_feed)

    cc0 = C_out / C0 if C0 > 0 else C_out

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=time_min,
        y=cc0,
        mode="lines",
        name="C/C\u2080",
        line=dict(color="steelblue", width=2),
    ))

    # Threshold horizontal lines
    thresholds = [(0.05, "5%", "green"), (0.10, "10%", "goldenrod"), (0.50, "50%", "tomato")]
    for frac, label, color in thresholds:
        fig.add_hline(
            y=frac,
            line_dash="dot",
            line_color=color,
            annotation_text=label,
            annotation_position="right",
            annotation_font_color=color,
        )

    # Vertical DBC annotations
    dbc_map = [
        (dbc_5, 0.05, "DBC\u2085", "green"),
        (dbc_10, 0.10, "DBC\u2081\u2080", "goldenrod"),
        (dbc_50, 0.50, "DBC\u2085\u2080", "tomato"),
    ]
    for dbc_val, frac, label, color in dbc_map:
        if dbc_val is not None and not np.isnan(dbc_val):
            # Find approximate breakthrough time from data
            above = np.where(cc0 >= frac)[0]
            if len(above) > 0:
                t_bt_min = time_min[above[0]]
                fig.add_vline(
                    x=t_bt_min,
                    line_dash="dash",
                    line_color=color,
                    opacity=0.6,
                    annotation_text=f"{label}={dbc_val:.1f} mol/m\u00b3",
                    annotation_position="top left",
                    annotation_font_color=color,
                    annotation_font_size=10,
                )

    fig.update_layout(
        title="Breakthrough Curve",
        xaxis_title="Time (min)",
        yaxis_title="C / C\u2080",
        yaxis=dict(range=[0, 1.05]),
        height=400,
    )
    return fig


# ─── 3. Peak Table ───────────────────────────────────────────────────────────

def plot_peak_table(peaks: list) -> go.Figure:
    """Plotly table of chromatographic peak statistics.

    Args:
        peaks: List of PeakInfo objects from GradientElutionResult.peaks.

    Returns:
        go.Figure containing a formatted table.
    """
    if not peaks:
        fig = go.Figure(go.Table(
            header=dict(values=["No peaks detected"]),
            cells=dict(values=[["—"]]),
        ))
        fig.update_layout(title="Peak Table", height=200)
        return fig

    components = [f"Component {p.component_idx + 1}" for p in peaks]
    tr_min = [
        f"{p.retention_time / 60.0:.2f}" if not np.isnan(p.retention_time) else "—"
        for p in peaks
    ]
    heights = [f"{p.peak_height:.3g}" for p in peaks]
    widths_s = [
        f"{p.peak_width_half:.1f}" if not np.isnan(p.peak_width_half) else "—"
        for p in peaks
    ]
    yields_pct = [f"{p.yield_fraction * 100:.1f}" for p in peaks]
    purities_pct = [f"{p.purity * 100:.1f}" for p in peaks]

    header_vals = ["Component", "tR (min)", "Height (mAU)", "Width (s)", "Yield (%)", "Purity (%)"]
    cell_vals = [components, tr_min, heights, widths_s, yields_pct, purities_pct]

    fig = go.Figure(go.Table(
        header=dict(
            values=header_vals,
            fill_color="steelblue",
            font=dict(color="white", size=12),
            align="center",
        ),
        cells=dict(
            values=cell_vals,
            fill_color=[["white", "aliceblue"] * (len(peaks) // 2 + 1)] * len(cell_vals),
            align="center",
            font=dict(size=11),
        ),
    ))
    fig.update_layout(
        title="Peak Summary Table",
        height=max(200, 80 + 30 * len(peaks)),
    )
    return fig


# ─── 4. Gradient Preview ─────────────────────────────────────────────────────

def plot_gradient_preview(gradient: GradientProgram, total_time: float) -> go.Figure:
    """Small line chart previewing the gradient program shape.

    Args:
        gradient: GradientProgram instance.
        total_time: Total simulation time [s] for x-axis range.

    Returns:
        Compact go.Figure (height=200) showing gradient value vs time.
    """
    t_eval = np.linspace(0.0, total_time, 200)
    g_vals = gradient.value_at_time(t_eval)
    t_min = t_eval / 60.0

    fig = go.Figure(go.Scatter(
        x=t_min,
        y=g_vals,
        mode="lines",
        line=dict(color="tomato", width=1.5),
        name="Gradient",
        fill="tozeroy",
        fillcolor="rgba(255, 99, 71, 0.12)",
    ))
    fig.update_layout(
        title="Gradient Preview",
        xaxis_title="Time (min)",
        yaxis_title="Gradient value",
        height=200,
        margin=dict(t=40, b=40, l=60, r=20),
        showlegend=False,
    )
    return fig


# ─── 5. Michaelis-Menten Curve ───────────────────────────────────────────────

def plot_michaelis_menten(
    S_range: np.ndarray,
    V_max: float,
    K_m: float,
    eta: float = 1.0,
) -> go.Figure:
    """Michaelis-Menten velocity curves with diffusion-limitation shading.

    Shows intrinsic (dashed) and effective (solid) reaction rate vs
    substrate concentration with a shaded zone between them when eta < 1.

    Args:
        S_range: Substrate concentration array [mol/m^3].
        V_max: Maximum reaction rate [mol/(m^3*s)].
        K_m: Michaelis constant [mol/m^3].
        eta: Effectiveness factor (0-1).  Default 1.0 (no limitation).

    Returns:
        go.Figure with dual-rate overlay.
    """
    S = np.asarray(S_range, dtype=float)
    v_intrinsic = V_max * S / (K_m + S)
    v_effective = eta * v_intrinsic

    fig = go.Figure()

    # Shaded diffusion-limitation zone between curves
    if eta < 1.0:
        fig.add_trace(go.Scatter(
            x=np.concatenate([S, S[::-1]]),
            y=np.concatenate([v_intrinsic, v_effective[::-1]]),
            fill="toself",
            fillcolor="rgba(255, 160, 0, 0.18)",
            line=dict(color="rgba(255,255,255,0)"),
            name="Diffusion limitation",
            hoverinfo="skip",
        ))

    # Intrinsic rate (dashed)
    fig.add_trace(go.Scatter(
        x=S,
        y=v_intrinsic,
        mode="lines",
        name="v intrinsic",
        line=dict(color="gray", width=2, dash="dash"),
    ))

    # Effective rate (solid)
    fig.add_trace(go.Scatter(
        x=S,
        y=v_effective,
        mode="lines",
        name=f"v effective (\u03b7={eta:.2f})",
        line=dict(color="darkorange", width=2),
    ))

    # K_m marker
    fig.add_vline(
        x=K_m,
        line_dash="dot",
        line_color="steelblue",
        annotation_text=f"K\u2098={K_m:.2g} mol/m\u00b3",
        annotation_position="top right",
    )

    fig.update_layout(
        title="Michaelis-Menten Kinetics",
        xaxis_title="[S] (mol/m\u00b3)",
        yaxis_title="Reaction rate (mol m\u207b\u00b3 s\u207b\u00b9)",
        height=400,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


# ─── 6. Effectiveness Factor vs Thiele Modulus ───────────────────────────────

def plot_effectiveness_factor(phi_range: np.ndarray) -> go.Figure:
    """Effectiveness factor vs Thiele modulus (log x-axis) with regime zones.

    Args:
        phi_range: Array of Thiele modulus values (dimensionless).

    Returns:
        go.Figure with colored kinetic/transition/diffusion zones.
    """
    phi = np.asarray(phi_range, dtype=float)
    phi_safe = np.clip(phi, 1e-6, None)
    eta = effectiveness_factor(phi_safe)

    fig = go.Figure()

    # Zone boundaries
    phi_min = float(phi.min())
    phi_max = float(phi.max())

    zone_defs = [
        (phi_min, 0.3, "rgba(0,180,0,0.10)", "Kinetic regime (\u03a6<0.3)"),
        (0.3, 3.0, "rgba(255,200,0,0.15)", "Transition"),
        (3.0, phi_max, "rgba(220,50,50,0.12)", "Diffusion-limited (\u03a6>3)"),
    ]
    for x0, x1, color, label in zone_defs:
        if x0 < phi_max and x1 > phi_min:
            x0_clip = max(x0, phi_min)
            x1_clip = min(x1, phi_max)
            fig.add_vrect(
                x0=x0_clip, x1=x1_clip,
                fillcolor=color,
                line_width=0,
                annotation_text=label,
                annotation_position="top left",
                annotation_font_size=10,
            )

    fig.add_trace(go.Scatter(
        x=phi,
        y=eta,
        mode="lines",
        name="\u03b7(\u03a6)",
        line=dict(color="steelblue", width=2),
    ))

    fig.update_layout(
        title="Effectiveness Factor vs Thiele Modulus",
        xaxis=dict(title="\u03a6 (Thiele modulus)", type="log"),
        yaxis=dict(title="\u03b7 (effectiveness factor)", range=[0, 1.05]),
        height=400,
    )
    return fig


# ─── 7. Activity Decay ───────────────────────────────────────────────────────

def plot_activity_decay(
    time_hours: np.ndarray,
    activity: np.ndarray,
) -> go.Figure:
    """Enzyme activity fraction vs time with half-life annotation.

    Args:
        time_hours: Time array [hours].
        activity: Remaining activity fraction [-], normalised to 1 at t=0.

    Returns:
        go.Figure with exponential decay trace and t_1/2 annotation.
    """
    t_h = np.asarray(time_hours, dtype=float)
    a = np.asarray(activity, dtype=float)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=t_h,
        y=a,
        mode="lines",
        name="Activity",
        line=dict(color="crimson", width=2),
        fill="tozeroy",
        fillcolor="rgba(220, 20, 60, 0.08)",
    ))

    # Half-life annotation: find t when activity ~ 0.5
    above_half = np.where(a <= 0.5)[0]
    if len(above_half) > 0:
        idx = above_half[0]
        if idx > 0 and a[idx - 1] > 0.5:
            frac = (0.5 - a[idx - 1]) / (a[idx] - a[idx - 1])
            t_half = t_h[idx - 1] + frac * (t_h[idx] - t_h[idx - 1])
        else:
            t_half = t_h[idx]

        fig.add_vline(
            x=t_half,
            line_dash="dash",
            line_color="gray",
            opacity=0.7,
            annotation_text=f"t\u2082= {t_half:.1f} h",
            annotation_position="top right",
        )
        fig.add_hline(
            y=0.5,
            line_dash="dot",
            line_color="gray",
            opacity=0.5,
        )

    fig.update_layout(
        title="Enzyme Activity Decay",
        xaxis_title="Time (hours)",
        yaxis=dict(title="Activity fraction [-]", range=[0, 1.05]),
        height=400,
        showlegend=False,
    )
    return fig


# ─── 8. Conversion vs Time ───────────────────────────────────────────────────

def plot_conversion_vs_time(
    time_hours: np.ndarray,
    conversion: np.ndarray,
) -> go.Figure:
    """Substrate conversion vs reactor time.

    Args:
        time_hours: Time array [hours].
        conversion: Fractional conversion [-], values in [0, 1].

    Returns:
        go.Figure with conversion trace.
    """
    t_h = np.asarray(time_hours, dtype=float)
    X = np.asarray(conversion, dtype=float)
    X = np.clip(X, 0.0, 1.0)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=t_h,
        y=X,
        mode="lines",
        name="Conversion",
        line=dict(color="teal", width=2),
    ))

    # Final conversion annotation
    X_final = float(X[-1]) if len(X) > 0 else 0.0
    fig.add_annotation(
        x=float(t_h[-1]) if len(t_h) > 0 else 0,
        y=X_final,
        text=f"X\u209c={X_final:.1%}",
        showarrow=True,
        arrowhead=2,
        ax=-40,
        ay=-30,
        font=dict(size=11),
    )

    fig.update_layout(
        title="Substrate Conversion vs Time",
        xaxis_title="Time (hours)",
        yaxis=dict(title="Conversion [-]", range=[0, 1.05]),
        height=400,
        showlegend=False,
    )
    return fig


# ─── 9. ESI-MS Spectrum ──────────────────────────────────────────────────────

def plot_esi_spectrum(
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    molecular_weight: float,
) -> go.Figure:
    """ESI-MS charge-envelope spectrum as a stem plot.

    Args:
        mz_array: m/z values [Da/e].
        intensity_array: Relative intensity [-], normalised to max=1.
        molecular_weight: Protein molecular weight [Da] for MW annotation.

    Returns:
        go.Figure with vertical stems at each m/z point.
    """
    mz = np.asarray(mz_array, dtype=float)
    inten = np.asarray(intensity_array, dtype=float)

    # Build stem-plot as vertical lines
    stem_x: list[float] = []
    stem_y: list[float] = []
    for xi, yi in zip(mz, inten):
        stem_x.extend([xi, xi, None])  # type: ignore[list-item]
        stem_y.extend([0.0, yi, None])  # type: ignore[list-item]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=stem_x,
        y=stem_y,
        mode="lines",
        name="ESI envelope",
        line=dict(color="steelblue", width=1),
    ))

    # Marker at peak tops
    peak_mask = inten > 0.05
    if peak_mask.any():
        fig.add_trace(go.Scatter(
            x=mz[peak_mask],
            y=inten[peak_mask],
            mode="markers",
            marker=dict(color="steelblue", size=5),
            showlegend=False,
        ))

    # MW annotation
    mw_kda = molecular_weight / 1000.0
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.99, y=0.97,
        text=f"MW = {mw_kda:.1f} kDa",
        showarrow=False,
        xanchor="right",
        yanchor="top",
        bgcolor="rgba(240,248,255,0.85)",
        bordercolor="steelblue",
        borderwidth=1,
        font=dict(size=11),
    )

    fig.update_layout(
        title="ESI-MS Spectrum (Charge Envelope)",
        xaxis_title="m/z (Da/e)",
        yaxis=dict(title="Relative intensity [-]", range=[0, 1.08]),
        height=400,
        showlegend=False,
    )
    return fig


# ─── 10. Pressure-Flow Curve ─────────────────────────────────────────────────

def plot_pressure_flow_curve(
    column: ColumnGeometry,
    mu: float = 1e-3,
) -> go.Figure:
    """Pressure drop vs flow rate with safe/unsafe zone annotations.

    Plots dP = Kozeny-Carman(Q) from zero to 1.5 × the maximum safe flow rate.
    A vertical dashed line marks Q_max with shaded safe/unsafe regions.

    Args:
        column: ColumnGeometry instance (provides particle size, bed height,
            porosity, and E_star for the safe-flow estimate).
        mu: Dynamic viscosity [Pa.s]. Default: water at 20 C (1e-3).

    Returns:
        go.Figure with dP [bar] vs Q [mL/min].
    """
    Q_max = column.max_safe_flow_rate(mu=mu)
    Q_plot_max = max(Q_max * 1.5, 1e-9)

    Q_arr = np.linspace(0.0, Q_plot_max, 200)
    dP_arr = np.array([column.pressure_drop(q, mu) for q in Q_arr])

    # Convert to display units: Q in mL/min, dP in bar
    Q_mlmin = Q_arr * 1e6 * 60.0     # m^3/s -> mL/min
    dP_bar = dP_arr / 1e5            # Pa -> bar
    Q_max_mlmin = Q_max * 1e6 * 60.0
    dP_max_bar = column.pressure_drop(Q_max, mu) / 1e5

    fig = go.Figure()

    # Safe zone fill (left of Q_max)
    fig.add_trace(go.Scatter(
        x=[0, Q_max_mlmin, Q_max_mlmin, 0],
        y=[0, dP_max_bar, 0, 0],
        fill="toself",
        fillcolor="rgba(0, 200, 100, 0.10)",
        line=dict(color="rgba(0,0,0,0)"),
        name="Safe zone",
        hoverinfo="skip",
    ))

    # Unsafe zone fill (right of Q_max, up to dP curve)
    Q_unsafe = Q_mlmin[Q_mlmin >= Q_max_mlmin]
    dP_unsafe = dP_bar[Q_mlmin >= Q_max_mlmin]
    if len(Q_unsafe) > 0:
        fig.add_trace(go.Scatter(
            x=np.concatenate([[Q_max_mlmin], Q_unsafe, [Q_unsafe[-1]], [Q_max_mlmin]]),
            y=np.concatenate([[dP_max_bar], dP_unsafe, [0], [0]]),
            fill="toself",
            fillcolor="rgba(220, 50, 50, 0.12)",
            line=dict(color="rgba(0,0,0,0)"),
            name="Unsafe zone",
            hoverinfo="skip",
        ))

    # dP curve
    fig.add_trace(go.Scatter(
        x=Q_mlmin,
        y=dP_bar,
        mode="lines",
        name="\u0394P (Kozeny-Carman)",
        line=dict(color="darkorange", width=2),
    ))

    # Q_max vertical marker
    fig.add_vline(
        x=Q_max_mlmin,
        line_dash="dash",
        line_color="tomato",
        annotation_text=f"Q_max={Q_max_mlmin:.2f} mL/min",
        annotation_position="top left",
    )

    fig.update_layout(
        title=f"Pressure-Flow Curve (E\u22c6={column.E_star/1000:.0f} kPa, dp={column.particle_diameter*1e6:.0f} \u00b5m)",
        xaxis_title="Flow rate (mL/min)",
        yaxis_title="\u0394P (bar)",
        height=400,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig
