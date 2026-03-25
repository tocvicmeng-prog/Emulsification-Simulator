"""Visualization utilities for simulation results.

Launch the web UI:  streamlit run src/emulsim/visualization/app.py
"""

from .plots import (
    plot_droplet_size_distribution,
    plot_phase_field,
    plot_crosslinking_kinetics,
    plot_hertz_contact,
    plot_kav_curve,
    plot_results_dashboard,
)
