"""Unit conversion and display formatting for the EmulSim UI.

Conventions:
  - All simulation internals use SI units (Pa, m, m^3/s, mol/m^3, K).
  - Display functions convert to human-readable units for the UI only.
  - No conversion is applied to the simulation data itself.

Quick reference:
  1 mol/m^3 = 1 mM   (since 1 L = 1e-3 m^3, 1 mol/m^3 = 1 mmol/L = 1 mM)
  1 Pa      = 1e-3 kPa = 1e-5 bar
  1 m       = 1e6 um  = 1e9 nm
  1 m^3/s   = 60e6 mL/min
"""

from __future__ import annotations

from typing import Optional
import math


# ─── Concentration ────────────────────────────────────────────────────────────

def mol_m3_to_mM(val: float) -> float:
    """Convert mol/m^3 to millimolar [mM]. Conversion factor is 1.0."""
    return float(val)  # 1 mol/m^3 == 1 mmol/L == 1 mM


def mM_to_mol_m3(val: float) -> float:
    """Convert millimolar [mM] to mol/m^3."""
    return float(val)


# ─── Pressure ─────────────────────────────────────────────────────────────────

def Pa_to_kPa(val: float) -> float:
    """Convert Pascal to kilopascal."""
    return val / 1_000.0


def Pa_to_bar(val: float) -> float:
    """Convert Pascal to bar (1 bar = 1e5 Pa)."""
    return val / 1.0e5


def bar_to_Pa(val: float) -> float:
    """Convert bar to Pascal."""
    return val * 1.0e5


def kPa_to_Pa(val: float) -> float:
    """Convert kilopascal to Pascal."""
    return val * 1_000.0


# ─── Length ───────────────────────────────────────────────────────────────────

def m_to_um(val: float) -> float:
    """Convert metres to micrometres [um]."""
    return val * 1.0e6


def m_to_nm(val: float) -> float:
    """Convert metres to nanometres [nm]."""
    return val * 1.0e9


def um_to_m(val: float) -> float:
    """Convert micrometres to metres."""
    return val / 1.0e6


def nm_to_m(val: float) -> float:
    """Convert nanometres to metres."""
    return val / 1.0e9


# ─── Flow rate ────────────────────────────────────────────────────────────────

def m3_s_to_mL_min(val: float) -> float:
    """Convert m^3/s to mL/min (1 m^3/s = 60e6 mL/min)."""
    return val * 60.0e6


def mL_min_to_m3_s(val: float) -> float:
    """Convert mL/min to m^3/s."""
    return val / 60.0e6


# ─── Temperature ──────────────────────────────────────────────────────────────

def K_to_C(val: float) -> float:
    """Convert Kelvin to Celsius."""
    return val - 273.15


def C_to_K(val: float) -> float:
    """Convert Celsius to Kelvin."""
    return val + 273.15


# ─── Formatting helpers ───────────────────────────────────────────────────────

def format_sci(val: float, sig: int = 3) -> str:
    """Format a number in scientific notation with `sig` significant figures.

    Examples:
        format_sci(0.000123, 3)  -> '1.23e-04'
        format_sci(1234567, 3)   -> '1.23e+06'
        format_sci(0.0, 3)       -> '0.000e+00'
    """
    if not math.isfinite(val):
        return str(val)
    if val == 0.0:
        fmt = f"0.{'0' * (sig - 1)}e+00"
        return fmt
    # Python's built-in e-notation
    precision = max(sig - 1, 0)
    return f"{val:.{precision}e}"


def format_with_unit(val: float, unit: str, sig: int = 3) -> str:
    """Format a value with its display unit, choosing sci or fixed notation.

    Uses scientific notation for values outside [0.01, 9999]; fixed otherwise.

    Args:
        val:  Numeric value (already in display units).
        unit: Unit string appended after the number, e.g. 'mM', 'bar', 'nm'.
        sig:  Number of significant figures.

    Returns:
        Formatted string like '12.3 mM' or '1.23e-04 m^3/s'.

    Examples:
        format_with_unit(12.345, 'mM', 3)      -> '12.3 mM'
        format_with_unit(0.000123, 'm', 3)      -> '1.23e-04 m'
        format_with_unit(3.0e5, 'Pa', 4)        -> '3.000e+05 Pa'
    """
    if not math.isfinite(val):
        return f"{val} {unit}"
    abs_val = abs(val)
    if abs_val == 0.0 or (0.01 <= abs_val < 10_000.0):
        # Fixed notation: choose decimal places to hit sig figs
        if abs_val == 0.0:
            decimals = sig - 1
        else:
            magnitude = math.floor(math.log10(abs_val))
            decimals = max(sig - 1 - magnitude, 0)
        return f"{val:.{decimals}f} {unit}"
    else:
        return f"{format_sci(val, sig)} {unit}"


def format_pressure(Pa: float) -> str:
    """Format a pressure value, choosing bar for >= 100 Pa, kPa otherwise."""
    if abs(Pa) >= 1.0e5:
        return format_with_unit(Pa_to_bar(Pa), "bar")
    elif abs(Pa) >= 100.0:
        return format_with_unit(Pa_to_kPa(Pa), "kPa")
    else:
        return format_with_unit(Pa, "Pa")


def format_length(m: float) -> str:
    """Format a length in metres, auto-selecting um or nm display range."""
    if abs(m) >= 1.0e-4:
        return format_with_unit(m * 1e2, "cm")
    elif abs(m) >= 1.0e-6:
        return format_with_unit(m_to_um(m), "um")
    else:
        return format_with_unit(m_to_nm(m), "nm")


def format_concentration(mol_m3: float) -> str:
    """Format a concentration in mol/m^3 as mM."""
    return format_with_unit(mol_m3_to_mM(mol_m3), "mM")


def format_flow_rate(m3_s: float) -> str:
    """Format a volumetric flow rate in mL/min."""
    return format_with_unit(m3_s_to_mL_min(m3_s), "mL/min")
