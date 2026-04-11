"""Input and result validation for the EmulSim UI.

Returns ValidationResult objects containing:
  - valid:    False only when hard blockers exist (simulation should NOT run).
  - blockers: List of blocking error messages that prevent meaningful results.
  - warnings: List of advisory messages that degrade confidence but allow running.

26 rules total: 9 M1 + 11 M2 + 6 M3.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# ─── Result container ─────────────────────────────────────────────────────────

@dataclass
class ValidationResult:
    """Container for validation outcome.

    Attributes:
        valid:    True iff no blockers were raised. Warnings do not affect validity.
        blockers: Hard error messages; simulation output should not be trusted.
        warnings: Advisory messages; simulation can proceed but with reduced confidence.
    """
    valid: bool = True
    blockers: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def merge(self, other: "ValidationResult") -> "ValidationResult":
        """Return a new ValidationResult combining self with another."""
        return ValidationResult(
            valid=self.valid and other.valid,
            blockers=self.blockers + other.blockers,
            warnings=self.warnings + other.warnings,
        )

    def add_blocker(self, msg: str) -> None:
        self.blockers.append(msg)
        self.valid = False

    def add_warning(self, msg: str) -> None:
        self.warnings.append(msg)


# ─── M1 Validators (9 rules) ──────────────────────────────────────────────────

# Step types supported by M2.
_M2_SUPPORTED_STEP_TYPES = {
    "secondary_crosslinking", "activation",
    "ligand_coupling", "protein_coupling", "quenching",
}

# Temperature above which standard agarose gel re-melts [Celsius].
_T_GEL_REMELT_LIMIT_C = 85.0

# Typical T_gel for standard agarose [Celsius] — used for T_crosslink check.
_T_AGAROSE_GEL_DEFAULT_C = 37.0


def validate_m1_inputs(
    rpm: float,
    phi_d: float,
    c_agarose: float,
    c_chitosan: float,
    dda: float,
    crosslinker_key: str,
    crosslinker_conc: float,
    T_crosslink: float,
    T_oil: float,
) -> ValidationResult:
    """Validate M1 (emulsification + gelation + crosslinking) inputs.

    9 rules:

    R1.  phi_d > 0.50              -> BLOCKER   (emulsion inversion, model invalid)
    R2.  phi_d > 0.30              -> WARNING   (concentrated regime, coalescence risk)
    R3.  phi_d <= 0.0              -> BLOCKER   (non-physical)
    R4.  c_agarose < 2 or > 6     -> WARNING   (outside empirical pore model range)
    R5.  c_agarose <= 0            -> BLOCKER   (non-physical)
    R6.  c_chitosan < 1 or > 3    -> WARNING   (outside validated range)
    R7.  c_chitosan <= 0           -> BLOCKER   (non-physical, no amine sites)
    R8.  T_crosslink >= T_gel_remelt -> BLOCKER (agarose re-melts, microspheres dissolve)
    R9.  dda < 0.5 or dda > 1.0   -> WARNING   (DDA outside typical range)

    Args:
        rpm:             Mixing speed [RPM].
        phi_d:           Dispersed-phase volume fraction [-].
        c_agarose:       Agarose concentration [wt% as g/100 mL, i.e. % w/v].
        c_chitosan:      Chitosan concentration [mg/mL].
        dda:             Degree of deacetylation [-], dimensionless in [0, 1].
        crosslinker_key: Key into CROSSLINKERS registry.
        crosslinker_conc: Crosslinker concentration [mol/m^3].
        T_crosslink:     Crosslinking temperature [Celsius].
        T_oil:           Oil-phase temperature [Celsius].

    Returns:
        ValidationResult with blockers and warnings populated.
    """
    result = ValidationResult()

    # R1 + R3: phi_d range
    if phi_d <= 0.0:
        result.add_blocker(
            f"Dispersed-phase volume fraction phi_d={phi_d:.3f} is non-physical (must be > 0)."
        )
    elif phi_d > 0.50:
        result.add_blocker(
            f"phi_d={phi_d:.2f} exceeds 0.50 — emulsion inversion regime. "
            "Breakage model is not valid above phi_d = 0.50."
        )
    elif phi_d > 0.30:
        result.add_warning(
            f"phi_d={phi_d:.2f} is in the concentrated regime (> 0.30). "
            "Coalescence risk increases significantly; results may underestimate d32."
        )

    # R5 + R4: c_agarose range
    if c_agarose <= 0.0:
        result.add_blocker(
            f"Agarose concentration c_agarose={c_agarose:.2f} is non-physical (must be > 0)."
        )
    elif c_agarose < 2.0 or c_agarose > 6.0:
        result.add_warning(
            f"c_agarose={c_agarose:.1f}% is outside empirical pore model calibration range "
            "(2–6%). Pore size estimates may be unreliable."
        )

    # R7 + R6: c_chitosan range
    if c_chitosan <= 0.0:
        result.add_blocker(
            f"Chitosan concentration c_chitosan={c_chitosan:.2f} mg/mL is non-physical. "
            "No amine sites available for crosslinking or functionalization."
        )
    elif c_chitosan < 1.0 or c_chitosan > 3.0:
        result.add_warning(
            f"c_chitosan={c_chitosan:.1f} mg/mL is outside validated range (1–3 mg/mL). "
            "Crosslinking kinetics may be extrapolated."
        )

    # R8: T_crosslink vs gel re-melt
    if T_crosslink >= _T_GEL_REMELT_LIMIT_C:
        result.add_blocker(
            f"T_crosslink={T_crosslink:.0f} C is at or above the agarose re-melt temperature "
            f"({_T_GEL_REMELT_LIMIT_C:.0f} C). Microspheres will dissolve. "
            "Lower T_crosslink below 85 C."
        )
    elif T_crosslink < 0.0:
        result.add_blocker(
            f"T_crosslink={T_crosslink:.0f} C is below 0 C — non-physical."
        )

    # R9: DDA range
    if dda < 0.5 or dda > 1.0:
        result.add_warning(
            f"DDA={dda:.2f} is outside the typical range (0.50–1.00). "
            "Amine site density calculations may be unreliable."
        )

    # crosslinker_conc sanity
    if crosslinker_conc < 0.0:
        result.add_blocker(
            f"Crosslinker concentration={crosslinker_conc:.3f} mol/m^3 is negative — non-physical."
        )

    return result


# ─── M2 Validators (11 rules) ─────────────────────────────────────────────────

def validate_m2_inputs(
    steps: list,
    acs_state,
    m1_trust_level: str,
) -> ValidationResult:
    """Validate M2 (functionalization) inputs.

    11 rules:

    R1.  m1_trust_level == 'UNRELIABLE' -> BLOCKER (cascade gate)
    R2.  steps is empty                 -> BLOCKER (no workflow defined)
    R3.  Any step_type not in SUPPORTED -> BLOCKER (unimplemented chemistry)
    R4.  LIGAND_COUPLING detected        -> BLOCKER with explicit 'Planned' message
    R5.  PROTEIN_COUPLING detected       -> BLOCKER with explicit 'Planned' message
    R6.  QUENCHING detected              -> BLOCKER with explicit 'Planned' message
    R7.  acs_state is None               -> BLOCKER (no surface sites to modify)
    R8.  Total accessible sites <= 0     -> BLOCKER (nothing to react)
    R9.  More than 6 steps               -> WARNING (combinatorial complexity)
    R10. Duplicate step types in sequence -> WARNING (may indicate copy-paste error)
    R11. acs_state.validate() returns violations -> WARNING (conservation violated)

    Args:
        steps:          List of modification step objects (must have .step_type attribute
                        with a .value string, or bare strings).
        acs_state:      ACSProfile object from M1/surface area calculation, or None.
        m1_trust_level: String trust level from M1 TrustAssessment ('TRUSTWORTHY',
                        'CAUTION', or 'UNRELIABLE').

    Returns:
        ValidationResult with blockers and warnings populated.
    """
    result = ValidationResult()

    # R1: Upstream trust gate — must check first
    if m1_trust_level == "UNRELIABLE":
        result.add_blocker(
            "M1 trust level is UNRELIABLE. M2 simulation is blocked until M1 inputs "
            "are corrected. Fix M1 blockers before proceeding."
        )

    # R7: ACS state availability
    if acs_state is None:
        result.add_blocker(
            "No ACS surface-site inventory available from M1. "
            "Run M1 simulation first to generate ACSProfile."
        )

    # R2: Empty workflow
    if not steps:
        result.add_blocker(
            "No modification steps defined. Add at least one step to define a workflow."
        )
        return result  # Cannot check step types if list is empty

    # Extract step type strings
    step_type_values = []
    for s in steps:
        if hasattr(s, "step_type"):
            val = s.step_type.value if hasattr(s.step_type, "value") else str(s.step_type)
        else:
            val = str(s)
        step_type_values.append(val)

    # R3–R6: Check each step type is supported
    seen_types = []
    for val in step_type_values:
        if val not in _M2_SUPPORTED_STEP_TYPES:
            result.add_blocker(
                f"Unknown step type '{val}'. "
                f"Supported types: {sorted(_M2_SUPPORTED_STEP_TYPES)}."
            )
        seen_types.append(val)

    # R9: Too many steps
    if len(steps) > 6:
        result.add_warning(
            f"{len(steps)} modification steps defined — this is a large workflow. "
            "Verify step ordering and that no steps are duplicated."
        )

    # R10: Duplicate step types
    seen = set()
    duplicates = set()
    for val in seen_types:
        if val in seen:
            duplicates.add(val)
        seen.add(val)
    if duplicates:
        result.add_warning(
            f"Duplicate step types detected: {sorted(duplicates)}. "
            "Verify this is intentional (e.g., sequential crosslinking passes)."
        )

    # R8 + R11: ACS state quality (only if acs_state is provided)
    if acs_state is not None:
        # R8: accessible sites
        accessible = getattr(acs_state, "accessible", None)
        if accessible is not None and accessible <= 0.0:
            result.add_blocker(
                "ACS accessible site count is zero or negative. "
                "No surface sites available for modification."
            )

        # R11: conservation check
        if hasattr(acs_state, "validate"):
            violations = acs_state.validate()
            if violations:
                for v in violations:
                    result.add_warning(f"ACS conservation violation: {v}")

    return result


# ─── M3 Validators (6 rules) ──────────────────────────────────────────────────

def validate_m3_chromatography(
    flow_rate: float,
    column,
    isotherm_type: str,
    gradient_enabled: bool,
) -> ValidationResult:
    """Validate M3 (chromatography) pre-simulation inputs.

    6 rules:

    R1.  flow_rate <= 0                          -> BLOCKER
    R2.  pressure_drop > max_safe_pressure        -> BLOCKER (structural risk)
    R3.  Re > 1                                   -> WARNING (inertial flow)
    R4.  Pe < 5                                   -> WARNING (dispersive regime)
    R5.  gradient + competitive_langmuir          -> WARNING ('diagnostic only')
    R6.  isotherm_type == 'linear' + gradient     -> WARNING (gradient has no effect)

    Args:
        flow_rate:       Volumetric flow rate [m^3/s].
        column:          ColumnGeometry object (must have radius, length, porosity attributes).
        isotherm_type:   String isotherm identifier ('langmuir', 'competitive_langmuir',
                         'linear', etc.).
        gradient_enabled: True if a gradient elution program is active.

    Returns:
        ValidationResult with blockers and warnings populated.
    """
    result = ValidationResult()

    # R1: flow rate must be positive
    if flow_rate <= 0.0:
        result.add_blocker(
            f"Flow rate={flow_rate:.3e} m^3/s is non-positive — non-physical."
        )
        return result  # Cannot compute Re/Pe without valid flow rate

    # Compute hydrodynamic numbers if column geometry is available
    if column is not None:
        import math

        radius = getattr(column, "radius", None)
        length = getattr(column, "length", None)
        porosity = getattr(column, "porosity", None)
        particle_diameter = getattr(column, "particle_diameter", None)

        if radius and length and porosity and particle_diameter:
            # Cross-sectional area and superficial velocity
            area = math.pi * radius ** 2
            u_sup = flow_rate / area  # [m/s]

            # R3: Reynolds number (using water viscosity ~1e-3 Pa.s, density ~1000 kg/m^3)
            mu_water = 1.0e-3   # Pa.s
            rho_water = 1000.0  # kg/m^3
            d_p = particle_diameter
            Re = rho_water * u_sup * d_p / mu_water
            if Re > 1.0:
                result.add_warning(
                    f"Reynolds number Re={Re:.2f} > 1 — inertial flow regime. "
                    "Kozeny-Carman pressure model assumes creeping flow (Re << 1). "
                    "Pressure drop estimate may be underestimated."
                )

            # R4: Peclet number (axial dispersion; D_axial ~ D_mol * (1 + Pe*Re*Sc/30))
            # Rough Pe = u_sup * d_p / D_mol, using D_mol ~ 1e-10 m^2/s for proteins
            D_mol = 1.0e-10  # m^2/s, typical protein diffusivity
            Pe = u_sup * d_p / D_mol
            if Pe < 5.0:
                result.add_warning(
                    f"Column Peclet number Pe={Pe:.1f} < 5 — dispersive regime. "
                    "Axial dispersion is significant; LRM plug-flow assumption is less accurate."
                )

            # R2: Pressure drop estimate (Kozeny-Carman)
            #     dP = 180 * mu * u_sup * L * (1-eps)^2 / (d_p^2 * eps^3)
            eps = porosity
            dP_estimate = (
                180.0 * mu_water * u_sup * length
                * (1.0 - eps) ** 2
                / (d_p ** 2 * eps ** 3)
            )
            max_safe_pressure = getattr(column, "max_safe_pressure", 3.0e5)  # 3 bar default
            if dP_estimate > max_safe_pressure:
                result.add_blocker(
                    f"Estimated pressure drop {dP_estimate/1e5:.2f} bar exceeds safe limit "
                    f"{max_safe_pressure/1e5:.2f} bar. Reduce flow rate or increase column diameter."
                )

    # R5: Gradient + competitive Langmuir — diagnostic only
    if gradient_enabled and isotherm_type in ("competitive_langmuir",):
        result.add_warning(
            "Gradient elution with competitive Langmuir isotherm is in diagnostic mode. "
            "Gradient salt concentration does not yet affect binding affinity (BF-2 not deployed). "
            "Results show elution order only — quantitative yields are not reliable."
        )

    # R6: Gradient + linear isotherm — gradient has no mechanistic effect
    if gradient_enabled and isotherm_type == "linear":
        result.add_warning(
            "Linear isotherm does not depend on mobile-phase composition. "
            "Gradient elution will have no effect on retention — use Langmuir or SMA isotherm."
        )

    return result


def validate_m3_result(
    mass_balance_error: float,
    pressure_drop: float,
    max_safe_pressure: float = 3.0e5,
) -> ValidationResult:
    """Post-simulation validation of M3 result quality.

    Rules:
        R1. mass_balance_error > 5%   -> BLOCKER  (numerical result is not trustworthy)
        R2. mass_balance_error 2-5%   -> WARNING  (treat with caution)
        R3. mass_balance_error < 2%   -> OK       (acceptable)
        R4. pressure_drop > max_safe  -> BLOCKER  (exceeded safe operating limit)

    Args:
        mass_balance_error: Relative mass balance error (dimensionless, e.g. 0.03 = 3%).
        pressure_drop:      Actual computed pressure drop [Pa].
        max_safe_pressure:  Safe pressure limit [Pa]; default 3 bar.

    Returns:
        ValidationResult with blockers and warnings populated.
    """
    result = ValidationResult()

    # R1–R3: Mass balance gating
    mb_pct = abs(mass_balance_error) * 100.0
    if mb_pct > 5.0:
        result.add_blocker(
            f"Mass balance error = {mb_pct:.1f}% exceeds 5% threshold. "
            "Simulation results are numerically unreliable. "
            "Reduce time step or increase spatial resolution."
        )
    elif mb_pct > 2.0:
        result.add_warning(
            f"Mass balance error = {mb_pct:.1f}% (2–5% range). "
            "Results should be treated with caution. "
            "Consider tightening solver tolerances."
        )
    # else: < 2% is acceptable, no message needed

    # R4: Pressure safety post-simulation
    if pressure_drop > max_safe_pressure:
        result.add_blocker(
            f"Computed pressure drop {pressure_drop/1e5:.2f} bar exceeds safe limit "
            f"{max_safe_pressure/1e5:.2f} bar. Column bed integrity may be compromised."
        )

    return result
