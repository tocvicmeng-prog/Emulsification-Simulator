"""Output-level metadata for scientific UI trust labeling.

Every numerical output displayed in the UI must carry one of four trust tiers:
  - MECHANISTIC:           validated governing equations + calibrated parameters
  - EMPIRICAL_CALIBRATED:  fitted correlation within calibration range (+/-20-30%)
  - SEMI_QUANTITATIVE:     uncalibrated defaults, relative ordering only
  - RANKING_ONLY:          directional trends only, no quantitative meaning
  - NOT_PREDICTED:         placeholder, model not implemented
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ModelBasis(Enum):
    MECHANISTIC = "mechanistic"
    EMPIRICAL_CALIBRATED = "empirical_calibrated"
    SEMI_QUANTITATIVE = "semi_quantitative"
    RANKING_ONLY = "ranking_only"
    NOT_PREDICTED = "not_predicted"


class ConfidenceLevel(Enum):
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    NONE = "none"


_ICONS = {
    "high":     "[OK]",
    "moderate": "[~]",
    "low":      "[!]",
    "none":     "[X]",
}


@dataclass(frozen=True)
class OutputMetadata:
    """Immutable trust descriptor attached to every simulated output.

    Attributes:
        model_basis:          Scientific basis of the model producing this output.
        confidence:           Overall confidence tier (high/moderate/low/none).
        validity_range:       Human-readable description of calibration range.
        calibration_required: True if user calibration is needed before trusting numbers.
        warnings:             Tuple of advisory strings to surface in the UI.
        recommended_use:      Short description of how/when this output should be used.
        source_module:        Which module produces this output (e.g., "M1", "M2", "M3").
    """

    model_basis: ModelBasis
    confidence: ConfidenceLevel
    validity_range: str = ""
    calibration_required: bool = False
    warnings: tuple = ()
    recommended_use: str = ""
    source_module: str = ""

    @property
    def icon(self) -> str:
        """Short text icon suitable for Streamlit badges."""
        return _ICONS[self.confidence.value]

    @property
    def label(self) -> str:
        """One-line label combining icon and basis, e.g. '[~] empirical_calibrated'."""
        return f"{self.icon} {self.model_basis.value}"


# ─── Pre-built metadata for every major output ────────────────────────────────

# M1: Emulsification outputs
M1_D50_META = OutputMetadata(
    model_basis=ModelBasis.MECHANISTIC,
    confidence=ConfidenceLevel.MODERATE,
    validity_range="d32 = 0.5–100 um; phi_d < 0.30; RPM 100–3000",
    calibration_required=False,
    warnings=("Kolmogorov-scale breakage model; not validated below 0.5 um.",),
    recommended_use="Formulation comparison and trend analysis.",
    source_module="M1",
)

M1_SPAN_META = OutputMetadata(
    model_basis=ModelBasis.MECHANISTIC,
    confidence=ConfidenceLevel.MODERATE,
    validity_range="phi_d < 0.30",
    calibration_required=False,
    warnings=("Span > 2 indicates broad distribution not suitable for chromatography.",),
    recommended_use="Polydispersity screening.",
    source_module="M1",
)

M1_PORE_META = OutputMetadata(
    model_basis=ModelBasis.EMPIRICAL_CALIBRATED,
    confidence=ConfidenceLevel.LOW,
    validity_range="Agarose 1–8 wt%; pore size 20–1000 nm",
    calibration_required=True,
    warnings=(
        "Empirical pore correlation; calibrated against literature agarose data only.",
        "Outside 1–8% agarose range, extrapolation may be unreliable.",
    ),
    recommended_use="Relative comparison of pore size across formulations.",
    source_module="M1",
)

M1_GELATION_META = OutputMetadata(
    model_basis=ModelBasis.MECHANISTIC,
    confidence=ConfidenceLevel.MODERATE,
    validity_range="Cooling rate 0.1–10 K/min; T_gel from literature",
    calibration_required=False,
    warnings=("Gelation extent < 80% indicates microspheres may not form properly.",),
    recommended_use="Verify process conditions for complete gelation.",
    source_module="M1",
)

M1_CROSSLINK_META = OutputMetadata(
    model_basis=ModelBasis.MECHANISTIC,
    confidence=ConfidenceLevel.MODERATE,
    validity_range="Genipin: c = 0.1–10 mol/m3; T = 25–60 C",
    calibration_required=True,
    warnings=(
        "Crosslinking kinetics calibrated for genipin only.",
        "Other crosslinkers use simplified rate models.",
    ),
    recommended_use="Crosslinking extent and time optimization.",
    source_module="M1",
)

M1_MECHANICAL_META = OutputMetadata(
    model_basis=ModelBasis.SEMI_QUANTITATIVE,
    confidence=ConfidenceLevel.LOW,
    validity_range="Rubber-elasticity IPN model; phenomenological coupling",
    calibration_required=True,
    warnings=(
        "G_DN uses phenomenological coupling formula (G1 + G2 + eta*sqrt(G1*G2)).",
        "Suitable for formulation ranking, not absolute mechanical prediction.",
    ),
    recommended_use="Ranking formulations by stiffness, not absolute modulus values.",
    source_module="M1",
)

# M2: Functionalization outputs
M2_ACS_META = OutputMetadata(
    model_basis=ModelBasis.SEMI_QUANTITATIVE,
    confidence=ConfidenceLevel.LOW,
    validity_range=(
        "All 9 backend step types: SECONDARY_CROSSLINKING, ACTIVATION, "
        "LIGAND_COUPLING, PROTEIN_COUPLING, QUENCHING, SPACER_ARM, "
        "METAL_CHARGING, PROTEIN_PRETREATMENT, WASHING"
    ),
    calibration_required=True,
    warnings=(
        "ACS inventory uses simplified site-density model.",
        "Protein coupling outputs are ranking_only unless calibrated.",
        "pH and temperature effects use Arrhenius approximation only.",
    ),
    recommended_use="Relative comparison of surface site availability between formulations.",
    source_module="M2",
)

M2_STEP_CONVERSION_META = OutputMetadata(
    model_basis=ModelBasis.EMPIRICAL_CALIBRATED,
    confidence=ConfidenceLevel.LOW,
    validity_range="Second-order kinetics; validated for genipin and ECH only",
    calibration_required=True,
    warnings=(
        "Conversion fraction uses second-order consumption model.",
        "Default rate constants are illustrative — user calibration required.",
    ),
    recommended_use="Trend analysis; not for absolute yield prediction.",
    source_module="M2",
)

# M3: Chromatography outputs
M3_BREAKTHROUGH_META = OutputMetadata(
    model_basis=ModelBasis.MECHANISTIC,
    confidence=ConfidenceLevel.MODERATE,
    validity_range="Single-component Langmuir LRM; Pe > 10; Re < 1",
    calibration_required=True,
    warnings=(
        "Default isotherm parameters are illustrative — user calibration required.",
        "LRM assumes instantaneous local equilibrium between mobile and stationary phases.",
    ),
    recommended_use="DBC estimation and loading strategy optimization.",
    source_module="M3",
)

M3_DBC_META = OutputMetadata(
    model_basis=ModelBasis.MECHANISTIC,
    confidence=ConfidenceLevel.MODERATE,
    validity_range="DBC at 5/10/50% breakthrough; single-component only",
    calibration_required=True,
    warnings=(
        "Multi-component competition effects are not captured in single-component mode.",
    ),
    recommended_use="Column loading capacity design.",
    source_module="M3",
)

M3_PRESSURE_META = OutputMetadata(
    model_basis=ModelBasis.EMPIRICAL_CALIBRATED,
    confidence=ConfidenceLevel.MODERATE,
    validity_range="Kozeny-Carman; Re < 1 (creeping flow)",
    calibration_required=False,
    warnings=(
        "Assumes rigid, incompressible bed — soft gels will compress under load.",
        "Pressure above 3 bar may damage agarose microspheres.",
    ),
    recommended_use="Pressure safety screening before column packing.",
    source_module="M3",
)

M3_GRADIENT_META = OutputMetadata(
    model_basis=ModelBasis.SEMI_QUANTITATIVE,
    confidence=ConfidenceLevel.LOW,
    validity_range="Gradient-sensitive SMA, HIC, IMAC, Protein A, and lectin adapters",
    calibration_required=True,
    warnings=(
        "Gradient-sensitive isotherms update binding during elution; plain competitive "
        "Langmuir remains diagnostic only.",
        "Default isotherm parameters are illustrative — user calibration required.",
    ),
    recommended_use="Qualitative elution order prediction; not for yield optimization.",
    source_module="M3",
)

M3_MASS_BALANCE_META = OutputMetadata(
    model_basis=ModelBasis.MECHANISTIC,
    confidence=ConfidenceLevel.HIGH,
    validity_range="Mass balance error computed from numerical integration",
    calibration_required=False,
    warnings=(
        "Error > 5% is a BLOCKER — results are unreliable.",
        "Error 2–5% is a WARNING — treat results with caution.",
    ),
    recommended_use="Quality gate: must pass before reporting other M3 outputs.",
    source_module="M3",
)
