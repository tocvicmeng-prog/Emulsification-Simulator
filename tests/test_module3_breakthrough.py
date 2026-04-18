"""Tests for Module 3: Chromatographic Performance (Phase C).

Covers:
  C.1  Kozeny-Carman pressure drop + compressibility
  C.2  Langmuir isotherm equilibrium + jacobian
  C.3  LRM transport: mass balance, breakthrough shape
  C.4  UV detection: Beer-Lambert
  C.5  Orchestrator: run_breakthrough completeness
"""

from __future__ import annotations


import numpy as np
import pytest

from emulsim.module3_performance.hydrodynamics import ColumnGeometry
from emulsim.module3_performance.isotherms.langmuir import LangmuirIsotherm
from emulsim.module3_performance.transport.lumped_rate import (
    solve_lrm,
    LRMResult,
    classify_mass_balance,
    MassBalanceQuality,
)
from emulsim.module3_performance.detection.uv import (
    compute_uv_signal,
    apply_detector_broadening,
)
from emulsim.module3_performance.orchestrator import (
    BreakthroughResult,
    run_breakthrough,
)


# ─── Fixtures ─────────────────────────────────────────────────────────

@pytest.fixture
def default_column() -> ColumnGeometry:
    return ColumnGeometry()


@pytest.fixture
def default_isotherm() -> LangmuirIsotherm:
    return LangmuirIsotherm()


# ─── C.1: Hydrodynamics ──────────────────────────────────────────────

class TestPressureDrop:
    """Kozeny-Carman pressure drop tests."""

    def test_pressure_drop_kozeny_carman(self, default_column: ColumnGeometry):
        """dP matches hand calculation for known geometry.

        For default column:
            dp = 100e-6 m, eps = 0.38, L = 0.10 m, mu = 1e-3 Pa.s
            A = pi/4 * 0.01^2 = 7.854e-5 m^2
            Q = 1e-8 m^3/s => u = 1.274e-4 m/s

        dP = 150 * 1e-3 * 1.274e-4 * 0.10 * (0.62)^2 / (1e-8 * 0.38^3)
           = 150 * 1e-3 * 1.274e-4 * 0.10 * 0.3844 / (1e-8 * 0.054872)
        """
        Q = 1e-8  # m^3/s
        dP = default_column.pressure_drop(Q)

        # Hand calculation
        u = Q / default_column.cross_section_area
        eps = 0.38
        dp = 100e-6
        L = 0.10
        mu = 1e-3
        dP_expected = 150 * mu * u * L * (1 - eps) ** 2 / (dp ** 2 * eps ** 3)

        assert dP == pytest.approx(dP_expected, rel=1e-10)
        assert dP > 0

    def test_pressure_drop_scales_with_flow(self, default_column: ColumnGeometry):
        """dP is linear in flow rate (Kozeny-Carman)."""
        Q1 = 1e-8
        Q2 = 2e-8
        dP1 = default_column.pressure_drop(Q1)
        dP2 = default_column.pressure_drop(Q2)
        assert dP2 == pytest.approx(2 * dP1, rel=1e-10)

    def test_pressure_drop_zero_flow(self, default_column: ColumnGeometry):
        """Zero flow gives zero pressure drop."""
        assert default_column.pressure_drop(0.0) == 0.0


class TestMaxSafeFlowRate:
    """Max safe flow rate tests."""

    def test_max_safe_flow_rate_positive(self, default_column: ColumnGeometry):
        """Q_max is positive and finite."""
        Q_max = default_column.max_safe_flow_rate()
        assert Q_max > 0
        assert np.isfinite(Q_max)

    def test_max_safe_inversely_proportional_to_mu(self):
        """u_max is inversely proportional to mu."""
        col = ColumnGeometry()
        Q1 = col.max_safe_flow_rate(mu=1e-3)
        Q2 = col.max_safe_flow_rate(mu=2e-3)
        assert Q1 == pytest.approx(2 * Q2, rel=1e-10)

    def test_max_safe_inversely_proportional_to_L(self):
        """u_max is inversely proportional to bed height."""
        col1 = ColumnGeometry(bed_height=0.10)
        col2 = ColumnGeometry(bed_height=0.20)
        Q1 = col1.max_safe_flow_rate()
        Q2 = col2.max_safe_flow_rate()
        assert Q1 == pytest.approx(2 * Q2, rel=1e-10)

    def test_pressure_at_max_flow_equals_safety_E_star(self, default_column: ColumnGeometry):
        """At Q_max, dP = safety * E_star."""
        safety = 0.8
        Q_max = default_column.max_safe_flow_rate(safety=safety)
        dP = default_column.pressure_drop(Q_max)
        assert dP == pytest.approx(safety * default_column.E_star, rel=1e-6)


class TestBedCompression:
    """Bed compression tests."""

    def test_compression_linear(self, default_column: ColumnGeometry):
        """Compression is linear in dP."""
        c1 = default_column.bed_compression_fraction(1000.0)
        c2 = default_column.bed_compression_fraction(2000.0)
        assert c2 == pytest.approx(2 * c1, rel=1e-10)

    def test_zero_pressure_zero_compression(self, default_column: ColumnGeometry):
        assert default_column.bed_compression_fraction(0.0) == 0.0


# ─── C.2: Langmuir Isotherm ──────────────────────────────────────────

class TestLangmuirEquilibrium:
    """Langmuir isotherm equilibrium loading tests."""

    def test_equilibrium_at_high_C(self, default_isotherm: LangmuirIsotherm):
        """q* approaches q_max at C >> K_d."""
        C_high = 1000.0 / default_isotherm.K_L  # 1000 * K_d
        q = default_isotherm.equilibrium_loading(C_high)
        # q*/q_max = K_L*C/(1+K_L*C) = 1000/1001 ~ 0.999
        assert q == pytest.approx(default_isotherm.q_max, rel=0.002)

    def test_equilibrium_at_zero(self, default_isotherm: LangmuirIsotherm):
        """q*(0) = 0."""
        assert default_isotherm.equilibrium_loading(0.0) == 0.0

    def test_equilibrium_at_Kd(self, default_isotherm: LangmuirIsotherm):
        """q*(K_d) = q_max / 2."""
        K_d = default_isotherm.K_d
        q = default_isotherm.equilibrium_loading(K_d)
        assert q == pytest.approx(default_isotherm.q_max / 2.0, rel=1e-10)

    def test_equilibrium_array(self, default_isotherm: LangmuirIsotherm):
        """Works with numpy arrays."""
        C = np.array([0.0, 0.001, 0.01, 0.1, 1.0])
        q = default_isotherm.equilibrium_loading(C)
        assert q.shape == (5,)
        assert np.all(q >= 0)
        # Monotonically increasing
        assert np.all(np.diff(q) >= 0)

    def test_negative_C_clipped(self, default_isotherm: LangmuirIsotherm):
        """Negative concentrations give q* = 0."""
        assert default_isotherm.equilibrium_loading(-1.0) == 0.0


class TestLangmuirJacobian:
    """Langmuir isotherm Jacobian (dq*/dC) tests."""

    def test_jacobian_numerical_vs_analytical(self, default_isotherm: LangmuirIsotherm):
        """Numerical derivative matches analytical Jacobian."""
        C_vals = np.array([0.0001, 0.001, 0.01, 0.1])
        dC = 1e-8
        for C in C_vals:
            analytical = default_isotherm.jacobian(C)
            numerical = (
                default_isotherm.equilibrium_loading(C + dC)
                - default_isotherm.equilibrium_loading(C - dC)
            ) / (2.0 * dC)
            assert analytical == pytest.approx(float(numerical), rel=1e-4), (
                f"Jacobian mismatch at C={C}"
            )

    def test_jacobian_at_zero(self, default_isotherm: LangmuirIsotherm):
        """dq*/dC at C=0 = q_max * K_L (initial slope)."""
        j0 = default_isotherm.jacobian(0.0)
        assert j0 == pytest.approx(
            default_isotherm.q_max * default_isotherm.K_L, rel=1e-10
        )

    def test_jacobian_positive(self, default_isotherm: LangmuirIsotherm):
        """Jacobian is always positive for positive C."""
        C = np.linspace(0, 1.0, 100)
        j = default_isotherm.jacobian(C)
        assert np.all(j > 0)


# ─── C.3: LRM Transport ──────────────────────────────────────────────

class TestMassBalanceQuality:
    """Mass balance quality classification tests."""

    def test_mass_balance_quality_classification(self):
        """Verify quality thresholds map to correct enum members."""
        assert classify_mass_balance(0.01) == MassBalanceQuality.ACCEPTABLE
        assert classify_mass_balance(0.02) == MassBalanceQuality.ACCEPTABLE   # boundary: <= 2% is ACCEPTABLE
        assert classify_mass_balance(0.03) == MassBalanceQuality.CAUTION
        assert classify_mass_balance(0.05) == MassBalanceQuality.CAUTION      # boundary: <= 5% is CAUTION
        assert classify_mass_balance(0.06) == MassBalanceQuality.UNRELIABLE

    def test_lrm_result_has_quality_field(self):
        """LRMResult includes mass_balance_quality as a string."""
        col = ColumnGeometry(bed_height=0.05, particle_diameter=50e-6)
        iso = LangmuirIsotherm(q_max=20.0, K_L=100.0)
        result = solve_lrm(
            column=col,
            isotherm=iso,
            C_feed=1.0,
            feed_duration=300.0,
            flow_rate=1e-8,
            total_time=600.0,
            n_z=60,
        )
        assert hasattr(result, "mass_balance_quality"), (
            "LRMResult is missing the mass_balance_quality field"
        )
        assert isinstance(result.mass_balance_quality, str)
        assert result.mass_balance_quality in {q.value for q in MassBalanceQuality}, (
            f"Unexpected quality value: {result.mass_balance_quality!r}"
        )
        # For a well-resolved run the quality should be acceptable
        assert result.mass_balance_quality == MassBalanceQuality.ACCEPTABLE.value, (
            f"Expected 'acceptable', got {result.mass_balance_quality!r} "
            f"(error={result.mass_balance_error:.4%})"
        )


class TestLRMMassBalance:
    """LRM mass balance tests."""

    def test_lrm_mass_balance(self):
        """Injected = eluted + bound + remaining, error < 2%."""
        # Use moderate binding so breakthrough partially occurs,
        # giving a non-trivial mass balance test with both eluted and bound.
        col = ColumnGeometry(bed_height=0.05, particle_diameter=50e-6)
        iso = LangmuirIsotherm(q_max=20.0, K_L=100.0)
        result = solve_lrm(
            column=col,
            isotherm=iso,
            C_feed=1.0,
            feed_duration=300.0,
            flow_rate=1e-8,
            total_time=600.0,
            n_z=60,
        )
        assert result.mass_balance_error < 0.02, (
            f"Mass balance error {result.mass_balance_error:.4f} exceeds 2%"
        )

    def test_lrm_result_fields(self, default_column: ColumnGeometry, default_isotherm: LangmuirIsotherm):
        """LRMResult has all expected fields with correct shapes."""
        result = solve_lrm(
            column=default_column,
            isotherm=default_isotherm,
            C_feed=1.0,
            feed_duration=100.0,
            flow_rate=1e-8,
            total_time=200.0,
            n_z=20,
        )
        assert isinstance(result, LRMResult)
        assert result.time.ndim == 1
        assert result.z.ndim == 1
        assert len(result.z) == 20
        assert result.C_outlet.shape == result.time.shape
        assert result.q_average.shape == result.time.shape
        assert result.mass_injected > 0
        assert np.isfinite(result.mass_balance_error)


class TestLRMBreakthroughShape:
    """Breakthrough curve shape tests."""

    def test_breakthrough_monotonic_during_loading(self):
        """Breakthrough curve is monotonically increasing during loading."""
        col = ColumnGeometry(bed_height=0.05, particle_diameter=50e-6)
        iso = LangmuirIsotherm(q_max=50.0, K_L=500.0)
        result = solve_lrm(
            column=col,
            isotherm=iso,
            C_feed=1.0,
            feed_duration=600.0,
            flow_rate=1e-8,
            total_time=600.0,
            n_z=30,
        )
        # During loading (after initial transient), outlet should be non-decreasing
        # Skip first ~10% of points (initial zero region)
        n_skip = len(result.time) // 10
        C_loading = result.C_outlet[n_skip:]
        diffs = np.diff(C_loading)
        # Allow tiny numerical noise (1e-10)
        assert np.all(diffs >= -1e-10), (
            f"Breakthrough curve is not monotonic: min diff = {diffs.min():.2e}"
        )

    def test_no_loading_means_pulse(self):
        """Short injection produces a peak that eventually decays.

        Uses a non-binding column (q_max~0) so the pulse travels through
        as a dispersed peak and the outlet returns to near zero.
        """
        col = ColumnGeometry(bed_height=0.05, particle_diameter=50e-6)
        # Very weak binding so pulse passes through
        iso = LangmuirIsotherm(q_max=0.01, K_L=1.0)
        result = solve_lrm(
            column=col,
            isotherm=iso,
            C_feed=1.0,
            feed_duration=5.0,         # very short injection
            flow_rate=1e-8,
            total_time=600.0,
            n_z=30,
        )
        # Outlet should eventually return near zero
        # Check last 10% of time
        n_tail = len(result.C_outlet) // 10
        C_tail = result.C_outlet[-n_tail:]
        C_peak = result.C_outlet.max()
        assert C_peak > 0, "Pulse should produce a peak at outlet"
        assert np.max(C_tail) < 0.1 * C_peak, (
            f"Pulse tail too high — outlet did not return to zero. "
            f"Tail max={np.max(C_tail):.4e}, peak={C_peak:.4e}"
        )

    def test_outlet_initially_zero(self, default_column: ColumnGeometry, default_isotherm: LangmuirIsotherm):
        """Outlet concentration starts near zero."""
        result = solve_lrm(
            column=default_column,
            isotherm=default_isotherm,
            C_feed=1.0,
            feed_duration=300.0,
            flow_rate=1e-8,
            total_time=600.0,
            n_z=40,
        )
        # First few time points should be ~0
        assert result.C_outlet[0] == pytest.approx(0.0, abs=1e-6)

    def test_positivity(self, default_column: ColumnGeometry, default_isotherm: LangmuirIsotherm):
        """All concentrations remain non-negative."""
        result = solve_lrm(
            column=default_column,
            isotherm=default_isotherm,
            C_feed=1.0,
            feed_duration=300.0,
            flow_rate=1e-8,
            total_time=600.0,
            n_z=40,
        )
        assert np.all(result.C_outlet >= 0.0)
        assert np.all(result.q_average >= 0.0)


# ─── C.4: UV Detection ───────────────────────────────────────────────

class TestUVBeerLambert:
    """Beer-Lambert UV detection tests."""

    def test_uv_beer_lambert(self):
        """A = epsilon * c * l for known concentration.

        C = 1.0 mol/m^3 = 0.001 M
        epsilon = 36000 1/(M*cm)
        l = 0.01 m = 1 cm
        A = 36000 * 0.001 * 1 = 36.0 AU = 36000 mAU
        """
        C = np.array([1.0])
        signal = compute_uv_signal(C, extinction_coeff=36000.0, path_length=0.01)
        assert signal[0] == pytest.approx(36000.0, rel=1e-6)

    def test_uv_zero_concentration(self):
        """Zero concentration gives zero signal."""
        C = np.array([0.0, 0.0])
        signal = compute_uv_signal(C)
        assert np.all(signal == 0.0)

    def test_uv_linearity(self):
        """Signal is linear in concentration."""
        C1 = np.array([1.0])
        C2 = np.array([2.0])
        s1 = compute_uv_signal(C1)
        s2 = compute_uv_signal(C2)
        assert s2[0] == pytest.approx(2.0 * s1[0], rel=1e-10)

    def test_broadening_preserves_area(self):
        """Gaussian broadening approximately preserves total area."""
        time = np.linspace(0, 100, 1000)
        # Sharp peak
        signal = np.zeros(1000)
        signal[450:550] = 100.0

        broadened = apply_detector_broadening(signal, time, sigma_detector=2.0)
        area_original = np.trapezoid(signal, time)
        area_broadened = np.trapezoid(broadened, time)
        assert area_broadened == pytest.approx(area_original, rel=0.05)

    def test_broadening_zero_sigma(self):
        """Zero sigma returns copy of original."""
        signal = np.array([1.0, 2.0, 3.0])
        time = np.array([0.0, 1.0, 2.0])
        result = apply_detector_broadening(signal, time, sigma_detector=0.0)
        np.testing.assert_array_equal(result, signal)


# ─── C.5: Orchestrator ───────────────────────────────────────────────

class TestBreakthroughResult:
    """run_breakthrough integration tests."""

    def test_breakthrough_result_complete(self):
        """run_breakthrough returns all required fields."""
        col = ColumnGeometry()
        result = run_breakthrough(
            column=col,
            C_feed=1.0,
            flow_rate=1e-8,
            feed_duration=300.0,
            total_time=600.0,
            n_z=30,
        )
        assert isinstance(result, BreakthroughResult)
        assert result.time.ndim == 1
        assert result.uv_signal.shape == result.time.shape
        assert result.C_outlet.shape == result.time.shape
        assert np.isfinite(result.dbc_5pct)
        assert np.isfinite(result.dbc_10pct)
        assert np.isfinite(result.dbc_50pct)
        assert result.pressure_drop > 0
        assert np.isfinite(result.mass_balance_error)

    def test_dbc_ordering(self):
        """DBC_5% <= DBC_10% <= DBC_50%."""
        col = ColumnGeometry(bed_height=0.05)
        iso = LangmuirIsotherm(q_max=50.0, K_L=500.0)
        result = run_breakthrough(
            column=col,
            C_feed=1.0,
            flow_rate=1e-8,
            feed_duration=600.0,
            total_time=800.0,
            isotherm=iso,
            n_z=30,
        )
        # DBC at lower threshold should be <= DBC at higher threshold
        # because more protein is loaded by the time higher breakthrough occurs
        assert result.dbc_5pct <= result.dbc_10pct + 1e-10, (
            f"DBC_5%={result.dbc_5pct:.4f} > DBC_10%={result.dbc_10pct:.4f}"
        )
        assert result.dbc_10pct <= result.dbc_50pct + 1e-10, (
            f"DBC_10%={result.dbc_10pct:.4f} > DBC_50%={result.dbc_50pct:.4f}"
        )

    def test_pressure_drop_matches_column(self):
        """Pressure drop in result matches column.pressure_drop()."""
        col = ColumnGeometry()
        Q = 1e-8
        result = run_breakthrough(
            column=col,
            C_feed=1.0,
            flow_rate=Q,
            feed_duration=100.0,
            total_time=200.0,
            n_z=20,
        )
        expected_dP = col.pressure_drop(Q)
        assert result.pressure_drop == pytest.approx(expected_dP, rel=1e-10)


# ─── Node 5 (v6.1): M3 ModelManifest evidence wiring ──────────────────────


class TestM3ModelManifest:
    """Verifies that ModelManifest fields are populated end-to-end through M3.

    Acceptance for Node 5:
      - run_breakthrough attaches a manifest to BreakthroughResult.
      - When an FMC with a downgraded tier is supplied, the breakthrough
        manifest cannot be stronger than the FMC tier (evidence inheritance).
      - run_gradient_elution attaches a manifest with mass-balance gating.
      - solve_packed_bed (catalysis) attaches a manifest with regime
        classification (reaction_limited / borderline / diffusion_limited).
    """

    def test_breakthrough_has_manifest_default_semi(self, default_column):
        """No FMC, default isotherm -> SEMI_QUANTITATIVE manifest."""
        from emulsim.datatypes import ModelEvidenceTier, ModelManifest

        result = run_breakthrough(
            column=default_column,
            C_feed=1.0,
            flow_rate=1e-8,
            feed_duration=200.0,
            total_time=400.0,
            n_z=20,
        )
        assert result.model_manifest is not None
        assert isinstance(result.model_manifest, ModelManifest)
        assert result.model_manifest.evidence_tier == ModelEvidenceTier.SEMI_QUANTITATIVE
        assert result.model_manifest.model_name == "M3.breakthrough.LRM"
        diag = result.model_manifest.diagnostics
        assert diag["fmc_provided"] is False
        # Diagnostics include the isotherm class so RunReport consumers can
        # tell which equilibrium model fed the breakthrough numbers.
        assert "Langmuir" in diag["isotherm_class"]
        assert diag["mass_balance_status"] in {"ok", "caution", "blocker"}

    def test_breakthrough_inherits_fmc_qualitative_tier(self, default_column):
        """FMC tagged QUALITATIVE_TREND -> breakthrough cannot be stronger."""
        from emulsim.datatypes import ModelEvidenceTier, ModelManifest
        from emulsim.module2_functionalization.orchestrator import (
            FunctionalMediaContract,
        )

        # Hand-build a minimal FMC carrying a ranking-only manifest. The
        # breakthrough call inherits from FMC.model_manifest.evidence_tier;
        # other FMC fields are ignored by the chrom code path here because
        # the caller passed isotherm explicitly via the LangmuirIsotherm
        # default + fmc-driven adapter (we leave isotherm=None and rely on
        # FMC routing falling back to default Langmuir).
        ranking_manifest = ModelManifest(
            model_name="M2.FMC.affinity",
            evidence_tier=ModelEvidenceTier.QUALITATIVE_TREND,
            assumptions=["Ranking-only ligand class."],
            diagnostics={"fmc_confidence_tier": "ranking_only"},
        )
        fmc = FunctionalMediaContract(
            bead_d50=default_column.particle_diameter,
            porosity=default_column.particle_porosity,
            estimated_q_max=100.0,
            confidence_tier="ranking_only",
            model_manifest=ranking_manifest,
        )

        result = run_breakthrough(
            column=default_column,
            C_feed=1.0,
            flow_rate=1e-8,
            feed_duration=200.0,
            total_time=400.0,
            n_z=20,
            fmc=fmc,
        )

        assert result.model_manifest is not None
        # Tier must not be stronger (lower index in _ORDER) than the upstream FMC.
        _ORDER = list(ModelEvidenceTier)
        upstream_idx = _ORDER.index(ModelEvidenceTier.QUALITATIVE_TREND)
        result_idx = _ORDER.index(result.model_manifest.evidence_tier)
        assert result_idx >= upstream_idx, (
            f"M3 tier {result.model_manifest.evidence_tier} stronger than FMC tier "
            f"QUALITATIVE_TREND — evidence inheritance broken."
        )
        assert result.model_manifest.diagnostics["fmc_provided"] is True

    def test_gradient_has_manifest(self, default_column):
        """Gradient elution attaches a manifest with per-component mass balance."""
        from emulsim.datatypes import ModelEvidenceTier
        from emulsim.module3_performance.gradient import GradientProgram
        from emulsim.module3_performance.isotherms.competitive_langmuir import (
            CompetitiveLangmuirIsotherm,
        )
        from emulsim.module3_performance.orchestrator import run_gradient_elution

        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([100.0, 100.0]),
            K_L=np.array([1e3, 5e2]),
        )
        # Simple linear salt gradient that completes within the run window.
        grad = GradientProgram(
            segments=[(0.0, 0.0), (100.0, 0.0), (300.0, 500.0), (400.0, 500.0)]
        )

        result = run_gradient_elution(
            column=default_column,
            C_feed=np.array([1.0, 1.0]),
            gradient=grad,
            flow_rate=1e-8,
            total_time=400.0,
            feed_duration=80.0,
            isotherm=iso,
            n_z=15,
        )

        assert result.model_manifest is not None
        assert result.model_manifest.model_name == "M3.gradient_elution.LRM"
        diag = result.model_manifest.diagnostics
        assert diag["n_components"] == 2
        # tier valid against the enum
        assert isinstance(result.model_manifest.evidence_tier, ModelEvidenceTier)

    def test_catalytic_manifest_classifies_regime(self):
        """Packed-bed enzyme reactor attaches a manifest with regime label."""
        from emulsim.datatypes import ModelEvidenceTier, ModelManifest
        from emulsim.module3_performance.catalysis.packed_bed import (
            solve_packed_bed,
        )

        # Reaction-limited regime: small particles + fast diffusion -> phi << 1.
        result = solve_packed_bed(
            bed_length=0.10,
            bed_diameter=0.01,
            particle_diameter=50e-6,
            bed_porosity=0.38,
            particle_porosity=0.5,
            V_max=10.0,
            K_m=1.0,
            S_feed=10.0,
            flow_rate=1e-8,
            D_eff=1e-9,
            total_time=600.0,
            n_z=20,
        )

        assert result.model_manifest is not None
        assert isinstance(result.model_manifest, ModelManifest)
        assert result.model_manifest.model_name == "M3.catalysis.packed_bed.MM"
        diag = result.model_manifest.diagnostics
        assert diag["regime"] in {"reaction_limited", "borderline", "diffusion_limited"}
        assert diag["mass_balance_status"] in {"ok", "caution", "blocker"}
        # Default tier is SEMI unless gates fire
        assert result.model_manifest.evidence_tier in {
            ModelEvidenceTier.SEMI_QUANTITATIVE,
            ModelEvidenceTier.QUALITATIVE_TREND,
            ModelEvidenceTier.UNSUPPORTED,
        }
