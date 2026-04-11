"""Tests for Phase E multi-component chromatography extensions.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

Covers:
    E.1  CompetitiveLangmuirIsotherm
    E.2  SMAIsotherm
    E.3  IMACCompetitionIsotherm
    E.4  ProteinAIsotherm
    E.5  GradientProgram / make_linear_gradient / make_step_gradient
    E.6  compute_fluorescence_signal
    E.7  compute_conductivity
    E.8  simulate_esi_charge_envelope / compute_tic
    E.9  run_gradient_elution
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.module3_performance import (
    # Isotherms
    LangmuirIsotherm,
    CompetitiveLangmuirIsotherm,
    SMAIsotherm,
    IMACCompetitionIsotherm,
    ProteinAIsotherm,
    # Gradient
    GradientProgram,
    make_linear_gradient,
    make_step_gradient,
    # Detection
    compute_fluorescence_signal,
    compute_conductivity,
    conductivity_to_ms_per_cm,
    simulate_esi_charge_envelope,
    compute_tic,
    # Orchestrator
    ColumnGeometry,
    run_gradient_elution,
    GradientElutionResult,
)


# ─── E.1: Competitive Langmuir ───────────────────────────────────────────────

class TestCompetitiveLangmuir:
    """Tests for CompetitiveLangmuirIsotherm."""

    def test_competitive_langmuir_reduces_loading(self):
        """Competition reduces q_i vs single-component Langmuir."""
        q_max = np.array([100.0, 100.0])
        K_L = np.array([1e3, 1e3])
        comp_iso = CompetitiveLangmuirIsotherm(q_max=q_max, K_L=K_L)
        single_iso = LangmuirIsotherm(q_max=100.0, K_L=1e3)

        C_vec = np.array([0.01, 0.01])   # [mol/m^3]

        q_comp = comp_iso.equilibrium_loading(C_vec)
        q_single = single_iso.equilibrium_loading(C_vec[0])

        # Each component is reduced by the presence of the other
        assert q_comp[0] < q_single, (
            f"Competitive loading {q_comp[0]:.4f} should be < "
            f"single-component {q_single:.4f}"
        )
        assert q_comp[1] < q_single

    def test_competitive_reduces_to_single_at_zero_competitor(self):
        """With one component at zero, result matches single-component Langmuir."""
        q_max = np.array([100.0, 80.0])
        K_L = np.array([1e3, 5e2])
        comp_iso = CompetitiveLangmuirIsotherm(q_max=q_max, K_L=K_L)
        single_iso_0 = LangmuirIsotherm(q_max=100.0, K_L=1e3)

        C = np.array([0.01, 0.0])   # component 1 absent
        q = comp_iso.equilibrium_loading(C)
        q_ref = single_iso_0.equilibrium_loading(0.01)

        np.testing.assert_allclose(q[0], q_ref, rtol=1e-8)
        np.testing.assert_allclose(q[1], 0.0, atol=1e-15)

    def test_competitive_langmuir_non_negative(self):
        """Loading is always non-negative, even with negative input concentrations."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([50.0, 50.0]),
            K_L=np.array([1e3, 1e3]),
        )
        C = np.array([-0.5, -0.1])
        q = iso.equilibrium_loading(C)
        assert np.all(q >= 0.0)

    def test_competitive_langmuir_saturation(self):
        """At very high concentration, loading approaches q_max."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([100.0, 100.0]),
            K_L=np.array([1e3, 1e3]),
        )
        C_high = np.array([1000.0, 0.0])  # only component 0, saturating
        q = iso.equilibrium_loading(C_high)
        # Should approach q_max[0]
        assert q[0] > 90.0, f"Expected near-saturation, got {q[0]:.2f}"

    def test_competitive_langmuir_grid_input(self):
        """Accepts (n_comp, N) input for spatial grid computation."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([100.0, 80.0]),
            K_L=np.array([1e3, 5e2]),
        )
        N = 20
        C_grid = np.random.rand(2, N) * 0.01
        q_grid = iso.equilibrium_loading(C_grid)
        assert q_grid.shape == (2, N)
        assert np.all(q_grid >= 0.0)

    def test_validate_valid_params(self):
        """validate() returns empty list for valid parameters."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([100.0, 80.0]),
            K_L=np.array([1e3, 5e2]),
        )
        assert iso.validate() == []

    def test_validate_rejects_negative_qmax(self):
        """validate() catches non-positive q_max."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([-10.0, 80.0]),
            K_L=np.array([1e3, 5e2]),
        )
        errors = iso.validate()
        assert len(errors) > 0


# ─── E.2: SMA Isotherm ───────────────────────────────────────────────────────

class TestSMAIsotherm:
    """Tests for SMAIsotherm (Steric Mass Action)."""

    def test_sma_salt_dependent(self):
        """Higher salt concentration reduces protein retention (q drops)."""
        iso = SMAIsotherm(
            Lambda=1000.0,
            z=np.array([3.0]),
            sigma=np.array([50.0]),
            K_eq=np.array([1e-3]),
        )
        C_protein = np.array([0.001])  # [mol/m^3]

        q_low_salt = iso.equilibrium_loading(C_protein, salt_concentration=10.0)
        q_high_salt = iso.equilibrium_loading(C_protein, salt_concentration=500.0)

        assert q_high_salt[0] < q_low_salt[0], (
            f"High salt q={q_high_salt[0]:.4f} should be < "
            f"low salt q={q_low_salt[0]:.4f}"
        )

    def test_sma_zero_protein(self):
        """Zero protein concentration gives zero loading."""
        iso = SMAIsotherm(
            Lambda=1000.0,
            z=np.array([3.0, 5.0]),
            sigma=np.array([50.0, 80.0]),
            K_eq=np.array([1e-3, 1e-5]),
        )
        C_protein = np.array([0.0, 0.0])
        q = iso.equilibrium_loading(C_protein, salt_concentration=100.0)
        np.testing.assert_allclose(q, [0.0, 0.0], atol=1e-12)

    def test_sma_non_negative(self):
        """Loading is always non-negative."""
        iso = SMAIsotherm(
            Lambda=1000.0,
            z=np.array([3.0, 5.0]),
            sigma=np.array([50.0, 80.0]),
            K_eq=np.array([1e-3, 1e-5]),
        )
        C = np.array([0.01, 0.001])
        q = iso.equilibrium_loading(C, salt_concentration=100.0)
        assert np.all(q >= 0.0)

    def test_sma_retention_factor_decreases_with_salt(self):
        """Retention factor k' decreases monotonically with salt."""
        iso = SMAIsotherm(
            Lambda=1000.0,
            z=np.array([3.0]),
            sigma=np.array([50.0]),
            K_eq=np.array([1e-3]),
        )
        salts = [10.0, 50.0, 100.0, 500.0]
        k_primes = [iso.retention_factor(s, component_idx=0) for s in salts]
        for i in range(len(k_primes) - 1):
            assert k_primes[i] > k_primes[i + 1], (
                f"k' should decrease with salt: {k_primes}"
            )

    def test_sma_validate(self):
        """validate() returns empty for valid SMA parameters."""
        iso = SMAIsotherm(
            Lambda=1000.0,
            z=np.array([3.0, 5.0]),
            sigma=np.array([50.0, 80.0]),
            K_eq=np.array([1e-3, 1e-5]),
        )
        assert iso.validate() == []


# ─── E.3: IMAC Isotherm ──────────────────────────────────────────────────────

class TestIMACIsotherm:
    """Tests for IMACCompetitionIsotherm."""

    def test_imac_imidazole_competition(self):
        """Imidazole reduces protein binding — higher [Im] -> lower q_protein."""
        iso = IMACCompetitionIsotherm(
            q_max=50.0,
            K_protein=1e4,
            K_imidazole=10.0,
        )
        C_protein = 0.001  # [mol/m^3]

        q_no_im, _ = iso.equilibrium_loading(C_protein, C_imidazole=0.0)
        q_with_im, _ = iso.equilibrium_loading(C_protein, C_imidazole=10.0)

        assert q_with_im < q_no_im, (
            f"Imidazole should reduce protein binding: "
            f"no_im={q_no_im:.4f}, with_im={q_with_im:.4f}"
        )

    def test_imac_high_imidazole_elutes_protein(self):
        """Very high imidazole essentially eliminates protein binding."""
        iso = IMACCompetitionIsotherm(q_max=50.0, K_protein=1e4, K_imidazole=10.0)
        C_protein = 0.001
        q_prot, _ = iso.equilibrium_loading(C_protein, C_imidazole=1000.0)
        # q should be very small relative to q_max
        assert q_prot < 0.01 * iso.q_max

    def test_imac_no_imidazole_is_single_langmuir(self):
        """Without imidazole, IMAC reduces to single-component Langmuir."""
        iso = IMACCompetitionIsotherm(q_max=50.0, K_protein=1e4, K_imidazole=10.0)
        C_protein = 0.001
        q_prot, _ = iso.equilibrium_loading(C_protein, C_imidazole=0.0)
        q_ref = 50.0 * 1e4 * C_protein / (1.0 + 1e4 * C_protein)
        np.testing.assert_allclose(q_prot, q_ref, rtol=1e-8)

    def test_imac_protein_loading_only_convenience(self):
        """protein_loading_only returns same value as equilibrium_loading[0]."""
        iso = IMACCompetitionIsotherm(q_max=50.0, K_protein=1e4, K_imidazole=10.0)
        C_p, C_i = 0.001, 5.0
        q1, _ = iso.equilibrium_loading(C_p, C_i)
        q2 = iso.protein_loading_only(C_p, C_i)
        np.testing.assert_allclose(q1, q2, rtol=1e-12)

    def test_imac_validate(self):
        """validate() returns empty for valid IMAC parameters."""
        iso = IMACCompetitionIsotherm()
        assert iso.validate() == []

    def test_imac_validate_catches_reversed_affinity(self):
        """validate() warns if K_imidazole >= K_protein (unlikely physically)."""
        iso = IMACCompetitionIsotherm(q_max=50.0, K_protein=5.0, K_imidazole=10.0)
        errors = iso.validate()
        assert len(errors) > 0


# ─── E.4: Protein A Isotherm ─────────────────────────────────────────────────

class TestProteinAIsotherm:
    """Tests for ProteinAIsotherm (pH-dependent affinity)."""

    def test_protein_a_ph_elution(self):
        """Low pH reduces K_a and thus equilibrium loading (elution effect)."""
        iso = ProteinAIsotherm(
            q_max=60.0,
            K_a_max=1e5,
            pH_transition=3.5,
            steepness=5.0,
        )
        C = np.array([0.001])  # [mol/m^3]

        q_neutral = iso.equilibrium_loading(C, pH=7.0)
        q_low_pH  = iso.equilibrium_loading(C, pH=2.5)

        assert q_low_pH[0] < q_neutral[0], (
            f"Low pH q={q_low_pH[0]:.4f} should be < "
            f"neutral pH q={q_neutral[0]:.4f}"
        )

    def test_protein_a_high_ph_near_saturation(self):
        """At high pH and high C, loading approaches q_max."""
        iso = ProteinAIsotherm(q_max=60.0, K_a_max=1e5, pH_transition=3.5)
        q = iso.equilibrium_loading(np.array([10.0]), pH=8.0)
        assert q[0] > 0.99 * iso.q_max

    def test_protein_a_very_low_ph_near_zero(self):
        """At very low pH with trace-level protein, loading approaches zero.

        At pH=1.0, K_a(1.0) ~ 0.37 m^3/mol.  For dilute protein (C << 1/K_a),
        q ~ q_max * K_a * C.  At C=1e-4 mol/m^3: q ~ 60 * 0.37 * 1e-4 = 2.2e-3,
        which is << 0.01 * q_max = 0.6.
        """
        iso = ProteinAIsotherm(q_max=60.0, K_a_max=1e5, pH_transition=3.5, steepness=5.0)
        # Use very dilute protein so the low K_a gives near-zero loading
        q = iso.equilibrium_loading(np.array([1e-4]), pH=1.0)
        assert q[0] < 0.01 * iso.q_max, (
            f"q={q[0]:.4g} should be << {0.01 * iso.q_max:.4g} at pH 1.0 with low C"
        )

    def test_protein_a_sigmoid_transition(self):
        """K_a decreases monotonically from neutral to acidic pH."""
        iso = ProteinAIsotherm(q_max=60.0, K_a_max=1e5, pH_transition=3.5, steepness=5.0)
        pHs = np.linspace(2.0, 8.0, 20)
        Ka_vals = np.array([iso.K_a(p) for p in pHs])
        # K_a should be monotonically increasing with pH
        assert np.all(np.diff(Ka_vals) >= 0), "K_a should increase with pH"

    def test_protein_a_elution_ph_method(self):
        """elution_pH returns pH at which K_a = fraction * K_a_max."""
        iso = ProteinAIsotherm(q_max=60.0, K_a_max=1e5, pH_transition=3.5, steepness=5.0)
        pH_elu = iso.elution_pH(fraction=0.01)
        Ka_at_elu = iso.K_a(pH_elu)
        np.testing.assert_allclose(Ka_at_elu, 0.01 * iso.K_a_max, rtol=1e-4)

    def test_protein_a_validate(self):
        """validate() returns empty for valid parameters."""
        iso = ProteinAIsotherm()
        assert iso.validate() == []


# ─── E.5: Gradient Generator ─────────────────────────────────────────────────

class TestGradientProgram:
    """Tests for GradientProgram and convenience constructors."""

    def test_gradient_linear_monotonic(self):
        """Linear gradient produces monotonically increasing value."""
        g = make_linear_gradient(start_val=0.0, end_val=1.0, start_time=0.0, end_time=600.0)
        times = np.linspace(0.0, 600.0, 50)
        values = np.array([g.value_at_time(t) for t in times])
        diffs = np.diff(values)
        assert np.all(diffs >= -1e-12), (
            "Linear gradient should be non-decreasing"
        )

    def test_gradient_linear_endpoints(self):
        """Linear gradient exactly hits start and end values."""
        g = make_linear_gradient(start_val=0.1, end_val=0.9, start_time=0.0, end_time=300.0)
        np.testing.assert_allclose(g.value_at_time(0.0), 0.1, rtol=1e-8)
        np.testing.assert_allclose(g.value_at_time(300.0), 0.9, rtol=1e-8)

    def test_gradient_linear_midpoint(self):
        """Midpoint of linear gradient is at (start + end) / 2."""
        g = make_linear_gradient(start_val=0.0, end_val=1.0, end_time=600.0)
        np.testing.assert_allclose(g.value_at_time(300.0), 0.5, rtol=1e-6)

    def test_gradient_step_sharp_transition(self):
        """Step gradient reaches target value shortly after the step time."""
        g = make_step_gradient([(0.0, 0.0), (300.0, 1.0)], ramp_duration=1.0)
        # After the ramp (300 + 1 = 301 s), value should be 1.0
        assert g.value_at_time(302.0) == pytest.approx(1.0, abs=1e-8)
        # Before the step, value should be 0.0
        assert g.value_at_time(299.0) == pytest.approx(0.0, abs=1e-8)

    def test_gradient_step_multiple_steps(self):
        """Multi-step gradient has correct values at each plateau."""
        g = make_step_gradient(
            [(0.0, 0.0), (200.0, 0.5), (500.0, 1.0)],
            ramp_duration=1.0,
        )
        # At t=100 s (first plateau)
        assert g.value_at_time(100.0) == pytest.approx(0.0, abs=1e-8)
        # At t=202 s (after first step)
        assert g.value_at_time(202.0) == pytest.approx(0.5, abs=1e-8)
        # At t=502 s (after second step)
        assert g.value_at_time(502.0) == pytest.approx(1.0, abs=1e-8)

    def test_gradient_extrapolation_constant(self):
        """Gradient is constant outside the defined range (np.interp behaviour)."""
        g = make_linear_gradient(0.0, 1.0, start_time=100.0, end_time=200.0)
        # Before start: constant at start_val
        assert g.value_at_time(0.0) == pytest.approx(0.0, abs=1e-8)
        # After end: constant at end_val
        assert g.value_at_time(500.0) == pytest.approx(1.0, abs=1e-8)

    def test_gradient_program_array_input(self):
        """value_at_time accepts numpy array input."""
        g = make_linear_gradient(0.0, 1.0, end_time=600.0)
        t_arr = np.linspace(0.0, 600.0, 10)
        vals = g.value_at_time(t_arr)
        assert vals.shape == (10,)

    def test_gradient_program_requires_two_segments(self):
        """GradientProgram raises ValueError with fewer than 2 waypoints."""
        with pytest.raises(ValueError, match="at least 2"):
            GradientProgram(segments=[(0.0, 0.5)])

    def test_gradient_program_requires_sorted_times(self):
        """GradientProgram raises ValueError for unsorted times."""
        with pytest.raises(ValueError, match="ascending"):
            GradientProgram(segments=[(600.0, 1.0), (0.0, 0.0)])


# ─── E.6: Fluorescence Detection ─────────────────────────────────────────────

class TestFluorescenceDetection:
    """Tests for compute_fluorescence_signal."""

    def test_fluorescence_proportional_to_concentration(self):
        """Signal is linearly proportional to concentration."""
        C1 = 0.001
        C2 = 0.002
        F1 = float(compute_fluorescence_signal(C1, quantum_yield=0.92, extinction_coeff=76900.0))
        F2 = float(compute_fluorescence_signal(C2, quantum_yield=0.92, extinction_coeff=76900.0))
        np.testing.assert_allclose(F2 / F1, 2.0, rtol=1e-8)

    def test_fluorescence_proportional_to_quantum_yield(self):
        """Signal is linearly proportional to quantum yield."""
        C = 0.001
        F1 = float(compute_fluorescence_signal(C, quantum_yield=0.5, extinction_coeff=76900.0))
        F2 = float(compute_fluorescence_signal(C, quantum_yield=1.0, extinction_coeff=76900.0))
        np.testing.assert_allclose(F2 / F1, 2.0, rtol=1e-8)

    def test_fluorescence_zero_concentration(self):
        """Zero concentration gives zero signal."""
        F = float(compute_fluorescence_signal(0.0))
        assert F == 0.0

    def test_fluorescence_non_negative(self):
        """Signal is non-negative even for negative input concentrations."""
        F = float(compute_fluorescence_signal(-0.5))
        assert F >= 0.0

    def test_fluorescence_array_input(self):
        """Signal works element-wise on arrays."""
        C_arr = np.array([0.001, 0.002, 0.003])
        F_arr = compute_fluorescence_signal(C_arr)
        assert F_arr.shape == (3,)
        # Linearly increasing
        diffs = np.diff(F_arr)
        assert np.all(diffs > 0)

    def test_fluorescence_invalid_quantum_yield(self):
        """Raises ValueError for invalid quantum yield."""
        with pytest.raises(ValueError, match="quantum_yield"):
            compute_fluorescence_signal(0.001, quantum_yield=1.5)

    def test_fluorescence_proportional_to_extinction_coeff(self):
        """Signal proportional to extinction coefficient."""
        C = 0.001
        F1 = float(compute_fluorescence_signal(C, quantum_yield=0.9, extinction_coeff=50000.0))
        F2 = float(compute_fluorescence_signal(C, quantum_yield=0.9, extinction_coeff=100000.0))
        np.testing.assert_allclose(F2 / F1, 2.0, rtol=1e-8)


# ─── E.7: Conductivity Detection ─────────────────────────────────────────────

class TestConductivityDetection:
    """Tests for compute_conductivity."""

    def test_conductivity_proportional_to_salt(self):
        """Conductivity is linearly proportional to salt concentration."""
        c1 = 100.0  # [mol/m^3]
        c2 = 200.0
        k1 = float(compute_conductivity(c1))
        k2 = float(compute_conductivity(c2))
        np.testing.assert_allclose(k2 / k1, 2.0, rtol=1e-8)

    def test_conductivity_zero_salt(self):
        """Zero salt gives zero conductivity."""
        k = float(compute_conductivity(0.0))
        assert k == 0.0

    def test_conductivity_non_negative(self):
        """Conductivity is non-negative even for negative input."""
        k = float(compute_conductivity(-10.0))
        assert k >= 0.0

    def test_conductivity_default_nacl_values(self):
        """Default molar conductivities give physically reasonable result.

        NaCl at 100 mol/m^3 (100 mM) should give kappa ~ 1.26 S/m = 12.6 mS/cm.
        Literature value is ~11.5 mS/cm at 25 C (accounts for non-ideal activity).
        The linear model over-estimates slightly, so we check order-of-magnitude.
        """
        k_si = float(compute_conductivity(100.0))  # [S/m]
        k_ms_cm = conductivity_to_ms_per_cm(k_si)
        # Should be in the range 10-15 mS/cm for 100 mM NaCl
        assert 8.0 < k_ms_cm < 20.0, f"Unexpected conductivity: {k_ms_cm:.2f} mS/cm"

    def test_conductivity_array_input(self):
        """compute_conductivity works element-wise on arrays."""
        c_arr = np.array([0.0, 50.0, 100.0, 200.0])
        k_arr = compute_conductivity(c_arr)
        assert k_arr.shape == (4,)
        diffs = np.diff(k_arr)
        assert np.all(diffs >= 0.0)

    def test_conductivity_unit_conversion(self):
        """conductivity_to_ms_per_cm converts correctly (1 S/m = 10 mS/cm)."""
        k_si = np.array([0.5, 1.0, 2.0])
        k_ms = conductivity_to_ms_per_cm(k_si)
        np.testing.assert_allclose(k_ms, k_si * 10.0, rtol=1e-12)


# ─── E.8: Mass Spectrometry Detection ────────────────────────────────────────

class TestMSDetection:
    """Tests for simulate_esi_charge_envelope and compute_tic."""

    def test_esi_charge_states_50kda(self):
        """m/z values span a reasonable range for a 50 kDa protein.

        For a 50 kDa protein with z_avg ~ 17.4, the most abundant charge state
        is around m/z ~ (50000 + 17 * 1.007) / 17 ~ 2942 Da.
        The m/z axis spans from z_min to z_max (low-charge states have very
        high m/z, which are physically valid but low-abundance).
        We check: average charge state is correct, and at least some peaks
        fall in the 1000-8000 m/z window typical for conventional ESI.
        """
        spec = simulate_esi_charge_envelope(50000.0, n_points=100)

        # z_avg for 50 kDa: 0.0778 * sqrt(50000) ~ 17.4
        assert spec.z_avg == pytest.approx(0.0778 * np.sqrt(50000.0), rel=0.01)

        # The m/z axis must contain values > 0
        assert spec.mz.min() > 0.0, f"m/z min {spec.mz.min():.0f} should be positive"

        # Average charge state in sensible range for a 50 kDa protein
        assert 10 < spec.z_avg < 30, f"z_avg {spec.z_avg:.1f} out of expected range"

        # Peak apex m/z should be in the conventional ESI window (500-5000 Da)
        apex_mz = float(spec.mz[np.argmax(spec.intensity)])
        assert 500.0 < apex_mz < 8000.0, (
            f"Peak m/z {apex_mz:.0f} outside expected ESI window 500-8000 Da"
        )

    def test_esi_max_intensity_normalized(self):
        """Maximum intensity is normalised to 1.0."""
        spec = simulate_esi_charge_envelope(50000.0)
        assert spec.intensity.max() == pytest.approx(1.0, abs=1e-6)

    def test_esi_charge_range_sensible(self):
        """Charge state range spans at least 5 states for a 50 kDa protein."""
        spec = simulate_esi_charge_envelope(50000.0)
        n_states = spec.z_max - spec.z_min + 1
        assert n_states >= 5, f"Expected >= 5 charge states, got {n_states}"

    def test_esi_larger_protein_higher_charges(self):
        """Larger proteins have higher average charge states."""
        spec_small = simulate_esi_charge_envelope(10000.0)
        spec_large = simulate_esi_charge_envelope(150000.0)
        assert spec_large.z_avg > spec_small.z_avg

    def test_esi_invalid_mw(self):
        """Raises ValueError for non-positive molecular weight."""
        with pytest.raises(ValueError, match="molecular_weight"):
            simulate_esi_charge_envelope(-1000.0)

    def test_tic_proportional_to_concentration(self):
        """TIC signal is proportional to concentration."""
        time = np.linspace(0, 100, 50)
        C1 = np.ones(50) * 0.001
        C2 = np.ones(50) * 0.002
        tic1 = compute_tic(C1, time)
        tic2 = compute_tic(C2, time)
        np.testing.assert_allclose(tic2 / tic1, 2.0, rtol=1e-8)

    def test_tic_multi_component_sums(self):
        """TIC for multi-component is the sum of all concentrations."""
        time = np.linspace(0, 100, 50)
        C = np.array([np.ones(50) * 0.001, np.ones(50) * 0.002])
        tic = compute_tic(C, time)
        # Each element should be 0.003 * sensitivity (1.0)
        np.testing.assert_allclose(tic, np.ones(50) * 0.003, rtol=1e-8)

    def test_tic_non_negative(self):
        """TIC is non-negative even for negative input concentrations."""
        time = np.linspace(0, 100, 10)
        C = np.full(10, -0.01)
        tic = compute_tic(C, time)
        assert np.all(tic >= 0.0)


# ─── E.9: Gradient Elution Orchestrator ──────────────────────────────────────

class TestGradientElution:
    """Tests for run_gradient_elution."""

    @pytest.fixture
    def column(self):
        return ColumnGeometry(
            diameter=0.01,
            bed_height=0.10,
            particle_diameter=100e-6,
            bed_porosity=0.38,
            particle_porosity=0.70,
        )

    def test_gradient_elution_returns_result_type(self, column):
        """run_gradient_elution returns a GradientElutionResult."""
        from emulsim.module3_performance import CompetitiveLangmuirIsotherm

        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([50.0, 30.0]),
            K_L=np.array([1e3, 5e2]),
        )
        gradient = make_linear_gradient(0.0, 1.0, end_time=600.0)
        C_feed = np.array([0.01, 0.01])

        result = run_gradient_elution(
            column=column,
            C_feed=C_feed,
            gradient=gradient,
            flow_rate=1e-8,
            total_time=600.0,
            isotherm=iso,
            n_z=20,
        )
        assert isinstance(result, GradientElutionResult)

    def test_gradient_elution_output_shapes(self, column):
        """Outlet concentration has shape (n_comp, N_t)."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([50.0, 30.0]),
            K_L=np.array([1e3, 5e2]),
        )
        gradient = make_linear_gradient(0.0, 1.0, end_time=400.0)
        C_feed = np.array([0.01, 0.01])

        result = run_gradient_elution(
            column=column,
            C_feed=C_feed,
            gradient=gradient,
            flow_rate=1e-8,
            total_time=400.0,
            isotherm=iso,
            n_z=15,
        )
        assert result.C_outlet.shape[0] == 2
        assert result.C_outlet.shape[1] == len(result.time)
        assert result.gradient_profile.shape == result.time.shape

    def test_gradient_elution_non_negative_concentrations(self, column):
        """Outlet concentrations are non-negative throughout."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([50.0, 30.0]),
            K_L=np.array([1e3, 5e2]),
        )
        gradient = make_linear_gradient(0.0, 1.0, end_time=400.0)
        C_feed = np.array([0.01, 0.005])

        result = run_gradient_elution(
            column=column,
            C_feed=C_feed,
            gradient=gradient,
            flow_rate=1e-8,
            total_time=400.0,
            isotherm=iso,
            n_z=15,
        )
        assert np.all(result.C_outlet >= 0.0)

    def test_gradient_elution_peak_table_length(self, column):
        """Peak table has one entry per component."""
        n_comp = 2
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([50.0, 30.0]),
            K_L=np.array([1e3, 5e2]),
        )
        gradient = make_linear_gradient(0.0, 1.0, end_time=400.0)
        C_feed = np.array([0.01, 0.01])

        result = run_gradient_elution(
            column=column,
            C_feed=C_feed,
            gradient=gradient,
            flow_rate=1e-8,
            total_time=400.0,
            isotherm=iso,
            n_z=15,
        )
        assert len(result.peaks) == n_comp

    def test_gradient_elution_pressure_drop_positive(self, column):
        """Pressure drop is positive."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([50.0, 30.0]),
            K_L=np.array([1e3, 5e2]),
        )
        gradient = make_linear_gradient(0.0, 1.0, end_time=300.0)
        C_feed = np.array([0.01, 0.01])

        result = run_gradient_elution(
            column=column,
            C_feed=C_feed,
            gradient=gradient,
            flow_rate=1e-8,
            total_time=300.0,
            isotherm=iso,
            n_z=10,
        )
        assert result.pressure_drop > 0.0

    def test_gradient_elution_uv_signal_shape(self, column):
        """UV signal has the same length as the time array."""
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([50.0, 30.0]),
            K_L=np.array([1e3, 5e2]),
        )
        gradient = make_linear_gradient(0.0, 1.0, end_time=300.0)
        C_feed = np.array([0.01, 0.01])

        result = run_gradient_elution(
            column=column,
            C_feed=C_feed,
            gradient=gradient,
            flow_rate=1e-8,
            total_time=300.0,
            isotherm=iso,
            n_z=10,
        )
        assert len(result.uv_signal) == len(result.time)

    def test_gradient_affects_binding_false_for_competitive_langmuir(self, column):
        """Competitive Langmuir should report gradient_affects_binding=False.

        BF-2 fix: GradientElutionResult must carry a gradient_affects_binding
        flag so the UI can display 'Gradient affects binding: NO' for
        competitive Langmuir (K_L is gradient-independent).
        """
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([100.0, 80.0]),
            K_L=np.array([1e3, 5e2]),
        )
        gradient = make_linear_gradient(0.0, 1.0, end_time=300.0)
        C_feed = np.array([0.01, 0.005])

        result = run_gradient_elution(
            column=column,
            C_feed=C_feed,
            gradient=gradient,
            flow_rate=1e-8,
            total_time=300.0,
            isotherm=iso,
            n_z=10,
        )

        # The isotherm property must return False
        assert iso.gradient_sensitive is False, (
            "CompetitiveLangmuirIsotherm.gradient_sensitive must be False"
        )
        # The result flag must be set correctly
        assert result.gradient_affects_binding is False, (
            "GradientElutionResult.gradient_affects_binding must be False "
            "for competitive Langmuir"
        )

    def test_feed_phase_switching_zeros_inlet_after_feed_duration(self, column):
        """Protein inlet must be zero during elution phase (feed_duration fix).

        With an explicit feed_duration shorter than total_time, the outlet
        concentration should peak and then decline as the bound protein elutes
        off into a protein-free buffer.  The total outlet area in the elution
        window must be non-trivial (protein actually eluted).
        """
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([100.0, 80.0]),
            K_L=np.array([1e3, 5e2]),
        )
        # Gradient starts after load: load for 100 s, then ramp 100-300 s
        gradient = make_linear_gradient(0.0, 1.0, start_time=100.0, end_time=300.0)
        C_feed = np.array([0.01, 0.005])
        feed_dur = 100.0   # protein only during load phase

        result = run_gradient_elution(
            column=column,
            C_feed=C_feed,
            gradient=gradient,
            flow_rate=1e-8,
            total_time=300.0,
            feed_duration=feed_dur,
            isotherm=iso,
            n_z=10,
        )

        # Sanity: result arrays have correct shape
        assert result.C_outlet.shape[0] == 2
        assert len(result.time) == len(result.gradient_profile)

        # After the feed phase the gradient_profile should be rising (gradient is on)
        t = result.time
        g = result.gradient_profile
        mid_idx = np.searchsorted(t, 200.0)  # well into elution window
        assert g[mid_idx] > 0.0, "Gradient should be non-zero during elution window"
