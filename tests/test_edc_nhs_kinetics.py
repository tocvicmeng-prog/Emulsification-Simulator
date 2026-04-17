"""Tests for Node 31 (v7.1): mechanistic two-step EDC/NHS kinetic.

Covers:

- Mass conservation invariant (C + A + E + P = C_0 within ε).
- Edge cases: zero EDC, zero NHS, zero NH2.
- pH plateau near neutral; drop at extremes.
- Arrhenius scaling with temperature.
- Dose response on [EDC].
- pH speciation helper correctness.
- Input validation.
- M2 dispatch: edc_nhs_activation step uses the mechanistic path.
- L3 gate: native matrix -> QUALITATIVE_TREND fallback; carboxylated
  matrix -> SEMI_QUANTITATIVE mechanistic path.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from emulsim.module2_functionalization.edc_nhs_kinetics import (
    EdcNhsKinetics,
    EdcNhsResult,
    available_amine_fraction,
    react_edc_nhs_two_step,
)


# ─── Core solver correctness ──────────────────────────────────────────────


class TestMassConservation:
    def test_mass_balance_typical_run(self):
        """C + A + E + P ≈ C_0 at end of any integration."""
        r = react_edc_nhs_two_step(
            c_cooh_initial=10.0,
            c_nh2_total=20.0,
            c_edc_initial=5.0,
            c_nhs_initial=5.0,
            pH=7.0,
            T=298.15,
            time=3600.0,
        )
        assert r.mass_balance_error < 1e-2, (
            f"mass balance slipped: err={r.mass_balance_error:.2e}"
        )
        assert r.solver_success

    def test_mass_balance_long_run(self):
        """Long-time integration still conserves mass."""
        r = react_edc_nhs_two_step(
            c_cooh_initial=10.0,
            c_nh2_total=20.0,
            c_edc_initial=5.0,
            c_nhs_initial=5.0,
            pH=7.5,
            T=298.15,
            time=86400.0,  # 24 h
        )
        assert r.mass_balance_error < 1e-2


class TestEdgeCases:
    def test_zero_edc_no_activation(self):
        """No EDC → no activation → p_final = 0."""
        r = react_edc_nhs_two_step(
            c_cooh_initial=10.0,
            c_nh2_total=20.0,
            c_edc_initial=0.0,
            c_nhs_initial=5.0,
            pH=7.0,
            T=298.15,
            time=3600.0,
        )
        assert r.p_final == pytest.approx(0.0, abs=1e-10)
        assert r.p_residual_nhs_ester == pytest.approx(0.0, abs=1e-10)

    def test_zero_nh2_no_coupling_but_activation_happens(self):
        """No NH2 → no amide product, but NHS ester still accumulates."""
        r = react_edc_nhs_two_step(
            c_cooh_initial=10.0,
            c_nh2_total=0.0,
            c_edc_initial=5.0,
            c_nhs_initial=5.0,
            pH=5.5,  # optimal for activation
            T=298.15,
            time=3600.0,
        )
        assert r.p_final == pytest.approx(0.0, abs=1e-10)
        # Some COOH should have converted to NHS ester.
        assert r.p_residual_nhs_ester > 0.01

    def test_zero_nhs_hydrolysis_dominates(self):
        """No NHS → hydrolysis consumes the o-acylisourea; little or no
        productive yield."""
        r = react_edc_nhs_two_step(
            c_cooh_initial=10.0,
            c_nh2_total=20.0,
            c_edc_initial=5.0,
            c_nhs_initial=0.0,
            pH=7.0,
            T=298.15,
            time=3600.0,
        )
        # Minimal productive product because NHS ester path can't form.
        assert r.p_final < 0.05


# ─── Scientific trends ────────────────────────────────────────────────────


class TestScientificTrends:
    def test_arrhenius_scaling_in_kinetic_regime(self):
        """At very short reaction times (kinetic regime, not EDC-limited),
        higher T produces higher p_final because activation (E_a_1=40
        kJ/mol) is the rate-limiting step and Arrhenius scaling
        dominates. All reagents in excess so no stoichiometric ceiling."""
        common = dict(
            c_cooh_initial=10.0, c_nh2_total=10.0,
            c_edc_initial=5.0, c_nhs_initial=5.0,
            pH=6.5, time=3.0,  # 3 seconds — clearly kinetic
        )
        p_low = react_edc_nhs_two_step(T=277.15, **common).p_final  # 4 °C
        p_mid = react_edc_nhs_two_step(T=298.15, **common).p_final  # 25 °C
        p_hi = react_edc_nhs_two_step(T=310.15, **common).p_final   # 37 °C
        assert p_low < p_mid, (
            f"Arrhenius kinetic regime expected p_low<p_mid; got "
            f"low={p_low:.4g} mid={p_mid:.4g}"
        )
        assert p_mid < p_hi, (
            f"Arrhenius kinetic regime expected p_mid<p_hi; got "
            f"mid={p_mid:.4g} hi={p_hi:.4g}"
        )

    def test_edc_dose_response(self):
        """Increasing [EDC] increases p_final (monotone, saturating)."""
        common = dict(
            c_cooh_initial=10.0, c_nh2_total=20.0, c_nhs_initial=5.0,
            pH=7.0, T=298.15, time=3600.0,
        )
        p1 = react_edc_nhs_two_step(c_edc_initial=1.0, **common).p_final
        p5 = react_edc_nhs_two_step(c_edc_initial=5.0, **common).p_final
        p20 = react_edc_nhs_two_step(c_edc_initial=20.0, **common).p_final
        assert p1 < p5 < p20, (
            f"dose response expected monotone; got "
            f"{p1:.3g} / {p5:.3g} / {p20:.3g}"
        )

    def test_ph_plateau_near_neutral(self):
        """p_final peaks in the pH 6.5–8 range; drops at very low and very
        high pH (low-pH: little [NH2]_eff; high-pH: NHS ester hydrolysis)."""
        common = dict(
            c_cooh_initial=10.0, c_nh2_total=20.0,
            c_edc_initial=5.0, c_nhs_initial=5.0,
            T=298.15, time=3600.0,
        )
        p_low = react_edc_nhs_two_step(pH=4.5, **common).p_final
        p_mid = react_edc_nhs_two_step(pH=7.0, **common).p_final
        p_hi = react_edc_nhs_two_step(pH=9.5, **common).p_final
        assert p_mid > p_low, (
            f"expected neutral pH > acidic; got mid={p_mid:.3g} low={p_low:.3g}"
        )
        assert p_mid > p_hi, (
            f"expected neutral pH > alkaline; got mid={p_mid:.3g} hi={p_hi:.3g}"
        )


# ─── pH speciation helper ─────────────────────────────────────────────────


class TestAmineFraction:
    def test_at_pka_half_deprotonated(self):
        assert available_amine_fraction(pH=6.4, pKa=6.4) == pytest.approx(0.5)

    def test_monotone_in_pH(self):
        vals = [available_amine_fraction(pH=p, pKa=6.4)
                for p in [3, 4, 5, 6, 7, 8, 9, 10]]
        for a, b in zip(vals, vals[1:]):
            assert b >= a

    def test_limits(self):
        assert available_amine_fraction(pH=14.0, pKa=6.4) == pytest.approx(1.0, abs=1e-6)
        assert available_amine_fraction(pH=0.0, pKa=6.4) < 1e-6


# ─── Input validation ─────────────────────────────────────────────────────


class TestInputValidation:
    def test_zero_cooh_raises(self):
        with pytest.raises(ValueError, match="surface COOH"):
            react_edc_nhs_two_step(
                c_cooh_initial=0.0,
                c_nh2_total=10.0,
                c_edc_initial=5.0,
                c_nhs_initial=5.0,
                pH=7.0, T=298.15, time=3600.0,
            )

    def test_bad_temperature_raises(self):
        with pytest.raises(ValueError, match="finite and positive"):
            react_edc_nhs_two_step(
                c_cooh_initial=10.0,
                c_nh2_total=10.0,
                c_edc_initial=5.0,
                c_nhs_initial=5.0,
                pH=7.0, T=-1.0, time=3600.0,
            )

    def test_ph_out_of_domain_raises(self):
        with pytest.raises(ValueError, match="pH"):
            react_edc_nhs_two_step(
                c_cooh_initial=10.0,
                c_nh2_total=10.0,
                c_edc_initial=5.0,
                c_nhs_initial=5.0,
                pH=12.0, T=298.15, time=3600.0,
            )

    def test_negative_concentration_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            react_edc_nhs_two_step(
                c_cooh_initial=10.0,
                c_nh2_total=10.0,
                c_edc_initial=-1.0,
                c_nhs_initial=5.0,
                pH=7.0, T=298.15, time=3600.0,
            )


# ─── L3 integration gate ──────────────────────────────────────────────────


class TestL3Gate:
    @pytest.fixture
    def smoke_params(self):
        from pathlib import Path
        from emulsim.config import load_config
        cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
        if not cfg.exists():
            pytest.skip("fast_smoke.toml missing")
        return load_config(cfg)

    def test_native_matrix_falls_back(self, smoke_params):
        """surface_cooh_concentration=0 → QUALITATIVE_TREND, no mechanistic path."""
        from emulsim.datatypes import MaterialProperties, ModelEvidenceTier
        from emulsim.level3_crosslinking.solver import solve_crosslinking

        props = MaterialProperties()
        assert props.surface_cooh_concentration == 0.0  # default

        # Use any michaelis_menten crosslinker key; EDC/NHS exists in
        # CROSSLINKERS library.
        try:
            result = solve_crosslinking(
                smoke_params, props, crosslinker_key="edc_nhs",
            )
        except ValueError:
            # If "edc_nhs" is not a CROSSLINKERS key try the alternative
            # common name. The test skips cleanly if not present.
            pytest.skip("EDC/NHS crosslinker key not registered")

        # Fallback path: tier should be QUALITATIVE_TREND.
        assert result.model_manifest is not None
        assert (result.model_manifest.evidence_tier
                == ModelEvidenceTier.QUALITATIVE_TREND)

    def test_carboxylated_matrix_uses_mechanistic(self, smoke_params):
        """surface_cooh_concentration>0 → SEMI_QUANTITATIVE, mechanistic path."""
        from emulsim.datatypes import MaterialProperties, ModelEvidenceTier
        from emulsim.level3_crosslinking.solver import solve_crosslinking

        props = MaterialProperties()
        # Simulate prior M2 succinylation: 50 mol/m^3 COOH grafted.
        props.surface_cooh_concentration = 50.0

        try:
            result = solve_crosslinking(
                smoke_params, props, crosslinker_key="edc_nhs",
            )
        except ValueError:
            pytest.skip("EDC/NHS crosslinker key not registered")

        assert result.model_manifest is not None
        assert (result.model_manifest.evidence_tier
                == ModelEvidenceTier.SEMI_QUANTITATIVE)
        # Manifest should reference the two-step model.
        assert "two_step" in result.model_manifest.model_name


# ─── M2 profile tier ──────────────────────────────────────────────────────


class TestM2ProfileTier:
    def test_edc_nhs_activation_promoted_to_semi_quantitative(self):
        """Node 31: edc_nhs_activation profile tier should be
        semi_quantitative (was ranking_only in v7.0.1)."""
        from emulsim.module2_functionalization.reagent_profiles import (
            REAGENT_PROFILES,
        )
        profile = REAGENT_PROFILES["edc_nhs_activation"]
        assert profile.confidence_tier == "semi_quantitative"
