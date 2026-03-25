"""Trust gates -- assess reliability of simulation predictions."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .datatypes import FullResult, SimulationParameters, MaterialProperties


@dataclass
class TrustAssessment:
    """Assessment of prediction reliability."""
    trustworthy: bool
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)

    @property
    def level(self) -> str:
        if self.blockers:
            return "UNRELIABLE"
        if self.warnings:
            return "CAUTION"
        return "TRUSTWORTHY"

    def summary(self) -> str:
        lines = [f"Trust Level: {self.level}"]
        for b in self.blockers:
            lines.append(f"  [BLOCKER] {b}")
        for w in self.warnings:
            lines.append(f"  [WARNING] {w}")
        if not self.blockers and not self.warnings:
            lines.append("  All checks passed.")
        return "\n".join(lines)


def assess_trust(result: FullResult, params: SimulationParameters,
                 props: MaterialProperties) -> TrustAssessment:
    """Evaluate whether the simulation results should be trusted.

    Checks for conditions where the model assumptions break down.
    """
    warnings = []
    blockers = []

    e = result.emulsification
    g = result.gelation
    x = result.crosslinking
    m = result.mechanical

    # 1. Droplet size outside model validity
    if e.d32 < 0.5e-6:
        blockers.append(
            f"d32={e.d32*1e6:.2f} um is below 0.5 um -- sub-Kolmogorov regime, "
            "PBE breakage kernels not validated at this scale"
        )
    if e.d32 > 100e-6:
        warnings.append(
            f"d32={e.d32*1e6:.0f} um -- very large droplets, may not form stable microspheres"
        )

    # 2. PBE convergence
    if not e.converged:
        warnings.append("PBE solver did not reach steady state -- increase emulsification time")

    # 3. Size distribution too broad
    if e.span > 2.0:
        warnings.append(
            f"Span={e.span:.1f} -- very polydisperse, not suitable for chromatographic media"
        )

    # 4. Pore size outside useful range for chromatography
    if g.pore_size_mean < 20e-9:
        warnings.append(
            f"Pore size={g.pore_size_mean*1e9:.0f} nm -- too small for protein access"
        )
    if g.pore_size_mean > 1000e-9:
        warnings.append(
            f"Pore size={g.pore_size_mean*1e9:.0f} nm -- too large for SEC resolution"
        )

    # 5. Gelation incomplete
    if g.alpha_final < 0.8:
        blockers.append(
            f"Gelation only {g.alpha_final:.0%} complete -- microspheres may not form properly. "
            "Check cooling rate and T_gel."
        )

    # 6. Crosslinking stoichiometry-limited
    if x.p_final < 0.01:
        warnings.append(
            f"Crosslinking fraction={x.p_final:.1%} -- very low, increase genipin concentration"
        )

    # 7. Genipin is limiting reagent
    NH2_total = params.formulation.c_chitosan * 1000 * props.DDA / props.M_GlcN
    if NH2_total > 0 and params.formulation.c_genipin / NH2_total < 0.05:
        warnings.append(
            f"Genipin/NH2 ratio = {params.formulation.c_genipin/NH2_total:.3f} -- "
            "genipin-limited, increasing crosslinking time will not help"
        )

    # 8. Mechanical properties unreasonable
    if m.G_DN < 100:
        blockers.append(
            f"G_DN={m.G_DN:.0f} Pa -- too soft for column packing"
        )
    if m.G_DN > 1e6:
        warnings.append(
            f"G_DN={m.G_DN/1000:.0f} kPa -- unusually stiff, check agarose concentration"
        )

    # 9. Viscosity ratio extreme
    if props.mu_d > 0 and props.mu_oil > 0:
        ratio = props.mu_d / props.mu_oil
        if ratio > 100:
            warnings.append(
                f"Viscosity ratio mu_d/mu_c={ratio:.0f} -- extreme, breakage predictions uncertain. "
                "Consider enabling C3 viscous correction."
            )

    # 10. Empirical pore model outside calibration range
    c_agar_pct = params.formulation.c_agarose / 10.0
    if c_agar_pct < 1.0 or c_agar_pct > 8.0:
        warnings.append(
            f"Agarose={c_agar_pct:.1f}% -- outside empirical pore model calibration range (1-8%)"
        )

    trustworthy = len(blockers) == 0

    return TrustAssessment(
        trustworthy=trustworthy,
        warnings=warnings,
        blockers=blockers,
    )
