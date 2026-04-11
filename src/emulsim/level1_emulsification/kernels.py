"""Breakage and coalescence kernel functions for the PBE solver.

Implements:
- Alopaeus et al. (2002) viscosity-corrected breakage kernel
- Coulaloglou & Tavlarides (1977) coalescence kernel
- Daughter size distribution models
"""

from __future__ import annotations

import numpy as np

from ..datatypes import BreakageModel, CoalescenceModel, KernelConfig


# ─── Breakage Kernels ─────────────────────────────────────────────────────

def breakage_rate_alopaeus(d: np.ndarray, epsilon: float, sigma: float,
                            rho_c: float, mu_d: float,
                            C1: float = 0.986, C2: float = 0.0115,
                            C3: float = 0.0, nu_c: float = None) -> np.ndarray:
    """Alopaeus et al. (2002) breakage rate with viscosity correction [1/s].

    g(d) = C₁·(ε/ν_c)^(1/2)·exp(-C₂·σ/(ρ_c·ε^(2/3)·d^(5/3))
                                  - C₃·Vi)
    where Vi = µ_d / sqrt(ρ_c·σ·d) is the dimensionless viscosity group.

    Parameters
    ----------
    d : np.ndarray
        Droplet diameters [m].
    epsilon : float
        Energy dissipation rate [m²/s³].
    sigma : float
        Interfacial tension [N/m].
    rho_c : float
        Continuous phase density [kg/m³].
    mu_d : float
        Dispersed phase dynamic viscosity [Pa·s].
    C1, C2, C3 : float
        Model constants.
    nu_c : float, optional
        Kinematic viscosity of continuous phase [m²/s].
        If None, estimated as mu_oil / rho_c with mu_oil=0.005 Pa·s.

    Returns
    -------
    np.ndarray
        Breakage rate for each diameter [1/s].
    """
    d = np.asarray(d, dtype=float)

    if nu_c is None:
        nu_c = 0.005 / rho_c  # default estimate

    if epsilon <= 0:
        return np.zeros_like(d)

    # Prefactor
    prefactor = C1 * np.sqrt(epsilon / nu_c)

    # Surface tension resistance (Coulaloglou-Tavlarides term)
    with np.errstate(divide='ignore', invalid='ignore'):
        exp_arg1 = -C2 * sigma / (rho_c * epsilon**(2.0/3.0) * d**(5.0/3.0))

    # Viscous resistance (Alopaeus extension)
    # Vi = mu_d / sqrt(rho_c * sigma * d) is the dimensionless viscosity group
    # that measures the ratio of viscous to capillary resistance to deformation.
    exp_arg2 = np.zeros_like(d)
    if C3 > 0 and mu_d > 0:
        with np.errstate(divide='ignore', invalid='ignore'):
            Vi = mu_d / np.sqrt(rho_c * sigma * np.maximum(d, 1e-15))
            # Cap Vi to prevent exp(-C3*Vi) from killing breakage at small d.
            # Without the cap, sub-micron droplets can have Vi > 50, making
            # exp(-C3*Vi) ~ 0 and creating a nonmonotonic RPM->d32 feedback
            # loop where breakage dies but coalescence continues (see F1 audit).
            Vi = np.minimum(Vi, 100.0)
            exp_arg2 = -C3 * Vi

    # Clip total exponent to prevent extreme underflow (g < 1e-100 is effectively 0)
    exp_total = np.maximum(exp_arg1 + exp_arg2, -200.0)

    g = prefactor * np.exp(exp_total)

    # Ensure non-negative and handle numerical issues
    return np.where(np.isfinite(g), np.maximum(g, 0.0), 0.0)


def breakage_rate_coulaloglou(d: np.ndarray, epsilon: float, sigma: float,
                               rho_c: float,
                               C1: float = 0.00481, C2: float = 0.08) -> np.ndarray:
    """Coulaloglou & Tavlarides (1977) breakage rate [1/s].

    g(d) = C₁·ε^(1/3)·d^(-2/3)·exp(-C₂·σ/(ρ_c·ε^(2/3)·d^(5/3)))

    Classical model without viscous correction.
    """
    d = np.asarray(d, dtype=float)

    if epsilon <= 0:
        return np.zeros_like(d)

    prefactor = C1 * epsilon**(1.0/3.0) * d**(-2.0/3.0)

    with np.errstate(divide='ignore', invalid='ignore'):
        exp_arg = -C2 * sigma / (rho_c * epsilon**(2.0/3.0) * d**(5.0/3.0))

    g = prefactor * np.exp(exp_arg)
    return np.where(np.isfinite(g), np.maximum(g, 0.0), 0.0)


# ─── Daughter Size Distribution ───────────────────────────────────────────

def daughter_uniform_binary(v_daughter: float, v_parent: float) -> float:
    """Uniform binary breakage: parent breaks into two equal fragments.

    β(v|v') = 2·δ(v - v'/2)

    In discretized form, this is handled by the fixed-pivot redistribution.
    This function returns the daughter volume for binary breakage.
    """
    return v_parent / 2.0


def daughter_beta_distribution(v_ratio: np.ndarray, alpha: float = 2.0,
                                beta_param: float = 2.0) -> np.ndarray:
    """Beta distribution for daughter size fraction.

    P(f) = f^(α-1)·(1-f)^(β-1) / B(α,β)
    where f = v_daughter / v_parent.

    Parameters
    ----------
    v_ratio : np.ndarray
        Volume fraction f = v_daughter / v_parent.
    alpha, beta_param : float
        Beta distribution parameters. alpha=beta=2 gives symmetric bell shape.

    Returns
    -------
    np.ndarray
        Probability density for each fraction.
    """
    from scipy.stats import beta as beta_dist
    return beta_dist.pdf(v_ratio, alpha, beta_param)


# ─── Coalescence Kernels ──────────────────────────────────────────────────

def coalescence_rate_ct(d_i: np.ndarray, d_j: np.ndarray,
                        epsilon: float, sigma: float,
                        rho_c: float, mu_c: float,
                        phi_d: float = 0.05,
                        C4: float = 2.17e-4, C5: float = 2.28e13) -> np.ndarray:
    """Coulaloglou & Tavlarides (1977) coalescence rate [m³/s].

    q(d_i, d_j) = h(d_i, d_j) · λ(d_i, d_j)

    h = C₄·ε^(1/3)·(d_i² + d_j²)·(d_i^(2/3) + d_j^(2/3))^(1/2) / (1 + φ_d)
    λ = exp(-C₅·µ_c·ρ_c·ε/σ²·(d_i·d_j/(d_i+d_j))⁴)

    Parameters
    ----------
    d_i, d_j : np.ndarray
        Droplet diameters [m]. Can be 1D arrays; will be broadcast.
    epsilon : float
        Energy dissipation rate [m²/s³].
    sigma : float
        Interfacial tension [N/m].
    rho_c : float
        Continuous phase density [kg/m³].
    mu_c : float
        Continuous phase dynamic viscosity [Pa·s].
    phi_d : float
        Dispersed phase volume fraction.
    C4, C5 : float
        Model constants.

    Returns
    -------
    np.ndarray
        Coalescence rate [m³/s].
    """
    d_i = np.asarray(d_i, dtype=float)
    d_j = np.asarray(d_j, dtype=float)

    if epsilon <= 0:
        return np.zeros(np.broadcast_shapes(d_i.shape, d_j.shape))

    # Collision frequency following Coulaloglou & Tavlarides (1977), Eq. 17.
    # h(d_i, d_j) = C4 * eps^(1/3) * (d_i^2 + d_j^2) * (d_i^(2/3) + d_j^(2/3))^(1/2) / (1+phi_d)
    # The (d_i^2 + d_j^2) term is the collision cross-section and
    # (d_i^(2/3) + d_j^(2/3))^(1/2) is the relative turbulent velocity.
    h = (C4 * epsilon**(1.0/3.0) / (1.0 + phi_d)
         * (d_i**2 + d_j**2)
         * np.sqrt(d_i**(2.0/3.0) + d_j**(2.0/3.0)))

    # Film drainage efficiency
    d_sum = d_i + d_j
    d_harmonic = np.where(d_sum > 0, d_i * d_j / d_sum, 0.0)

    exp_arg = -C5 * mu_c * rho_c * epsilon / sigma**2 * d_harmonic**4
    # Clip to prevent underflow
    exp_arg = np.maximum(exp_arg, -500.0)
    lam = np.exp(exp_arg)

    q = h * lam
    return np.where(np.isfinite(q), np.maximum(q, 0.0), 0.0)


# ─── Kolmogorov-Hinze Prediction ──────────────────────────────────────────

def hinze_dmax(epsilon: float, sigma: float, rho_c: float,
               C1: float = 0.725) -> float:
    """Classical Kolmogorov-Hinze maximum stable droplet size [m].

    d_max = C₁·(σ/ρ_c)^(3/5)·ε^(-2/5)

    Valid for inertial sub-range (d >> η_K).
    """
    if epsilon <= 0:
        return np.inf
    return C1 * (sigma / rho_c)**0.6 * epsilon**(-0.4)


def hinze_dmax_viscous(epsilon: float, sigma: float, rho_c: float,
                       mu_d: float, C1: float = 0.054,
                       C3: float = 1.38) -> float:
    """Calabrese et al. (1986) modified Hinze with viscous correction [m].

    Solves iteratively:
    d_max = C₁·(σ/ρ_c)^(3/5)·ε^(-2/5)·[1 + C₃·Vi]^(3/5)
    where Vi = µ_d / sqrt(ρ_c·σ·d_max)

    Parameters
    ----------
    epsilon : float
        Energy dissipation rate [m²/s³].
    sigma : float
        Interfacial tension [N/m].
    rho_c : float
        Continuous phase density [kg/m³].
    mu_d : float
        Dispersed phase viscosity [Pa·s].
    """
    if epsilon <= 0:
        return np.inf

    # Initial guess from classical Hinze
    d = C1 * (sigma / rho_c)**0.6 * epsilon**(-0.4)

    # Fixed-point iteration
    for _ in range(50):
        Vi = mu_d / np.sqrt(rho_c * sigma * max(d, 1e-15))
        d_new = C1 * (sigma / rho_c)**0.6 * epsilon**(-0.4) * (1.0 + C3 * Vi)**0.6
        if abs(d_new - d) / max(d, 1e-15) < 1e-8:
            break
        d = d_new

    return d


# ─── Dispatch Functions ──────────────────────────────────────────────────

def breakage_rate_dispatch(
    d: np.ndarray,
    epsilon: float,
    sigma: float,
    rho_c: float,
    mu_d: float,
    config: KernelConfig,
    nu_c: float | None = None,
) -> np.ndarray:
    """Select and evaluate the breakage rate kernel based on KernelConfig.

    Parameters
    ----------
    d : np.ndarray
        Droplet diameters [m].
    epsilon : float
        Energy dissipation rate [m²/s³].
    sigma : float
        Interfacial tension [N/m].
    rho_c : float
        Continuous phase density [kg/m³].
    mu_d : float
        Dispersed phase dynamic viscosity [Pa·s].
    config : KernelConfig
        Kernel configuration selecting the model and its constants.
    nu_c : float, optional
        Kinematic viscosity of continuous phase [m²/s].
        Required for Alopaeus model; if None, computed as mu_d/(rho_c) fallback.

    Returns
    -------
    np.ndarray
        Breakage rate for each diameter [1/s].
    """
    if getattr(config.breakage_model, 'value', config.breakage_model) == "alopaeus":
        kwargs = dict(
            C1=config.breakage_C1,
            C2=config.breakage_C2,
            C3=config.breakage_C3,
        )
        if nu_c is not None:
            kwargs["nu_c"] = nu_c
        return breakage_rate_alopaeus(
            d, epsilon, sigma, rho_c, mu_d, **kwargs,
        )
    elif getattr(config.breakage_model, 'value', config.breakage_model) == "coulaloglou_tavlarides":
        return breakage_rate_coulaloglou(
            d, epsilon, sigma, rho_c,
            C1=config.breakage_C1,
            C2=config.breakage_C2,
        )
    else:
        raise ValueError(f"Unknown breakage model: {config.breakage_model}")


def coalescence_rate_dispatch(
    d: np.ndarray,
    epsilon: float,
    sigma: float,
    rho_c: float,
    config: KernelConfig,
    phi_d: float = 0.05,
    mu_c: float = None,
) -> np.ndarray:
    """Build the coalescence rate matrix using the configured model.

    Parameters
    ----------
    d : np.ndarray
        Droplet diameters [m], 1-D array of length N.
    epsilon : float
        Energy dissipation rate [m²/s³].
    sigma : float
        Interfacial tension [N/m].
    rho_c : float
        Continuous phase density [kg/m³].
    config : KernelConfig
        Kernel configuration selecting the model and its constants.
    phi_d : float
        Dispersed phase volume fraction (used for crowding correction).
    mu_c : float, optional
        Continuous phase dynamic viscosity [Pa·s].
        If None, estimated as 0.001 Pa·s (water-like).

    Returns
    -------
    np.ndarray
        Coalescence rate matrix of shape (N, N) [m³/s].
    """
    # ASSUMPTION: default mu_c = 0.001 Pa·s (water at ~20 °C)
    if mu_c is None:
        mu_c = 0.001

    d = np.asarray(d, dtype=float)
    N = d.size

    # Build pairwise rate matrix using the underlying CT coalescence kernel.
    # ASSUMPTION: CoalescenceModel is always CT for now; extensible via elif.
    d_i = d[:, np.newaxis]  # (N, 1)
    d_j = d[np.newaxis, :]  # (1, N)

    Q = coalescence_rate_ct(
        d_i, d_j, epsilon, sigma, rho_c, mu_c,
        phi_d=phi_d,
        C4=config.coalescence_C4,
        C5=config.coalescence_C5,
    )

    # Concentrated-emulsion correction: dampen coalescence at high phi_d.
    # The CT kernel already has a 1/(1+phi_d) in the collision frequency.
    # The additional correction here accounts for crowding-induced drainage
    # inhibition not captured by the original CT model.
    # Total coalescence ~ 1/(1+phi_d)^(1+n) where n = coalescence_exponent.
    # At phi_d=0.40 with n=2: total damping = 1/(1.4)^(1+2) = 0.36
    if config.phi_d_correction and config.coalescence_exponent > 0:
        Q = Q / (1.0 + phi_d) ** config.coalescence_exponent

    return Q
