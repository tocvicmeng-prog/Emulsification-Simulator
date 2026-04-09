"""Population Balance Equation solver using the Fixed-Pivot method.

Implements Kumar & Ramkrishna (1996) fixed-pivot technique with:
- Logarithmically spaced bins
- Alopaeus breakage kernel with viscosity correction
- Coulaloglou-Tavlarides coalescence kernel
- BDF time integration via scipy.integrate.solve_ivp

Coalescence birth is precomputed via a target-bin mapping for O(n²) RHS.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..datatypes import (
    EmulsificationResult,
    MaterialProperties,
    SimulationParameters,
)
import logging

from .energy import (
    emulsion_density,
    gap_shear_rate,
    max_dissipation,
    average_dissipation,
    power_draw,
    impeller_reynolds_number,
    power_number_corrected,
)
from .kernels import (
    breakage_rate_alopaeus,
    coalescence_rate_ct,
    breakage_rate_dispatch,
    coalescence_rate_dispatch,
)
from .thermal import temperature_profile

logger = logging.getLogger(__name__)


class PBESolver:
    """Fixed-Pivot Population Balance Equation solver.

    Solves the spatially homogeneous (0D) PBE for droplet size distribution
    evolution under simultaneous breakage and coalescence.
    """

    def __init__(self, n_bins: int = 50, d_min: float = 0.1e-6,
                 d_max: float = 500e-6):
        self.n_bins = n_bins

        # Logarithmically spaced bin edges
        self.d_edges = np.logspace(np.log10(d_min), np.log10(d_max), n_bins + 1)

        # Pivot diameters (geometric mean of edges)
        self.d_pivots = np.sqrt(self.d_edges[:-1] * self.d_edges[1:])

        # Volume pivots
        self.v_pivots = np.pi / 6.0 * self.d_pivots**3

        # Bin widths (in diameter space)
        self.d_widths = self.d_edges[1:] - self.d_edges[:-1]

        # Precompute coalescence target-bin mapping
        self._coal_target = self._build_coalescence_map()

    def _build_coalescence_map(self) -> list:
        """Precompute which bin each (j,k) coalescence product lands in.

        Returns a list of (j, k, target_bin, weight) tuples for all pairs
        where the coalesced volume falls within the grid. Uses fixed-pivot
        redistribution to the two nearest bins.
        """
        n = self.n_bins
        entries = []
        for j in range(n):
            for k in range(j, n):
                v_sum = self.v_pivots[j] + self.v_pivots[k]
                d_sum = (6.0 * v_sum / np.pi) ** (1.0 / 3.0)

                # Skip if product is below grid
                if d_sum < self.d_edges[0]:
                    continue

                # Clamp oversized products to the last bin (conservative)
                if d_sum > self.d_edges[-1]:
                    entries.append((j, k, n - 1, 1.0))
                    continue

                # Find target bin via searchsorted on pivot diameters
                if d_sum <= self.d_pivots[0]:
                    entries.append((j, k, 0, 1.0))
                elif d_sum >= self.d_pivots[-1]:
                    entries.append((j, k, n - 1, 1.0))
                else:
                    idx = np.searchsorted(self.d_pivots, d_sum) - 1
                    idx = max(0, min(idx, n - 2))
                    i_lo = idx
                    i_hi = idx + 1
                    v_lo = self.v_pivots[i_lo]
                    v_hi = self.v_pivots[i_hi]
                    if v_hi > v_lo:
                        w_lo = (v_hi - v_sum) / (v_hi - v_lo)
                        w_hi = 1.0 - w_lo
                    else:
                        w_lo = 1.0
                        w_hi = 0.0
                    if w_lo > 1e-15:
                        entries.append((j, k, i_lo, w_lo))
                    if w_hi > 1e-15:
                        entries.append((j, k, i_hi, w_hi))

        return entries

    def _build_breakage_matrix(self, g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Build breakage birth matrix using fixed-pivot redistribution.

        Binary equal breakage is used (each daughter = v_parent/2).
        A beta-distribution daughter size is available in kernels.daughter_beta_distribution
        but requires additional redistribution logic. Binary breakage is the standard
        Kumar-Ramkrishna default and matches most published PBE studies.
        """
        n = self.n_bins
        birth_matrix = np.zeros((n, n))

        for j in range(n):
            if g[j] <= 0:
                continue

            v_daughter = self.v_pivots[j] / 2.0
            d_daughter = (6.0 * v_daughter / np.pi) ** (1.0 / 3.0)

            if d_daughter <= self.d_pivots[0]:
                birth_matrix[0, j] += 2.0 * g[j]
            elif d_daughter >= self.d_pivots[-1]:
                birth_matrix[-1, j] += 2.0 * g[j]
            else:
                idx = np.searchsorted(self.d_pivots, d_daughter) - 1
                idx = max(0, min(idx, n - 2))
                i_low, i_high = idx, idx + 1

                v_low = self.v_pivots[i_low]
                v_high = self.v_pivots[i_high]

                if v_high > v_low:
                    frac_low = (v_high - v_daughter) / (v_high - v_low)
                    frac_high = 1.0 - frac_low
                else:
                    frac_low, frac_high = 1.0, 0.0

                birth_matrix[i_low, j] += 2.0 * frac_low * g[j]
                birth_matrix[i_high, j] += 2.0 * frac_high * g[j]

        return birth_matrix, g.copy()

    def _compute_rhs(self, t: float, N: np.ndarray,
                     birth_matrix: np.ndarray, death_rate: np.ndarray,
                     Q: np.ndarray) -> np.ndarray:
        """Compute dN/dt — vectorized for performance."""
        N = np.maximum(N, 0.0)
        n = self.n_bins
        dNdt = np.zeros(n)

        # Breakage: birth - death  (fully vectorized)
        dNdt += birth_matrix @ N
        dNdt -= death_rate * N

        # Coalescence death: vectorized  (O(n) per bin via dot product)
        QN = Q @ N  # shape (n,)
        dNdt -= N * QN

        # Coalescence birth: use precomputed map (O(n²) total entries)
        for j, k, target, weight in self._coal_target:
            sym = 1.0 if j == k else 2.0
            dNdt[target] += 0.5 * sym * weight * Q[j, k] * N[j] * N[k]

        return dNdt

    def _initial_distribution(self, phi_d: float,
                              d32_premix: float = 100e-6,
                              sigma_premix: float = 0.5) -> np.ndarray:
        """Log-normal initial (premix) distribution normalised to phi_d."""
        log_d = np.log(self.d_pivots)
        log_d0 = np.log(d32_premix)

        pdf = np.exp(-0.5 * ((log_d - log_d0) / sigma_premix) ** 2)
        pdf /= (self.d_pivots * sigma_premix * np.sqrt(2 * np.pi))

        # N_i = n(d_i) * Δd_i  (total count per bin)
        N = pdf * self.d_widths
        total_vol = np.sum(N * self.v_pivots)
        if total_vol > 0:
            N *= phi_d / total_vol

        return N

    def _compute_d_mode(self, N: np.ndarray) -> float:
        """Modal diameter from volume-weighted distribution."""
        vol_density = N * self.v_pivots / self.d_widths  # volume density per µm
        if vol_density.max() <= 0:
            return 0.0
        return self.d_pivots[np.argmax(vol_density)]

    def solve(self, params: SimulationParameters,
              props: MaterialProperties,
              phi_d: float | None = None) -> EmulsificationResult:
        """Solve the PBE for given process conditions.

        Parameters
        ----------
        phi_d : float, optional
            Dispersed-phase volume fraction override.  If ``None``,
            defaults to ``formulation.phi_d_from_volumes`` (stirred-vessel)
            or ``formulation.phi_d`` (legacy).
        """
        if params.emulsification.mode == "stirred_vessel":
            return self.solve_stirred_vessel(params, props, phi_d=phi_d)

        # Legacy mode: resolve default
        if phi_d is None:
            phi_d = params.formulation.phi_d

        rpm = params.emulsification.rpm
        t_emul = params.emulsification.t_emulsification
        mixer = params.emulsification.mixer

        # Energy dissipation
        rho_emul = emulsion_density(props.rho_oil, props.rho_aq, phi_d)
        epsilon_max = max_dissipation(mixer, rpm, rho_emul)
        epsilon_avg = average_dissipation(mixer, rpm, rho_emul)

        # Note: equilibrium interfacial tension (props.sigma) is used rather than
        # the dynamic σ(t) from properties/interfacial.py:dynamic_interfacial_tension().
        # At steady state (t >> τ_ads), σ_dynamic → σ_equilibrium.  The adsorption
        # timescale τ_ads = Γ_∞² / (D · c²) ~ O(ms) for typical Span-80 concentrations,
        # which is much shorter than the emulsification time (~minutes).  The dynamic
        # model is relevant only for the initial transient (first few ms after interface
        # creation) and does not affect the steady-state size distribution solved here.

        # Dispersed phase viscosity:
        # Breakage uses ZERO-SHEAR viscosity (resistance to bulk deformation).
        # Coalescence uses mu_oil (continuous phase) for film drainage, which
        # does not require shear correction for Newtonian oils.

        # Breakage rates (use zero-shear mu_d for viscous resistance Vi)
        # epsilon_max is used: breakage occurs in the high-shear rotor-stator gap
        nu_c = props.mu_oil / props.rho_oil
        g = breakage_rate_alopaeus(
            self.d_pivots, epsilon_max, props.sigma, props.rho_oil,
            props.mu_d, C3=props.breakage_C3, nu_c=nu_c,
        )
        birth_matrix, death_rate = self._build_breakage_matrix(g)

        # Coalescence rate matrix (vectorised construction)
        # epsilon_avg is used: coalescence occurs in the bulk flow
        di_grid, dj_grid = np.meshgrid(self.d_pivots, self.d_pivots, indexing='ij')
        Q = coalescence_rate_ct(
            di_grid, dj_grid,
            epsilon_avg, props.sigma, props.rho_oil,
            props.mu_oil, phi_d=phi_d,
        )

        # Initial distribution
        N0 = self._initial_distribution(phi_d)

        # Integrate
        def rhs(t, N):
            return self._compute_rhs(t, N, birth_matrix, death_rate, Q)

        t_span = (0.0, t_emul)
        t_eval = np.linspace(0, t_emul, 101)

        sol = solve_ivp(
            rhs, t_span, N0,
            method='BDF',
            rtol=params.solver.l1_rtol,
            atol=params.solver.l1_atol,
            t_eval=t_eval,
            max_step=t_emul / 10,
        )

        if not sol.success:
            logger.warning("PBE solver did not converge: %s", sol.message)

        N_final = np.maximum(sol.y[:, -1], 0.0)
        d32, d43, d10, d50, d90, span = self._compute_statistics(N_final)

        # Convergence check — d32 stable over last 10 % of time
        n_check = max(1, len(sol.t) // 10)
        if sol.y.shape[1] > n_check:
            d32_late = [self._sauter_mean(np.maximum(sol.y[:, k], 0.0))
                        for k in range(-n_check, 0)]
            converged = (max(d32_late) - min(d32_late)) / max(d32, 1e-15) < 0.01
        else:
            converged = False

        total_vol = np.sum(N_final * self.v_pivots)

        # Convert from total count per bin back to number density for output
        n_d_output = N_final / self.d_widths
        n_d_history = sol.y.T / self.d_widths[np.newaxis, :]

        d_mode = self._compute_d_mode(N_final)

        return EmulsificationResult(
            d_bins=self.d_pivots.copy(),
            n_d=n_d_output,
            d32=d32, d43=d43, d10=d10, d50=d50, d90=d90, span=span,
            total_volume_fraction=total_vol,
            converged=converged,
            d_mode=d_mode,
            t_history=sol.t,
            n_d_history=n_d_history,
        )

    # ── Stirred-vessel solver ────────────────────────────────────────────

    def solve_stirred_vessel(
        self,
        params: SimulationParameters,
        props: MaterialProperties,
        phi_d: float | None = None,
    ) -> EmulsificationResult:
        """Solve the PBE for stirred-vessel mode with time-dependent temperature.

        Key differences from legacy solve():
        - Uses vessel/stirrer/heating/kernels from params.emulsification
        - Time-dependent T(t) -> time-varying mu, sigma, epsilon
        - Kernel dispatch via KernelConfig
        - Broader premix (d32=500 µm, sigma=0.8)
        - phi_d from caller or formulation volumes
        """
        emul = params.emulsification
        vessel = emul.vessel
        stirrer = emul.stirrer
        heating = emul.heating
        kernels = emul.kernels
        rpm = emul.rpm
        t_emul = emul.t_emulsification
        tank_vol = vessel.working_volume

        # Dispersed-phase volume fraction: honour caller override, else from volumes
        if phi_d is None:
            phi_d = params.formulation.phi_d_from_volumes

        # Emulsion density (constant — density weakly T-dependent)
        # ASSUMPTION: rho_emul is constant over the temperature range
        rho_emul = emulsion_density(props.rho_oil, props.rho_aq, phi_d)

        # ── Temperature profile ─────────────────────────────────────────
        t_eval = np.linspace(0.0, t_emul, 101)
        T_profile = temperature_profile(heating, vessel, t_eval)
        T_initial = T_profile[0]
        T_final = T_profile[-1]

        # ── Arrhenius viscosity constants ────────────────────────────────
        # ASSUMPTION: E_a_mu / R ≈ 2500 K for paraffin oil (literature range
        # 2000-3000 K for mineral oils).
        EA_R = 2500.0       # [K]
        # props.mu_oil and props.sigma are already interpolated to
        # params.formulation.T_oil by PropertyDatabase.update_for_conditions().
        # Use T_oil as the reference temperature for Arrhenius/linear correction.
        T_REF = params.formulation.T_oil  # [K] — the temperature props were computed at
        mu_oil_ref = props.mu_oil         # [Pa·s] at T_REF
        sigma_ref = props.sigma           # [N/m] at T_REF

        def _mu_oil_at_T(T: float) -> float:
            """Arrhenius viscosity: mu(T) = mu_ref * exp(Ea_R*(1/T - 1/T_ref))."""
            return mu_oil_ref * np.exp(EA_R * (1.0 / T - 1.0 / T_REF))

        def _sigma_at_T(T: float) -> float:
            """Linear T-dependence of interfacial tension relative to T_REF.

            ASSUMPTION: d(sigma)/dT ≈ -0.0001 N/(m·K), typical for oil-water
            interfaces with surfactant.
            """
            return max(sigma_ref - 0.0001 * (T - T_REF), 1e-4)  # floor prevents negative IFT

        def _epsilon_avg_at_T(T: float) -> float:
            """Recompute average dissipation with Re-dependent Np correction."""
            mu_oil_T = _mu_oil_at_T(T)
            Re = impeller_reynolds_number(
                rpm, stirrer.impeller_diameter, rho_emul, mu_oil_T,
            )
            Np_corr = power_number_corrected(stirrer, Re)
            N_rps = rpm / 60.0
            D = stirrer.impeller_diameter
            P = Np_corr * rho_emul * N_rps**3 * D**5
            return P / (rho_emul * tank_vol)

        # ── Precompute kernels at temperature checkpoints ────────────────
        # ASSUMPTION: 15 checkpoints provide sufficient resolution for
        # linear interpolation of breakage/coalescence rates vs temperature.
        N_CHECKPOINTS = 15
        if abs(T_final - T_initial) < 0.1:
            # Isothermal case — single checkpoint
            T_checkpoints = np.array([T_initial])
        else:
            T_checkpoints = np.linspace(
                min(T_initial, T_final),
                max(T_initial, T_final),
                N_CHECKPOINTS,
            )

        # Storage: breakage vectors and coalescence matrices at each checkpoint
        g_cache = np.zeros((len(T_checkpoints), self.n_bins))
        Q_cache = np.zeros((len(T_checkpoints), self.n_bins, self.n_bins))

        # ASSUMPTION: dispersed-phase (aqueous polymer) viscosity follows
        # Arrhenius with E_a_d/R ~ 2000 K (weaker T-dependence than oil).
        EA_R_D = 2000.0  # [K]
        mu_d_ref = props.mu_d  # at T_REF

        def _mu_d_at_T(T: float) -> float:
            return mu_d_ref * np.exp(EA_R_D * (1.0 / T - 1.0 / T_REF))

        for idx, T_ck in enumerate(T_checkpoints):
            sigma_T = _sigma_at_T(T_ck)
            eps_avg = _epsilon_avg_at_T(T_ck)
            eps_max = stirrer.dissipation_ratio * eps_avg
            mu_oil_T = _mu_oil_at_T(T_ck)
            mu_d_T = _mu_d_at_T(T_ck)

            # Breakage (use eps_max — breakage in high-shear impeller zone)
            nu_c_T = mu_oil_T / props.rho_oil  # kinematic viscosity at T_ck
            g_cache[idx, :] = breakage_rate_dispatch(
                self.d_pivots, eps_max, sigma_T, props.rho_oil,
                mu_d_T, kernels, nu_c=nu_c_T,
            )

            # Coalescence (use eps_avg — coalescence in bulk flow)
            Q_cache[idx, :, :] = coalescence_rate_dispatch(
                self.d_pivots, eps_avg, sigma_T, props.rho_oil,
                kernels, phi_d=phi_d, mu_c=mu_oil_T,
            )

        def _interpolate_kernels(T: float):
            """Interpolate breakage vector and coalescence matrix at temperature T."""
            if len(T_checkpoints) == 1:
                return g_cache[0], Q_cache[0]

            # Clamp T to checkpoint range
            T_lo = T_checkpoints[0]
            T_hi = T_checkpoints[-1]
            T_clamped = np.clip(T, T_lo, T_hi)

            # Find bracketing indices
            idx = np.searchsorted(T_checkpoints, T_clamped) - 1
            idx = max(0, min(idx, len(T_checkpoints) - 2))

            T_a = T_checkpoints[idx]
            T_b = T_checkpoints[idx + 1]
            if T_b > T_a:
                frac = (T_clamped - T_a) / (T_b - T_a)
            else:
                frac = 0.0

            g_interp = g_cache[idx] * (1.0 - frac) + g_cache[idx + 1] * frac
            Q_interp = Q_cache[idx] * (1.0 - frac) + Q_cache[idx + 1] * frac
            return g_interp, Q_interp

        # ── Initial distribution (broader premix for stirred vessel) ─────
        # ASSUMPTION: stirred-vessel premix is coarser than rotor-stator:
        # d32_premix = 500 µm, sigma_premix = 0.8 (broader log-normal).
        d32_premix_sv = 500e-6
        if self.d_edges[-1] < d32_premix_sv * 2.0:
            logger.warning(
                "l1_d_max (%.0f um) is too small for stirred-vessel premix "
                "(d32=%.0f um). Clamping premix to grid maximum.",
                self.d_edges[-1] * 1e6, d32_premix_sv * 1e6,
            )
            d32_premix_sv = self.d_edges[-1] / 3.0  # fit within grid
        N0 = self._initial_distribution(
            phi_d, d32_premix=d32_premix_sv, sigma_premix=0.8,
        )

        # ── RHS with temperature-dependent kernels ───────────────────────
        def rhs(t, N):
            T_now = float(np.interp(t, t_eval, T_profile))
            g, Q = _interpolate_kernels(T_now)
            birth_matrix, death_rate = self._build_breakage_matrix(g)
            return self._compute_rhs(t, N, birth_matrix, death_rate, Q)

        # ── Integrate ────────────────────────────────────────────────────
        t_span = (0.0, t_emul)

        sol = solve_ivp(
            rhs, t_span, N0,
            method='BDF',
            rtol=params.solver.l1_rtol,
            atol=params.solver.l1_atol,
            t_eval=t_eval,
            max_step=t_emul / 10,
        )

        if not sol.success:
            logger.warning("PBE solver (stirred-vessel) did not converge: %s",
                           sol.message)

        N_final = np.maximum(sol.y[:, -1], 0.0)
        d32, d43, d10, d50, d90, span = self._compute_statistics(N_final)
        d_mode = self._compute_d_mode(N_final)

        # Breakage regime check: warn if d_mode/eta_K < 5
        # (CT kernel assumes inertial sub-range, valid for d >> eta_K)
        T_final_actual = float(np.interp(t_emul, t_eval, T_profile))
        mu_final = _mu_oil_at_T(T_final_actual)
        nu_c_final = mu_final / props.rho_oil
        eps_final = _epsilon_avg_at_T(T_final_actual) * stirrer.dissipation_ratio
        if eps_final > 0:
            from .energy import kolmogorov_length_scale
            eta_K = kolmogorov_length_scale(eps_final, nu_c_final)
            if d_mode > 0 and eta_K > 0:
                regime_ratio = d_mode / eta_K
                if regime_ratio < 5.0:
                    logger.warning(
                        "Breakage regime marginal: d_mode/eta_K = %.1f (< 5). "
                        "CT kernel assumes inertial sub-range (d >> eta_K). "
                        "Consider using Alopaeus kernel for better accuracy at "
                        "d_mode=%.0f um, eta_K=%.0f um.",
                        regime_ratio, d_mode * 1e6, eta_K * 1e6,
                    )

        # Convergence check — d32 stable over last 20% of time (relaxed for
        # stirred-vessel: larger droplets have noisier late-time dynamics)
        n_check = max(1, len(sol.t) // 5)   # 20% of history
        conv_tol = 0.05                       # 5% tolerance
        if sol.y.shape[1] > n_check:
            d32_late = [self._sauter_mean(np.maximum(sol.y[:, k], 0.0))
                        for k in range(-n_check, 0)]
            converged = (max(d32_late) - min(d32_late)) / max(d32, 1e-15) < conv_tol
        else:
            converged = False

        total_vol = np.sum(N_final * self.v_pivots)

        # Convert from total count per bin back to number density for output
        n_d_output = N_final / self.d_widths
        n_d_history = sol.y.T / self.d_widths[np.newaxis, :]

        return EmulsificationResult(
            d_bins=self.d_pivots.copy(),
            n_d=n_d_output,
            d32=d32, d43=d43, d10=d10, d50=d50, d90=d90, span=span,
            total_volume_fraction=total_vol,
            converged=converged,
            d_mode=d_mode,
            t_history=sol.t,
            n_d_history=n_d_history,
        )

    # ── Statistics helpers ────────────────────────────────────────────────

    def _sauter_mean(self, N: np.ndarray) -> float:
        num = np.sum(N * self.d_pivots**3)
        den = np.sum(N * self.d_pivots**2)
        return num / den if den > 0 else 0.0

    def _compute_statistics(self, N: np.ndarray) -> tuple:
        d32 = self._sauter_mean(N)

        num43 = np.sum(N * self.d_pivots**4)
        den43 = np.sum(N * self.d_pivots**3)
        d43 = num43 / den43 if den43 > 0 else 0.0

        vol_per_bin = N * self.v_pivots
        total_vol = np.sum(vol_per_bin)
        if total_vol <= 0:
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        cum_vol = np.cumsum(vol_per_bin) / total_vol
        d10 = np.interp(0.1, cum_vol, self.d_pivots)
        d50 = np.interp(0.5, cum_vol, self.d_pivots)
        d90 = np.interp(0.9, cum_vol, self.d_pivots)
        span = (d90 - d10) / d50 if d50 > 0 else 0.0

        return (d32, d43, d10, d50, d90, span)


def solve_emulsification(params: SimulationParameters,
                         props: MaterialProperties,
                         phi_d: float | None = None) -> EmulsificationResult:
    """Convenience function to solve Level 1 emulsification.

    If ``phi_d`` is None, resolved per mode (volumetric for stirred-vessel,
    formulation.phi_d for legacy).
    """
    solver = PBESolver(
        n_bins=params.solver.l1_n_bins,
        d_min=params.solver.l1_d_min,
        d_max=params.solver.l1_d_max,
    )
    return solver.solve(params, props, phi_d)
