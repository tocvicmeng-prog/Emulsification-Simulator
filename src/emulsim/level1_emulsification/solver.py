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
from .energy import emulsion_density, max_dissipation
from .kernels import (
    breakage_rate_alopaeus,
    coalescence_rate_ct,
)


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
        """Build breakage birth matrix using fixed-pivot redistribution."""
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

    def solve(self, params: SimulationParameters,
              props: MaterialProperties,
              phi_d: float = 0.05) -> EmulsificationResult:
        """Solve the PBE for given process conditions."""
        rpm = params.emulsification.rpm
        t_emul = params.emulsification.t_emulsification
        mixer = params.emulsification.mixer

        # Energy dissipation
        rho_emul = emulsion_density(props.rho_oil, props.rho_aq, phi_d)
        epsilon = max_dissipation(mixer, rpm, rho_emul)

        # Breakage rates
        nu_c = props.mu_oil / props.rho_oil
        g = breakage_rate_alopaeus(
            self.d_pivots, epsilon, props.sigma, props.rho_oil,
            props.mu_d, nu_c=nu_c,
        )
        birth_matrix, death_rate = self._build_breakage_matrix(g)

        # Coalescence rate matrix (vectorised construction)
        di_grid, dj_grid = np.meshgrid(self.d_pivots, self.d_pivots, indexing='ij')
        Q = coalescence_rate_ct(
            di_grid, dj_grid,
            epsilon, props.sigma, props.rho_oil,
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

        return EmulsificationResult(
            d_bins=self.d_pivots.copy(),
            n_d=n_d_output,
            d32=d32, d43=d43, d10=d10, d50=d50, d90=d90, span=span,
            total_volume_fraction=total_vol,
            converged=converged,
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
                         phi_d: float = 0.05) -> EmulsificationResult:
    """Convenience function to solve Level 1 emulsification."""
    solver = PBESolver(
        n_bins=params.solver.l1_n_bins,
        d_min=params.solver.l1_d_min,
        d_max=params.solver.l1_d_max,
    )
    return solver.solve(params, props, phi_d)
