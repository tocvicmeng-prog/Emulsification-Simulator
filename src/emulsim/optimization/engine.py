"""Bayesian optimisation engine using BoTorch for multi-objective optimisation.

Uses Expected Hypervolume Improvement (EHVI) acquisition function with
independent GP surrogates per objective.
"""

from __future__ import annotations

import copy
import json
import logging
import time
from pathlib import Path
from typing import Optional

import numpy as np
import torch
from botorch.acquisition.multi_objective import (
    qLogExpectedHypervolumeImprovement,
)
from botorch.models import SingleTaskGP
from botorch.models.model_list_gp_regression import ModelListGP
from botorch.models.transforms.outcome import Standardize
from botorch.optim import optimize_acqf
from botorch.utils.multi_objective.box_decompositions.non_dominated import (
    FastNondominatedPartitioning,
)
from botorch.utils.multi_objective.pareto import is_non_dominated
from botorch.utils.sampling import draw_sobol_samples
from botorch.utils.transforms import normalize, unnormalize
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.fit import fit_gpytorch_mll

from ..datatypes import (
    FullResult,
    OptimizationState,
    SimulationParameters,
)
from ..pipeline.orchestrator import PipelineOrchestrator
from ..properties.database import PropertyDatabase
from .objectives import (
    PARAM_BOUNDS,
    PARAM_NAMES,
    LOG_SCALE_INDICES,
    LOG_SCALE_INDICES_SV,
    compute_objectives,
    check_constraints,
    get_param_bounds,
)

logger = logging.getLogger(__name__)

# Reference point for hypervolume (worst acceptable objectives)
REF_POINT = torch.tensor([5.0, 5.0, 3.0], dtype=torch.double)

tkwargs = {"dtype": torch.double}


def _log_indices_for_mode(mode: str = "rotor_stator_legacy") -> list[int]:
    """Return log-scale parameter indices for the given mode."""
    if mode == "stirred_vessel":
        return LOG_SCALE_INDICES_SV
    return LOG_SCALE_INDICES


def _to_search_space(x: np.ndarray, mode: str = "rotor_stator_legacy") -> np.ndarray:
    """Transform parameters to search space (log-scale where appropriate)."""
    x_ss = x.copy()
    for i in _log_indices_for_mode(mode):
        x_ss[i] = np.log10(x_ss[i])
    return x_ss


def _from_search_space(x_ss: np.ndarray, mode: str = "rotor_stator_legacy") -> np.ndarray:
    """Transform from search space back to physical parameters."""
    x = x_ss.copy()
    for i in _log_indices_for_mode(mode):
        x[i] = 10.0 ** x[i]
    return x


def _get_search_bounds(mode: str = "rotor_stator_legacy",
                       stirrer_type: str = "pitched_blade") -> torch.Tensor:
    """Get parameter bounds in search space (2 x d)."""
    bounds = get_param_bounds(mode, stirrer_type)
    for i in _log_indices_for_mode(mode):
        bounds[i, 0] = np.log10(bounds[i, 0])
        bounds[i, 1] = np.log10(bounds[i, 1])
    return torch.tensor(bounds.T, **tkwargs)  # shape (2, d)


def _x_to_params(x: np.ndarray, template: SimulationParameters) -> SimulationParameters:
    """Convert physical parameter vector to SimulationParameters."""
    params = copy.deepcopy(template)
    params.emulsification.rpm = float(x[0])
    params.formulation.c_span80 = float(x[1])
    # Reconstruct agarose/chitosan from fraction
    frac = float(x[2])
    total = params.formulation.total_polymer
    params.formulation.c_agarose = frac * total
    params.formulation.c_chitosan = (1.0 - frac) * total
    params.formulation.T_oil = float(x[3])
    params.formulation.cooling_rate = float(x[4])
    params.formulation.c_genipin = float(x[5])
    params.formulation.t_crosslink = float(x[6])
    return params


def _build_gp_models(X: torch.Tensor, Y: torch.Tensor,
                     bounds: torch.Tensor) -> ModelListGP:
    """Build independent GP models for each objective."""
    X_norm = normalize(X, bounds)
    models = []
    for i in range(Y.shape[-1]):
        gp = SingleTaskGP(
            X_norm, Y[:, i:i+1],
            outcome_transform=Standardize(m=1),
        )
        mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
        fit_gpytorch_mll(mll)
        models.append(gp)
    return ModelListGP(*models)


class OptimizationEngine:
    """Multi-objective Bayesian optimisation for the emulsification pipeline.

    Optimises 7 process parameters to minimise 3 objectives:
    - Droplet size deviation from 2 µm
    - Pore size deviation from 80 nm
    - Modulus deviation from target
    """

    def __init__(
        self,
        n_initial: int = 15,
        max_iterations: int = 200,
        convergence_tol: float = 0.01,
        output_dir: Optional[Path] = None,
        db: Optional[PropertyDatabase] = None,
        template_params: Optional[SimulationParameters] = None,
    ):
        self.n_initial = n_initial
        self.max_iterations = max_iterations
        self.convergence_tol = convergence_tol
        self.output_dir = Path(output_dir) if output_dir else Path("output/optimization")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.orchestrator = PipelineOrchestrator(db=db, output_dir=self.output_dir / "runs")
        self.template_params = template_params or SimulationParameters()
        # Use mode-aware bounds from the template
        mode = self.template_params.emulsification.mode
        stirrer_type = "pitched_blade"
        if (self.template_params.emulsification.stirrer is not None):
            stirrer_type = self.template_params.emulsification.stirrer.stirrer_type.value
        self.bounds = _get_search_bounds(mode, stirrer_type)
        self._mode = mode

        # Storage
        self.X_observed: list[np.ndarray] = []
        self.Y_observed: list[np.ndarray] = []
        self.results: list[FullResult] = []
        self.hv_history: list[float] = []

    def _evaluate(self, x_physical: np.ndarray) -> tuple[np.ndarray, FullResult]:
        """Run the full pipeline for a parameter vector. Returns (objectives, result)."""
        params = _x_to_params(x_physical, self.template_params)
        params.run_id = f"opt_{len(self.X_observed):04d}"
        result = self.orchestrator.run_single(params)
        objectives = compute_objectives(result)

        # Check feasibility constraints and penalise infeasible points
        feasible, violations = check_constraints(result)
        if not feasible:
            logger.info("Infeasible point: %s", "; ".join(violations))
            objectives = objectives * 10.0 + 5.0

        return objectives, result

    def _generate_initial_points(self) -> np.ndarray:
        """Generate Sobol quasi-random initial parameter vectors in search space."""
        d = len(PARAM_NAMES)
        sobol = draw_sobol_samples(
            bounds=self.bounds,
            n=self.n_initial,
            q=1,
        ).squeeze(1)  # shape (n_initial, d)
        return sobol.numpy()

    def _compute_hypervolume(self, Y: torch.Tensor) -> float:
        """Compute dominated hypervolume of current Pareto front."""
        # Negate because BoTorch assumes maximisation for HV
        neg_Y = -Y
        ref = -REF_POINT
        try:
            partitioning = FastNondominatedPartitioning(ref_point=ref, Y=neg_Y)
            return float(partitioning.compute_hypervolume().item())
        except Exception:
            return 0.0

    def run(self, n_iterations: Optional[int] = None) -> OptimizationState:
        """Run the full optimisation campaign.

        Parameters
        ----------
        n_iterations : int, optional
            Override max_iterations.
        """
        max_iter = n_iterations or self.max_iterations
        total_evals = self.n_initial + max_iter

        logger.info("Starting optimisation: %d initial + %d BO iterations",
                     self.n_initial, max_iter)

        # ── Phase 1: Initial sampling ─────────────────────────────────
        logger.info("Phase 1: Sobol initial sampling (%d points)", self.n_initial)
        X_init_ss = self._generate_initial_points()

        for i in range(self.n_initial):
            x_ss = X_init_ss[i]
            x_phys = _from_search_space(x_ss)
            logger.info("Init %d/%d: %s", i + 1, self.n_initial,
                         {n: f"{v:.2f}" for n, v in zip(PARAM_NAMES, x_phys)})

            t0 = time.perf_counter()
            try:
                objectives, result = self._evaluate(x_phys)
            except Exception as e:
                logger.warning("Evaluation failed: %s — using penalty", e)
                objectives = np.array([10.0, 10.0, 10.0])
                result = None
            dt = time.perf_counter() - t0

            self.X_observed.append(x_ss)
            self.Y_observed.append(objectives)
            if result is not None:
                self.results.append(result)

            logger.info("  → f=[%.3f, %.3f, %.3f] (%.1fs)",
                         objectives[0], objectives[1], objectives[2], dt)

        # ── Phase 2: Bayesian Optimisation ────────────────────────────
        logger.info("Phase 2: Bayesian Optimisation (%d iterations)", max_iter)

        for iteration in range(max_iter):
            X_torch = torch.tensor(np.array(self.X_observed), **tkwargs)
            Y_torch = torch.tensor(np.array(self.Y_observed), **tkwargs)
            # Replace NaN/Inf with penalty values
            nan_mask = ~torch.isfinite(Y_torch)
            if nan_mask.any():
                Y_torch[nan_mask] = 10.0

            # Compute hypervolume
            hv = self._compute_hypervolume(Y_torch)
            self.hv_history.append(hv)

            # Convergence check
            if len(self.hv_history) >= 5:
                recent = self.hv_history[-5:]
                if max(recent) > 0:
                    rel_change = (max(recent) - min(recent)) / max(recent)
                    if rel_change < self.convergence_tol:
                        logger.info("Converged at iteration %d (HV change %.4f < %.4f)",
                                     iteration, rel_change, self.convergence_tol)
                        break

            # Fit GP models
            try:
                model = _build_gp_models(X_torch, Y_torch, self.bounds)
            except Exception as e:
                logger.warning("GP fitting failed: %s — using random point", e)
                x_next_ss = self._generate_initial_points()[0]
                x_next_phys = _from_search_space(x_next_ss)
            else:
                # Optimise acquisition function (qEHVI)
                X_norm = normalize(X_torch, self.bounds)
                neg_Y = -Y_torch  # BoTorch maximises

                partitioning = FastNondominatedPartitioning(
                    ref_point=-REF_POINT, Y=neg_Y,
                )

                acqf = qLogExpectedHypervolumeImprovement(
                    model=model,
                    ref_point=(-REF_POINT).tolist(),
                    partitioning=partitioning,
                    sampler=None,
                )

                # Find next candidate
                d = self.bounds.shape[1]
                standard_bounds = torch.stack([
                    torch.zeros(d, **tkwargs),
                    torch.ones(d, **tkwargs),
                ])
                candidates, _ = optimize_acqf(
                    acq_function=acqf,
                    bounds=standard_bounds,
                    q=1,
                    num_restarts=5,
                    raw_samples=128,
                )

                # Un-normalise
                x_next_norm = candidates.squeeze(0)
                x_next_ss = unnormalize(x_next_norm, self.bounds).detach().numpy()
                x_next_phys = _from_search_space(x_next_ss)

            # Evaluate
            logger.info("BO iter %d/%d: %s", iteration + 1, max_iter,
                         {n: f"{v:.2f}" for n, v in zip(PARAM_NAMES, x_next_phys)})

            t0 = time.perf_counter()
            try:
                objectives, result = self._evaluate(x_next_phys)
            except Exception as e:
                logger.warning("Evaluation failed: %s", e)
                objectives = np.array([10.0, 10.0, 10.0])
                result = None
            dt = time.perf_counter() - t0

            self.X_observed.append(x_next_ss)
            self.Y_observed.append(objectives)
            if result is not None:
                self.results.append(result)

            logger.info("  → f=[%.3f, %.3f, %.3f]  HV=%.4f (%.1fs)",
                         objectives[0], objectives[1], objectives[2], hv, dt)

        # ── Build final state ─────────────────────────────────────────
        X_all = np.array(self.X_observed)
        Y_all = np.array(self.Y_observed)

        # Find Pareto front, filtering out infeasible (penalised) points.
        # Penalised points have all objectives > 5.0; exclude them before
        # computing non-dominated set so the front only contains feasible designs.
        Y_torch = torch.tensor(Y_all, **tkwargs)
        feasible_mask = (Y_torch.max(dim=-1).values <= 5.0)
        if feasible_mask.any():
            Y_feasible = Y_torch[feasible_mask]
            neg_Y_f = -Y_feasible
            pareto_in_feasible = is_non_dominated(neg_Y_f)
            # Map back to full index space
            feasible_indices = torch.where(feasible_mask)[0]
            pareto_full = torch.zeros(len(Y_all), dtype=torch.bool)
            pareto_full[feasible_indices[pareto_in_feasible]] = True
            pareto_idx = pareto_full.numpy()
        else:
            # All points infeasible — fall back to full non-dominated set
            logger.warning("No feasible points found; Pareto front includes penalised points")
            neg_Y = -Y_torch
            pareto_mask = is_non_dominated(neg_Y)
            pareto_idx = pareto_mask.numpy()

        state = OptimizationState(
            X_observed=X_all,
            Y_observed=Y_all,
            pareto_X=X_all[pareto_idx],
            pareto_Y=Y_all[pareto_idx],
            iteration=len(self.hv_history),
            hypervolume=self.hv_history[-1] if self.hv_history else 0.0,
            hypervolume_history=self.hv_history,
            converged=len(self.hv_history) >= 5 and (
                (max(self.hv_history[-5:]) - min(self.hv_history[-5:]))
                / max(max(self.hv_history[-5:]), 1e-10) < self.convergence_tol
            ),
        )

        # Save state
        self._save_state(state)

        # Log Pareto summary
        logger.info("Optimisation complete: %d evaluations, %d Pareto points",
                     len(X_all), pareto_idx.sum())
        for i, idx in enumerate(np.where(pareto_idx)[0]):
            x_phys = _from_search_space(X_all[idx])
            logger.info("  Pareto %d: f=[%.3f, %.3f, %.3f]  RPM=%.0f Span80=%.1f",
                         i + 1, Y_all[idx, 0], Y_all[idx, 1], Y_all[idx, 2],
                         x_phys[0], x_phys[1])

        return state

    def _save_state(self, state: OptimizationState):
        """Save optimisation results to disk."""
        out = {
            "n_evaluations": len(state.X_observed),
            "n_pareto": len(state.pareto_X),
            "hypervolume": state.hypervolume,
            "converged": state.converged,
            "hv_history": state.hypervolume_history,
            "pareto_objectives": state.pareto_Y.tolist(),
            "pareto_params": [
                {n: float(v) for n, v in
                 zip(PARAM_NAMES, _from_search_space(x))}
                for x in state.pareto_X
            ],
        }
        with open(self.output_dir / "optimization_results.json", "w") as f:
            json.dump(out, f, indent=2)

        # Also save raw arrays
        np.savez(
            self.output_dir / "optimization_arrays.npz",
            X_observed=state.X_observed,
            Y_observed=state.Y_observed,
            pareto_X=state.pareto_X,
            pareto_Y=state.pareto_Y,
        )
