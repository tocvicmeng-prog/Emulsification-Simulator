"""Smoke gate for the optimization stack (ADR-002 pin verification).

This test is the operational anchor for the botorch / gpytorch / torch pin
declared in ADR-002. It verifies the actual call chain used by
`emulsim.optimization.engine`:

  1. FastNondominatedPartitioning constructs against a tiny synthetic Y.
  2. qLogExpectedHypervolumeImprovement accepts FastNondominatedPartitioning
     as the `partitioning=` argument (the duck-typed path that confuses mypy
     stubs but works at runtime).
  3. optimize_acqf returns a candidate without raising.

If any of the pinned packages bumps the relevant API, this test will fail
loudly before the engine is exercised end-to-end.
"""

from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")
botorch = pytest.importorskip("botorch")


@pytest.mark.smoke
def test_optimization_stack_construction():
    from botorch.acquisition.multi_objective.logei import (
        qLogExpectedHypervolumeImprovement,
    )
    from botorch.models.gp_regression import SingleTaskGP
    from botorch.optim.optimize import optimize_acqf
    from botorch.utils.multi_objective.box_decompositions.non_dominated import (
        FastNondominatedPartitioning,
    )

    ref = torch.tensor([0.0, 0.0], dtype=torch.double)
    Y = torch.tensor(
        [[1.0, 0.5], [0.5, 1.0], [0.7, 0.7]], dtype=torch.double
    )
    X = torch.tensor(
        [[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]], dtype=torch.double
    )

    partitioning = FastNondominatedPartitioning(ref_point=ref, Y=Y)
    assert partitioning is not None

    model = SingleTaskGP(X, Y)
    acqf = qLogExpectedHypervolumeImprovement(
        model=model,
        ref_point=ref.tolist(),
        partitioning=partitioning,
        sampler=None,
    )
    assert acqf is not None

    bounds = torch.stack(
        [torch.zeros(2, dtype=torch.double), torch.ones(2, dtype=torch.double)]
    )
    candidates, _ = optimize_acqf(
        acq_function=acqf,
        bounds=bounds,
        q=1,
        num_restarts=2,
        raw_samples=16,
    )
    assert candidates.shape == (1, 2)
