# CORDA2 performance and parallel invariants

## Exactness boundary

The reconstruction state machine remains the execution semantics of
`resendislab/corda` commit
`c02e06d50606bf93f23d8f2e6d6ade0e996ca70e` for `met_prod = NULL`.
Performance changes must not modify:

- directional variable order;
- confidence mutation order;
- target order;
- penalty values;
- target bounds;
- solver feasibility threshold;
- redundancy stopping rules;
- Stage 1, Stage 2a, Stage 2b or Stage 3 decisions;
- final `max(forward, reverse) == 3` reaction inclusion.

## Dominant runtime

The expensive operation is the sequence of LP optimizations inside one CORDA2
instance. The LP engine is already native C++ through HiGHS. Rewriting the R
control flow in Rust or a separate C++ implementation would not remove the
main optimization cost and would add a second solver integration that is harder
to keep source-identical.

The previous persistent path reused one HiGHS model and simplex basis, but sent
three complete vectors for every solve:

1. the complete objective vector;
2. the complete lower-bound vector;
3. the complete upper-bound vector.

For Human-GEM this repeatedly crossed the R/C++ boundary even though successive
CORDA2 solves usually change only one target bound and a small number of costs.

## Sparse native updates

The persistent solver now records its current native objective and bounds. For
each solve it computes exact changed indices and calls:

- `hi_solver_set_objective()` only for changed objective coefficients;
- `hi_solver_set_variable_bounds()` only for variables whose lower or upper
  bound changed.

All unchanged coefficients remain in the same native HiGHS model. This is
algebraically identical to complete-vector replacement. The same solver,
constraint matrix, variable order, bounds, objective, feasibility tolerance,
serial target order and simplex basis are retained.

Each reconstruction records:

- number of LP solves;
- number of objective coefficients updated;
- number of bound indices updated;
- numeric values that a complete-vector implementation would have transmitted;
- numeric values actually transmitted;
- avoided numeric transfers;
- transmitted fraction of the former full-vector path.

These fields are written into the cell-type-by-medium model-cache summary.

## Native resource release

A CORDA2 instance explicitly calls the HiGHS clear API when reconstruction
finishes. The external pointer finalizer remains as a secondary safeguard. Each
worker saves its completed union GEM atomically, drops the in-memory model, runs
full garbage collection, and then returns to the persistent BiocParallel pool
for another independent task.

If the persistent native API fails, the native model is cleared immediately and
that instance permanently switches to the existing one-shot fallback. A failed
pointer is not retried for every subsequent target.

## Parallel scope

Parallelism is allowed only across independent `cell type × medium` CORDA2
instances. A single instance remains serial because Python CORDA2 mutates
confidence and impossible-target state in target order. Target-level parallelism
would therefore change the algorithm.

The outer dispatcher now parallelizes whenever both conditions hold:

- at least two independent model tasks exist;
- at least two workers are available.

It no longer falls back to serial execution when the number of tasks is smaller
than the number of workers. BiocParallel `tasks` is set to the number of
independent models, so each completed worker receives the next available model
instead of holding a fixed multi-model chunk.

Actual concurrency is:

```text
min(number of cell-type × medium tasks, available workers)
```

One-thread solver settings remain active inside each worker to prevent nested
oversubscription and to preserve deterministic source-parity checks.

## Validation

`tests/corda-performance-check.R` verifies that:

- sparse persistent and one-shot solves return the same status, objective and
  primal solution on a nondegenerate network;
- only changed coefficients and bounds are transmitted;
- transmitted values are less than half of complete-vector transmission in the
  regression sequence;
- native solver state is explicitly cleared;
- `tasks < workers` still selects outer parallel execution;
- one BiocParallel task is created per independent model.

The Python oracle and all existing mathematical parity tests continue to run
after the performance regression. A full Human-GEM wall-time, scaling and peak
memory benchmark is still required before claiming a specific end-to-end speedup
factor.
