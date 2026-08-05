# CORDA2 performance and parallel invariants

## Mathematical invariants

Performance changes must not alter the original MATLAB CORDA2 contract:

- split only actively reversible reactions;
- close the opposite directional copy for each tested target;
- maximize and constrain the target before dependency minimization;
- preserve the HC, MC, NC and OT directional order;
- preserve Step-1, Step-2.1, Step-2.2 and Step-3 state transitions;
- preserve `MCxNCthresh`, `constraint`, `constrainby`, `om` and `ci` semantics;
- preserve the stable-dependency stopping rule.

## Solver reuse

The LP matrix is constant within one cell-type-by-medium reconstruction. The persistent HiGHS path therefore reuses one native model and updates only changed objective coefficients and variable bounds. If the persistent API fails, that reconstruction switches permanently to the one-shot solver path.

The solver is restricted to one thread per worker. Native solver state is explicitly cleared when reconstruction finishes.

## Parallel scope

Directional targets remain serial within one CORDA2 reconstruction because each step mutates dependency matrices and retained reaction groups in source order. Parallelism is restricted to independent `cell type × medium` reconstructions.

Actual structural concurrency is:

```text
min(number of cell-type × medium tasks, available workers)
```

Completed models are saved atomically. Primary multiome scoring and the RNA-only control reuse the same structural-model file and checksum.

## Recorded diagnostics

Each cache summary records the number of LP solves, objective updates, bound updates, solver runtime, reaction counts, confidence-class counts and reconstruction-stage counts. These diagnostics support runtime auditing without changing the mathematical program.
