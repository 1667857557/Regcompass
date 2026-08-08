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

The canonical `.rc_corda_build_three_stage_core()` contains both execution paths. `parallel = FALSE` or one effective worker keeps the original serial persistent-engine target order. Parallel execution changes scheduling only; the dependency LP, objective coefficients, target constraints, thresholds and confidence transitions are unchanged.

## Solver reuse

The LP matrix is constant within a CORDA2 step. Each parallel worker chunk creates one persistent HiGHS model and reuses it for the directional targets assigned to that chunk, updating only changed objective coefficients and variable bounds through the existing CORDA2 solver routines. If the persistent API is unavailable, the existing one-shot fallback remains in effect.

HiGHS is restricted to one thread per worker. This prevents solver-level oversubscription when the Layer-2 worker budget is used for directional-target parallelism.

## Step-barrier parallel scope

The four original CORDA2 state transitions remain strict barriers:

```text
Step 1 HC dependencies
  -> deterministic ordered reduce
  -> update HC / MC / NC
  -> release worker pool and HiGHS engines
Step 2.1 MC-to-NC dependencies
  -> deterministic ordered reduce
  -> build MCxNC
  -> release worker pool and HiGHS engines
Step 2.2 MC feasibility
  -> deterministic ordered reduce
  -> promote NC / retain feasible MC
  -> release worker pool and HiGHS engines
Step 3 HC-to-OT dependencies
  -> deterministic ordered reduce
  -> finalize included OT directions
  -> release worker pool and HiGHS engines
```

Every target inside a step sees the same immutable pre-step `split`, confidence classes, options and bounds. Worker results are written back to their original target indices before any dependency matrix, promotion or confidence state is changed. The next step cannot start until the preceding step has been completely reduced.

CORDA2 cell-type-by-medium models are therefore executed one at a time during structural reconstruction so that the current model can use the full Layer-2 worker budget inside each step. For a step with `N` candidate directional targets and `W` requested workers, effective target concurrency is bounded by `min(N, W)`. The caller-provided `BPPARAM` is treated as a backend and worker-count template; a fresh step-local pool is created and stopped for every mathematical step, followed by full garbage collection.

Completed cell-type-by-medium models are still saved atomically. Primary multiome scoring and the RNA-only control continue to reuse the same structural-model file and checksum.

## Progress reporting

Entering a CORDA2 step prints its name, candidate target count, remaining target count, effective worker count and chunk count. While the step is running, target completion markers provide a reaction-granular progress line such as:

```text
CORDA2 Step 2.1 MC-to-NC dependencies [========>                   ] 428/1964 (21.8%) remaining=1536 current=R123::forward
```

Step completion prints `remaining=0` and confirms that the worker pool was released. The same completed/total/remaining information is also written through the existing Layer-2 task-progress infrastructure with `scope = "corda2_stage"`.

## Recorded diagnostics

Each cache summary continues to record the number of LP solves, objective updates, bound updates, solver runtime, reaction counts, confidence-class counts and reconstruction-stage counts. Parallel reconstructions additionally record stage-level target counts, worker counts, chunk counts, ordered-reduce policy and worker lifecycle. These diagnostics audit scheduling without changing the mathematical program.
