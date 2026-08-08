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

The pre-parallel persistent-engine implementation remains as the canonical serial route. `parallel = FALSE` or an effective one-worker configuration uses that serial core directly.

## Solver reuse

The LP matrix and all CORDA2 mathematical controls are unchanged. Each directional target assessment receives the same split model, confidence snapshot, objective coefficients, bounds, target constraint and stable-dependency stopping rule as the preserved serial implementation.

HiGHS remains restricted to one solver thread per worker. Targets in each CORDA2 step are divided into a small number of stage-local chunks (up to four chunks per active worker). One persistent HiGHS engine is reused within each chunk, avoiding one solver-model construction per directional target while retaining dynamic scheduling across workers. Every chunk releases its native solver engine before the stage pool is released.

## Parallel scope

CORDA2 keeps the original Step 1 -> Step 2.1 -> Step 2.2 -> Step 3 state barriers. Directional targets inside one step read the same immutable stage snapshot and are evaluated in parallel. Results are placed back into the original directional-target positions and reduced in the original target order before any dependency matrix, confidence class or retained-reaction state is mutated.

Cell-type-by-medium reconstructions execute serially so each active CORDA2 step can consume the full Layer-2 worker budget. For one step, worker concurrency is:

```text
min(number of directional targets in the step, available Layer-2 workers)
```

The caller-supplied `BPPARAM` is a worker-count/backend template only; Layer 2 does not keep a long-lived CORDA2 pool alive. A fresh stage-local pool is created for every mathematical CORDA2 step. After all targets in the step finish, chunk-local HiGHS engines are released, the stage pool is stopped, full garbage collection runs, and only then does the next step start.

## Progress reporting

At entry to each CORDA2 step, Layer 2 prints the step name, directional-target count, remaining target count, active worker count and chunk count. During the step it prints a target-granular ASCII progress bar such as:

```text
CORDA2 Step 2.1 MC-to-NC dependencies [==========>                 ] 812/2134 (38.1%) remaining=1322 current=RXN::forward
```

The progress denominator is the directional-target count for that CORDA2 step, not the outer Layer-2 workflow count. Worker progress events also record `completed`, `remaining` and `current_target` under the `corda2_stage` scope in the Layer-2 task-progress files. At the barrier, Layer 2 records `remaining=0` and reports that the stage-local worker pool has been released before entering the next step.

## Recorded diagnostics

Each cache summary records the number of LP solves, objective updates, bound updates, solver runtime, reaction counts, confidence-class counts and reconstruction-stage counts. Reconstruction output additionally records, for every CORDA2 step, target count, worker count and chunk count, plus the stage-barrier parallel policy and worker lifecycle. The output contract distinguishes the preserved serial persistent-engine route from stage-parallel execution.
