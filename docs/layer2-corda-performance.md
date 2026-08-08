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

The LP matrix and all CORDA2 mathematical controls are unchanged. Each directional target assessment receives the same split model, confidence snapshot, objective coefficients, bounds, target constraint and stable-dependency stopping rule as the serial implementation.

The solver is restricted to one thread per worker. This prevents nested HiGHS oversubscription when multiple directional targets are evaluated concurrently.

## Parallel scope

CORDA2 keeps the original Step 1 -> Step 2.1 -> Step 2.2 -> Step 3 state barriers. Directional targets inside one step read the same immutable stage snapshot and are evaluated in parallel. Results are reduced in the original directional-target order before any dependency matrix, confidence class or retained-reaction state is mutated.

Cell-type-by-medium reconstructions execute serially so each active CORDA2 step can consume the full Layer-2 worker budget. For one step, structural concurrency is:

```text
min(number of directional targets in the step, available Layer-2 workers)
```

A fresh stage-local worker pool is created for every CORDA2 step. After all targets in the step finish, the pool is stopped, native worker/solver state is released, full garbage collection runs, and only then does the next step start.

## Progress reporting

At entry to each CORDA2 step, Layer 2 prints the step name, directional-target count, remaining target count and active worker count. The BiocParallel element-level progress bar advances as directional targets complete.

Worker progress events also record `completed`, `remaining` and `current_target` under the `corda2_stage` scope in the Layer-2 task-progress files. At the barrier, Layer 2 records `remaining=0` and reports that the stage-local worker pool has been released before entering the next step.

## Recorded diagnostics

Each cache summary records the number of LP solves, objective updates, bound updates, solver runtime, reaction counts, confidence-class counts and reconstruction-stage counts. Reconstruction output additionally records stage-barrier target parallelism and the stage-local worker lifecycle. These diagnostics support runtime auditing without changing the mathematical program.
