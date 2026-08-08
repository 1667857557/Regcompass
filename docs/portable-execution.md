# Portable execution, progress, timing, and worker cleanup

RegCompass automatically resolves the operating-system parallel backend and manages worker pools. Users do not need to construct `SnowParam` or `MulticoreParam` objects for the public workflow.

## One global worker budget

All public parallel stages use one requested worker budget:

```r
workers <- 10L
```

`10` is the default request. It can be changed per call or once per session:

```r
options(RegCompassR.workers = 60L)
```

RegCompass detects the scheduler/machine CPU allocation and always reserves two CPUs for the controller process, operating system, garbage collection, file I/O and other native work. The hard ceiling is therefore:

```text
worker_ceiling = max(1, available_cpus - 2)
resolved_workers = min(requested_workers, worker_ceiling)
```

Every concrete operation then uses:

```text
min(number_of_independent_tasks, resolved_workers)
```

Thus a 64-CPU allocation with `workers = 60L` uses at most 60 workers, while `workers = 64L` is capped at 62. A Stage 1 run with six broad-cell-type Pando jobs launches only six workers even if the resolved budget is larger. CORDA2 and large directional LP target sets can use the complete resolved budget.

Inspect the resolution without starting workers:

```r
rc_parallel_config(workers = 60L)
```

The returned object records `available_cpus`, `reserved_cpus`, `worker_ceiling`, `worker_budget` and the resolved backend.

## Automatic backend resolution

| Operating system | Resolved backend |
|---|---|
| Windows | `BiocParallel::SnowParam(type = "SOCK")` |
| Linux/macOS | `BiocParallel::MulticoreParam` |
| one resolved worker or unavailable BiocParallel | sequential |

No public workflow call needs a `BPPARAM` argument.

## Stage-specific use of the same budget

The budget is global, but each stage can consume a different number of workers because its independent task count differs.

- **Stage 1 Pando:** parallelized across independent broad cell types. Actual workers are `min(n_celltype_jobs, resolved_workers)`. Pando jobs do not start nested worker pools.
- **Stage 4 Layer 1:** reaction-support tasks use at most the same resolved budget.
- **Stage 5 CORDA2:** outer cell-type-by-medium construction is kept serial so an active CORDA2 mathematical step can use the full resolved budget for directional targets.
- **Stage 5 directional LP scoring:** large independent target sets can use the full resolved budget; smaller sets automatically shrink.

## CORDA2 stage lifecycle

CORDA2 preserves the original MATLAB Step 1 → Step 2.1 → Step 2.2 → Step 3 barriers. Within one mathematical step, independent directional targets are parallelized. Results are restored to original target order before any confidence/dependency state is mutated.

For each CORDA2 step RegCompass performs:

```text
print step name and target count
→ resolve effective workers <= global budget <= available CPUs - 2
→ create stage-local worker pool
→ keep HiGHS at one thread per worker
→ process directional-target chunks
→ update completed/remaining progress
→ reduce results in original target order
→ release all chunk solver engines
→ stop the stage pool
→ gc(full = TRUE)
→ enter the next CORDA2 step
```

No CORDA2 worker pool remains alive across mathematical steps.

## One outer worker equals one single-thread task

To prevent nested oversubscription, package-managed workers force internal numerical libraries and solvers to one thread. RegCompass temporarily sets common OpenMP/BLAS thread environment variables to 1 and keeps `mc.cores = 1L` inside a worker. HiGHS uses one solver thread per worker.

## Complete workflow example

```r
result <- rc_run_regcompass(
  object = A,
  gem = gem,
  outdir = "run",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  workers = 10L
)
```

The same `workers` value is passed through Stage 1, Stage 4 and Stage 5. Major workflow stages are separated by full garbage collection so an upstream pool or large temporary allocation is not intentionally retained into the next stage.

## Progress controls

Every public analysis stage accepts:

```r
progress = TRUE
```

Disable progress globally with:

```r
options(RegCompassR.progress = FALSE)
```

CORDA2 additionally records target-level `completed`, `remaining` and `current_target` information in the Layer 2 progress files. The outer Layer 2 12-part progress monitor remains separate from the mathematical CORDA2 step progress.

## Timing and execution provenance

Every stage writes `step_timing.tsv`. A complete run also records the resolved global worker budget and backend in result provenance. Relevant fields include:

```r
result$params$workers
result$params$parallel_backend
result$params$parallel_policy
```

Stage-specific outputs additionally record detected CPU capacity and worker ceilings where relevant, including `available_cpus`, `reserved_cpus` and `worker_ceiling`.
