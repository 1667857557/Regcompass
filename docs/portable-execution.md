# Portable execution, bundled GEMs, progress, timing, and worker cleanup

RegCompass removes three setup assumptions from the canonical workflow:

1. users do not need to prepare the default Human-GEM or Mouse-GEM;
2. users do not choose a platform-specific parallel backend;
3. users do not manually create or retain worker pools across stages.

## Offline default GEMs

The installed package contains validated RegCompass conversions of:

- Human-GEM 2.0.0;
- Mouse-GEM 1.8.0.

```r
human_gem <- rc_prepare_gem("human")
mouse_gem <- rc_prepare_gem("mouse")
rc_bundled_gem_manifest()
```

The default `source = "auto"` order is:

1. compatible user cache;
2. compatible installed bundled model;
3. official download only when no compatible local model exists.

Require offline bundled loading explicitly:

```r
human_gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)
```

## One worker cap for the complete workflow

The complete workflow exposes one parallel parameter only:

```r
workers <- 10L
```

The default is `10L`. Users may raise or lower it, for example:

```r
workers <- 60L
```

`workers` is an upper bound, not a fixed number of processes that must always be
started. Every dispatch uses:

```text
min(number of independent tasks, requested worker cap, system capacity)
```

Examples:

- 6 independent Pando cell-type jobs with `workers = 60L` use 6 workers;
- 30 independent Pando cell-type jobs with `workers = 60L` use 30 workers;
- a CORDA2 step with 2,000 directional targets and `workers = 60L` can use 60 workers;
- a directional LP batch with 8 tasks and `workers = 60L` uses 8 workers.

This replaces separate upstream and Layer-2 worker budgets. The same worker cap
is passed to Stage 1, Stage 4 and Stage 5, while each stage chooses the actual
number of processes required by its current independent-task count.

Set `workers = 1L` for fully serial execution:

```r
result <- rc_run_regcompass_one_shot(
  ...,
  workers = 1L
)
```

## Automatic backend and capacity resolution

RegCompass always requests `backend = "auto"` internally.

| Operating system | Resolved backend |
|---|---|
| Windows | `BiocParallel::SnowParam(type = "SOCK")` |
| Linux/macOS | `BiocParallel::MulticoreParam` |
| one effective worker or unavailable BiocParallel | sequential |

The requested worker cap is additionally limited by scheduler/cgroup/local CPU
capacity. For example, `workers = 60L` on a job allocation exposing only 32 CPUs
uses at most 32 workers.

Users therefore do not need to create `SnowParam` or `MulticoreParam` objects.
`rc_parallel_config()` remains available for diagnostics:

```r
rc_parallel_config(workers = 60L)
```

## One outer worker equals one single-thread task

RegCompass uses task-level parallelism. Every analysis running inside an outer
worker is constrained to one internal numerical/solver thread.

The workflow temporarily sets:

```text
OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
RCPP_PARALLEL_NUM_THREADS=1
HIGHS_THREADS=1
```

It also sets `mc.cores = 1L`. This prevents a pool of Pando, CORDA2 or LP tasks
from expanding into nested BLAS or solver thread pools.

## Stage- and CORDA2-step-scoped worker lifecycle

Parallel workers are not retained across unrelated stages.

For Stage 1 and Stage 4 RegCompass performs:

```text
resolve operating-system backend and worker cap
→ count independent tasks
→ create min(tasks, worker cap) workers
→ execute tasks
→ stop worker pool
→ release references
→ full garbage collection
```

Layer 2 follows the same principle. CORDA2 has an additional strict mathematical
barrier between Step 1, Step 2.1, Step 2.2 and Step 3. Each CORDA2 step creates
its own worker pool, completes and reduces all target results in deterministic
order, releases the pool and worker-local HiGHS engines, performs full garbage
collection, and only then starts the next step.

No Stage 1/4 worker pool remains active while Layer 2 is running, and no CORDA2
step pool remains active after that step ends.

## Dynamic scheduling examples

For a stepwise analysis, define one value and reuse it:

```r
workers <- 60L

step1 <- rc_regcompass_step_grn(
  ...,
  workers = workers
)

step4 <- rc_regcompass_step_layer1(
  ...,
  workers = workers
)

step5 <- rc_regcompass_step_layer2(
  ...,
  workers = workers
)
```

If `workers` is omitted from any of these functions, the default cap is 10.

## Linux numerical-library setup

The package applies single-thread settings during execution. For cluster batch
jobs, setting them before launching R remains recommended because it also covers
package loading and user-side preprocessing:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

## Progress controls

Every public analysis stage accepts:

```r
progress = TRUE
```

Disable progress for batch logs or non-interactive execution:

```r
options(RegCompassR.progress = FALSE)
```

or per call:

```r
result <- rc_run_regcompass_one_shot(..., progress = FALSE)
```

Layer 2 additionally reports CORDA2 step-level directional-target progress with
completed/total and remaining target counts.

## Timing and execution provenance

Every stage writes:

```text
<stage-output>/step_timing.tsv
```

The returned one-shot result records the resolved execution contract, including:

```r
result$params$workers
result$params$parallel_backend
result$params$parallel_worker_policy
```

Stage 1, Stage 4 and Stage 5 also record their worker cap/backend or parallel
contract in the corresponding stage object. Layer 2 records CORDA2 stage worker
counts and worker lifecycle in its completion and reconstruction diagnostics.
