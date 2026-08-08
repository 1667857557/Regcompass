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
started. RegCompass reserves two detected logical CPUs, so every task-level
dispatch is bounded by:

```text
min(number of independent tasks, workers, max(1, detected logical CPUs - 2))
```

Examples on a machine/allocation exposing at least 62 logical CPUs:

- 6 standard-Pando broad-cell-type jobs with `workers = 60L` use 6 workers;
- 30 condition-by-cell-type candidate/fit jobs with `workers = 60L` use 30 workers;
- a CORDA2 step with 2,000 directional targets and `workers = 60L` can use 60 workers;
- a directional LP batch with 8 tasks and `workers = 60L` uses 8 workers.

For a condition GRN, pooled-background and per-condition candidate-discovery jobs
are completed first. RegCompass then freezes one exact TF-peak-target dictionary
per broad cell type at a strict barrier before condition-by-cell-type fixed-
dictionary GLMs are dispatched. Different cell types keep separate Pando objects
and peak/motif feature spaces. Individual Pando jobs do not start nested pools.

Stage 2 fragment aggregation uses the same top-level cap. There is no separate
public `metacell_args$fragment_args$workers`; pass `workers` at the stage or
workflow level instead.

This replaces separate upstream and Layer-2 worker budgets. The same worker cap
is passed to every computational stage, while each stage chooses the actual
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
capacity after reserving two logical CPUs. For example, `workers = 60L` on a job
allocation exposing 32 CPUs uses at most 30 workers; on a 64-CPU allocation it
can use up to 60.

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

It also sets `mc.cores = 1L`. HiGHS and Gurobi LP calls are explicitly limited to
one solver thread per outer worker. This prevents a Pando, CORDA2 or LP task pool
from expanding into nested numerical thread pools.

## Stage- and CORDA2-step-scoped worker lifecycle

Parallel workers are not retained across unrelated stages.

For task-parallel stages RegCompass performs:

```text
resolve operating-system backend
→ reserve two logical CPUs
→ apply the user worker cap
→ count independent tasks
→ create min(tasks, protected worker cap) workers
→ execute tasks
→ stop worker pool
→ release references
→ full garbage collection
```

Stage 2 fragment aggregation receives the same protected cap through
`SuperCell::AggregateFragmentFile(nb_cl = ...)`.

Layer 2 follows the same principle. CORDA2 has an additional strict mathematical
barrier between Step 1, Step 2.1, Step 2.2 and Step 3. Each CORDA2 step creates
its own worker pool, completes and reduces all target results in deterministic
order, releases the pool, performs full garbage collection, and only then starts
the next step.

Every directional CORDA2 target starts from a fresh target-local solver engine.
Repeated maximize/dependency solves belonging to that same target reuse its
engine, but a target never inherits a simplex basis from a previous target or
from its chunk assignment. The target-local engine is released when that target
finishes.

No upstream worker pool remains active while Layer 2 is running, and no CORDA2
step pool remains active after that step ends.

## Dynamic scheduling examples

For a stepwise analysis, define one value and reuse it:

```r
workers <- 60L

step1 <- rc_regcompass_step_grn(..., workers = workers)
step2 <- rc_regcompass_step_metacells(..., workers = workers)
step4 <- rc_regcompass_step_layer1(..., workers = workers)
step5 <- rc_regcompass_step_layer2(..., workers = workers)
```

If `workers` is omitted from these functions, the default requested cap is 10.

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
result$params$requested_workers
result$params$detected_cpu_capacity
result$params$reserved_cpus
result$params$parallel_backend
result$params$parallel_worker_policy
```

Stage objects also record their worker cap/backend where relevant. Layer 2 records
CORDA2 stage worker counts, target-local solver-state scope, and worker lifecycle
in its completion and reconstruction diagnostics.
