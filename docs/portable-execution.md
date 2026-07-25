# Portable execution, bundled GEMs, progress, timing, and worker cleanup

RegCompassR 1.8.4 removes three setup assumptions from the canonical workflow:

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

Rebuild a pinned model from the official upstream repository:

```r
human_gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "download",
  force_download = TRUE
)
```

`rc_download_species_gem()` remains available for lower-level update and inspection workflows. `scripts/build-bundled-gems.R` reproduces the package assets. Model provenance and CC BY 4.0 attribution are recorded in `inst/extdata/gem/manifest.tsv`.

## Canonical two-layer worker model

The complete workflow exposes two worker counts only:

```r
upstream_workers <- 6L
layer2_workers <- 30L
```

`upstream_workers` applies to:

- condition-by-cell-type Pando GRNs;
- Layer 1 reaction-support calculation.

Stage 3 meta-module construction does not run FASTCORE and does not own a worker pool for feasibility completion.

`layer2_workers` applies to:

- one medium-specific union-GEM construction task per medium;
- the single global FASTCORE completion performed within each union-GEM build;
- directional LP scoring across metacells.

Set both values to one for fully serial execution:

```r
result <- rc_run_regcompass_one_shot(
  ...,
  upstream_workers = 1L,
  layer2_workers = 1L
)
```

## Automatic backend resolution

The canonical workflow always requests `backend = "auto"` internally.

| Operating system | Resolved backend |
|---|---|
| Windows | `BiocParallel::SnowParam(type = "SOCK")` |
| Linux/macOS | `BiocParallel::MulticoreParam` |
| one worker or unavailable BiocParallel | sequential |

The public complete-workflow interface therefore does not require `parallel_backend`. Low-level `rc_parallel_config()` remains available for diagnostics and package development.

## One outer worker equals one single-thread task

RegCompass uses task-level parallelism only. Every analysis running inside an outer worker is constrained to one internal thread.

The workflow temporarily sets:

```text
OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
RCPP_PARALLEL_NUM_THREADS=1
```

It also sets `mc.cores = 1L` and keeps Pando's internal `parallel = FALSE`. Package-managed child processes inherit the single-thread environment. This prevents an outer pool of LP tasks from expanding into nested multi-threaded solver or BLAS workloads.

HiGHS uses a one-thread control default. GLPK is used as a serial solver backend. Alternative low-level solver interfaces remain available, but the canonical tutorials use HiGHS.

## Stage-scoped worker lifecycle

Parallel workers are not retained across unrelated stages.

For each parallel stage RegCompass performs:

```text
resolve operating-system backend
→ set internal thread count to one
→ create worker pool
→ start worker pool
→ execute independent tasks
→ stop worker pool in guaranteed cleanup
→ remove pool reference
→ run gc(full = TRUE)
```

Cleanup is registered before pool startup. Stage 1, Stage 4, and Stage 5 each receive a fresh pool when parallel work is requested. Stage 3 is a catalogue-construction stage and does not create a local FASTCORE pool. No upstream worker pool remains active while Layer 2 is running.

## Global FASTCORE configuration

The removed interfaces must not be supplied:

```text
layer1_args$local_fastcore
layer1_args$local_fastcore_args
```

Configure the single medium-specific completion through:

```r
layer2_args = list(
  model_params = list(
    completion_time_limit = 600,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE
  )
)
```

Each medium scenario receives its own union GEM. All conditions and metacells evaluated under that medium reuse that one model file.

## Linux numerical-library setup

The package applies single-thread settings during execution. For cluster batch jobs, setting them before launching R remains recommended because it also covers package loading and user-side preprocessing:

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

The complete workflow reports progress across six stages. Each independently run stage reports its own start and completion status.

## Timing and execution provenance

Every stage writes:

```text
<stage-output>/step_timing.tsv
```

A complete run writes:

```text
<outdir>/00_execution_timing.tsv
```

and stores:

```r
result$timing$stages
result$timing$total
result$params$parallel_backend_resolved
result$params$upstream_workers
result$params$layer2_workers
result$params$internal_threads_per_task
result$params$parallel_worker_lifecycle
result$params$parallel_stage_groups
```

Timing columns include stage, status, timestamps, elapsed seconds, formatted elapsed time, OS type, and R version. Failed stages write an error-status timing row before propagating the original error.
