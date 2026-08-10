# Portable execution

RegCompass provides species-aware bundled GEM loading, one workflow-level worker cap, automatic platform backend selection, progress reporting, and restartable stage timing.

## Bundled default GEMs

```r
human_gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

mouse_gem <- rc_prepare_gem(
  species = "mouse",
  version = "1.8.0",
  source = "bundled"
)
```

With `source = "auto"`, `rc_prepare_gem()` checks a compatible user cache, then the installed bundled model, then uses the official download path when needed. Low-level download and bundled-manifest helpers are internal maintenance interfaces.

## One worker cap

```r
workers <- 10L
```

The same `workers` value is passed through computational stages. It is an upper bound: each dispatch uses no more workers than the number of independent tasks and the protected detected CPU capacity.

```r
step1 <- rc_regcompass_step_grn(..., workers = workers)
step2 <- rc_regcompass_step_metacells(..., workers = workers)
step4 <- rc_regcompass_step_layer1(..., workers = workers)
step5 <- rc_regcompass_step_layer2(..., workers = workers)
```

`workers = 1L` requests serial execution. Users do not need to create `SnowParam` or `MulticoreParam` objects.

## Backend selection

`rc_parallel_config()` reports the resolved configuration without starting workers:

```r
rc_parallel_config(workers = 10L)
```

Automatic execution uses SOCK/SnowParam workers on Windows and MulticoreParam on Linux/macOS when parallel execution is available. Numerical work inside each outer worker is constrained to a single internal thread to avoid nested oversubscription.

## Progress and timing

Every public workflow stage accepts `progress`. The package-wide default can be changed with:

```r
options(RegCompassR.progress = FALSE)
```

Each stage writes `step_timing.tsv` and persists stage-specific progress/diagnostic files in its output directory. Layer 2 records the resolved structural/scoring execution contract in the returned object.

The CORDA2 mathematical/state-machine contract is documented only in [mathematical-model.md](mathematical-model.md).
