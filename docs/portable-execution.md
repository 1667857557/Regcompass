# Portable execution, bundled GEMs, progress, and timing

RegCompassR 1.8.3 removes two setup assumptions from the canonical workflow:
users no longer need to prepare the default Human-GEM/Mouse-GEM themselves,
and they no longer need to choose a platform-specific parallel backend.

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

`rc_download_species_gem()` remains available for lower-level update and
inspection workflows. `scripts/build-bundled-gems.R` reproduces the package
assets. Model provenance and CC BY 4.0 attribution are recorded in
`inst/extdata/gem/manifest.tsv`.

## Automatic parallel backend

Use `parallel_backend = "auto"` unless a specific backend is required:

```r
rc_parallel_config(workers = 8L, backend = "auto")
```

Resolution rules:

| Operating system | Resolved backend |
|---|---|
| Windows | `BiocParallel::SnowParam(type = "SOCK")` |
| Linux/macOS | `BiocParallel::MulticoreParam` |
| one worker or unavailable BiocParallel | sequential |

Explicit `parallel_backend = "multicore"` is rejected on Windows rather than
silently creating an unsupported backend. The final result records requested
and actual backends and worker counts.

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

The complete workflow reports progress across six stages. Each independently
run stage reports its own start and completion status.

## Timing outputs

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
```

Timing columns include stage, status, start and finish timestamps, elapsed
seconds, formatted elapsed time, OS type, and R version. Failed stages write an
error-status timing row before propagating the original error.
