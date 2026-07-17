# Supported functions

Only four functions form the supported API. Internal helpers may change and
should not be called directly.

After installation, each supported function has a standard R help page. Open
it in RStudio with `?rc_run_regcompass`, `?rc_run_regcompass_one_shot`,
`?rc_prepare_human2_gem`, or `?rc_make_medium_scenarios`.

## `rc_run_regcompass_one_shot()`

The tutorial entry point. It prepares a Human-GEM model and shared medium when
needed, then calls `rc_run_regcompass()`.

```r
result <- rc_run_regcompass_one_shot(
  object, "RegCompass_result", motifs, genome, fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type"
)
```

## `rc_prepare_human2_gem()`

Downloads and converts one Human-GEM release.

```r
gem <- rc_prepare_human2_gem(version = "2.0.0")
```

## `rc_make_medium_scenarios()`

Creates the shared exchange constraints used by every condition.

```r
medium <- rc_make_medium_scenarios(gem, scenario = "compass_model_bounds")
```

Current named backgrounds include `normal_human_plasma`, `rpmi1640`,
`low_glucose`, and `high_lactate`. Retired names are not compatibility aliases.
Use `custom_medium` only when measured or justified bounds are available.

## `rc_run_regcompass()`

Runs the canonical workflow with an explicit GEM and medium:

```r
result <- rc_run_regcompass(
  object, gem, "RegCompass_result", motifs, genome, fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  medium_scenarios = medium
)
```

The main path is strict-stratum metacells and Pando inference, local FASTCORE,
global calibration and a shared GEM, followed by directional scoring. The
[parameter-selection guide](../README.md#choosing-analysis-parameters) explains
how to choose metacell resolution, Pando thresholds, calibration options, and
the LP solver. Project code should handle downstream reporting and statistics.
