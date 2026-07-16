# Public functions

RegCompassR exports five functions.

## `rc_prepare_human2_gem()`

Downloads a pinned Human-GEM release, converts it to the RegCompass format, and adds GPR and reaction annotations.

```r
gem <- rc_prepare_human2_gem(version = "2.0.0")
```

## `rc_make_medium_scenarios()`

Creates shared medium constraints for the GEM.

```r
medium <- rc_make_medium_scenarios(gem, scenario = "blood_like")
```

## `rc_run_regcompass()`

Runs the complete metacell, Pando, FASTCORE, global calibration and directional scoring workflow.

```r
result <- rc_run_regcompass(...)
```

## `rc_test_microcompass_differential()`

Aggregates scores to biological samples and tests condition effects within each cell type.

```r
differential <- rc_test_microcompass_differential(
  result$microcompass,
  method = "limma_continuous"
)
```

## `rc_export_microcompass()`

Writes score, penalty, feasibility and diagnostic outputs.

```r
rc_export_microcompass(result$microcompass, "export")
```
