# Public functions

RegCompassR exposes the main workflow plus two setup helpers used by the tutorial.

## `rc_prepare_human2_gem()`

Downloads and converts a Human-GEM release into the RegCompass GEM format.

```r
gem <- rc_prepare_human2_gem(version = "2.0.0")
```

## `rc_make_medium_scenarios()`

Builds shared medium constraints for all conditions before the main workflow.

```r
medium <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "blood_like"
)
```

## `rc_run_regcompass()`

Runs the complete RNA+ATAC metacell workflow:

```text
condition × sample × cell type
→ metacell construction
→ Pando GRN and reaction confidence
→ local FASTCORE-completed meta-modules
→ all-strata barrier
→ sample-balanced global calibration
→ shared union GEM
→ metacell-specific directional scoring
```

```r
result <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  medium_scenarios = medium
)
```

Downstream exports and statistics are intentionally left to project-specific code so the supported API remains focused on setup plus the canonical workflow.
