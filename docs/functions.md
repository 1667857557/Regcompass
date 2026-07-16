# Public functions

RegCompassR exposes a one-shot tutorial entry point, the canonical workflow,
and two setup helpers for users who want to run the workflow step by step.

## `rc_run_regcompass_one_shot()`

Runs setup and the complete workflow in a single call. If `gem` is omitted, it
prepares the requested Human-GEM release with `rc_prepare_human2_gem()`. If
`medium_scenarios` is omitted, it builds one shared medium table with
`rc_make_medium_scenarios()` before calling `rc_run_regcompass()`.

```r
result <- rc_run_regcompass_one_shot(
  object = object,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  humangem_version = "2.0.0",
  medium_scenario = "blood_like",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type"
)
```

Use this form for tutorials, smoke tests, and first-pass analyses where the
packaged Human-GEM and medium defaults are appropriate.

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

Runs the complete RNA+ATAC metacell workflow after explicit setup:

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

Use one fixed `metacell_args$gamma` across all strict strata. Strata that
produce fewer than `pando_args$min_metacells` are skipped before Pando and
excluded from downstream global calibration and scoring.

Downstream exports and statistics are intentionally left to project-specific code
so the supported API remains focused on setup plus the canonical workflow.
