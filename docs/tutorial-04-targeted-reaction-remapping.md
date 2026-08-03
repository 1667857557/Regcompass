# Tutorial 4: targeted reaction remapping

Use targeted remapping after a completed RegCompass run when a focused reaction list is required.

## Define targets

```r
targets <- data.frame(
  reaction_id = c("MAR00123", "MAR00456"),
  target_direction = c("forward", "both"),
  stringsAsFactors = FALSE
)
```

Reaction identifiers must map to the selected GEM. Keep direction explicit when reversible reactions are interpreted separately.

## Run targeted remapping

```r
targeted <- rc_regcompass_targeted_reactions(
  result = result,
  gem = gem,
  target_reactions = targets,
  outdir = "RegCompass_targeted",
  solver = "highs"
)
```

The function reuses compatible cached structural models and medium settings. A different GEM, medium table, cell-type catalogue, direction setting, or cache checksum requires rebuilding the affected model.

## Inspect outputs

```r
targeted$target_reactions
targeted$reaction_ranking
targeted$condition_summary
targeted$condition_contrast
targeted$model_cache_summary
```

Targeted remapping does not refit Pando or rebuild metacells.
