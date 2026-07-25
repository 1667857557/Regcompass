# Tutorial Level 5: compare reactions between conditions

Use this tutorial after Stage 6 or a complete one-shot run. Comparisons must use the same reaction, direction, medium scenario, and cell type.

Within one medium, all conditions share the same structural model. Condition differences therefore reflect the RNA+ATAC penalty matrix. Do not pool results from different media into one structural comparison.

## Inspect rankings and contrasts

```r
ranking <- result$reaction_ranking
contrast <- result$condition_contrast

head(ranking)
head(contrast)
```

The primary normalized score is:

```text
normalized_penalty = penalty / (omega × vmax)
```

Lower normalized penalty indicates stronger multiome support for the target reaction in the fixed model context.

## Select one reaction-direction-medium target

```r
reaction_id <- "MAR04324"
direction <- "forward"
medium_id <- "high_glucose"
cell_type <- "stem-cell_like"

one <- ranking[
  ranking$reaction_id == reaction_id &
    ranking$target_direction == direction &
    ranking$medium_scenario == medium_id &
    ranking$cell_type == cell_type,
  , drop = FALSE
]
```

## Run the comparison helper

```r
comparison <- rc_test_condition_reactions(
  result,
  reaction_ids = reaction_id,
  directions = direction,
  medium_scenarios = medium_id,
  cell_types = cell_type,
  condition_col = "Group"
)

comparison$summary
comparison$contrast
```

The statistical unit is one metacell. Condition-pooled metacells do not replace sample-level biological replication.

## Plot one reaction

```r
rc_plot_condition_reaction(
  result,
  reaction_id = reaction_id,
  target_direction = direction,
  medium_scenario = medium_id,
  cell_type = cell_type,
  condition_col = "Group"
)
```

## Inspect the cached model

```r
cache <- result$microcompass$model_cache_summary
cache[cache$medium_scenario == medium_id, c(
  "medium_scenario",
  "file",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "build_strategy"
)]
```

For a fixed target:

- positive `delta_support` means stronger support in the first condition;
- negative `delta_support` means stronger support in the second condition;
- interpretation is conditional on the selected medium and target-flux requirement.
