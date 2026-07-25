# Tutorial Level 5: compare reaction support between conditions

Use this tutorial after Stage 6 or a complete one-shot run to compare the same reaction target across conditions within a fixed cell type, target direction, and medium scenario.

## Structural comparison rule

All conditions evaluated under one medium scenario use the same medium-specific union GEM. That model is constructed once from the merged biological meta-module catalogue and one global FASTCORE completion. Condition differences therefore arise from the RNA+ATAC penalty matrix rather than condition-specific network structures.

Do not combine results across different medium scenarios as though they shared one structural model.

## Inspect descriptive rankings

```r
ranking <- result$reaction_ranking
contrast <- result$condition_contrast

head(ranking)
head(contrast)
```

The primary normalized quantity is the minimum evidence-discordance cost per unit required near-maximal target flux:

```text
normalized_penalty = penalty / (omega × vmax)
```

Lower normalized penalty indicates stronger multiome support for the target reaction in the fixed union-GEM context.

## Select one reaction-direction-medium target

```r
reaction_id <- "MAR04324"
direction <- "forward"
medium_id <- "high_glucose"
cell_type <- "stem-cell_like"

one <- subset(
  ranking,
  reaction_id == !!reaction_id &
    target_direction == !!direction &
    medium_scenario == !!medium_id &
    cell_type == !!cell_type
)
```

In base R, avoid tidy-evaluation syntax:

```r
one <- ranking[
  ranking$reaction_id == reaction_id &
    ranking$target_direction == direction &
    ranking$medium_scenario == medium_id &
    ranking$cell_type == cell_type,
  , drop = FALSE
]
```

## Run the package comparison helper

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

The statistical unit is one metacell. These tests describe within-dataset condition-associated separation; condition-pooled metacells are not independent biological replicates and do not replace sample-level replication.

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

## Verify the shared structural model

```r
cache <- result$microcompass$model_cache_summary
cache[cache$medium_scenario == medium_id, c(
  "medium_scenario",
  "file",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "build_strategy"
)]

union_gem <- readRDS(
  cache$file[match(medium_id, cache$medium_scenario)]
)

stopifnot(
  isTRUE(union_gem$is_union_gem),
  identical(
    union_gem$union_gem_medium_scenario,
    medium_id
  )
)
```

## Interpret contrasts

For a fixed reaction, direction, medium, and cell type:

- positive `delta_support` means stronger support in the first condition;
- negative `delta_support` means stronger support in the second condition;
- the comparison is conditional on the same union GEM and the same `vmax` target requirement.

The Stage 3 object is not involved in condition-specific LP structure. It supplies only the merged biological reaction catalogue and core target set used to build the Stage 5 union GEM.
