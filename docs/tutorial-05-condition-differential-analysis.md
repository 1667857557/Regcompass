# Tutorial Level 5: compare reactions between conditions

Use this tutorial after Stage 6, a complete one-shot run, or a target-union
second pass.

Current APIs are `rc_test_condition_reactions()`,
`rc_report_condition_directions()`, `rc_plot_condition_reaction()`, and
`rc_plot_condition_gene_reactions()`; see [functions.md](functions.md).

## Inspect rankings and descriptive contrasts

```r
ranking <- result$reaction_ranking
contrast <- result$condition_contrast

head(ranking)
head(contrast)
```

The primary normalized penalty is:

```text
normalized_penalty = penalty / (omega × vmax)
```

The condition-statistics support score is:

```text
support = -log(normalized_penalty + eps)
```

Lower normalized penalty and higher support indicate stronger multiome support
for the specified reaction direction.

## Select one reaction target

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

## Run direction-specific comparisons

```r
comparison <- rc_test_condition_reactions(
  result,
  reaction_ids = reaction_id,
  target_directions = direction,
  medium_scenarios = medium_id,
  cell_types = cell_type,
  condition_col = "Group",
  comparisons = list(
    c("Control", "JQ1"),
    c("Control", "MS177")
  ),
  include_scores = TRUE
)

comparison$pairwise
comparison$omnibus
```

Each row compares one fixed:

```text
reaction × direction × medium × cell type
```

A positive `delta_median_score_b_minus_a` means stronger support in
`condition_b`. These are direction-specific LP support results, not measured
net fluxes.

## Build the final direction-aware report

For a reversible reaction, forward and reverse are separate counterfactual LP
targets. They may be numerically identical when the shared GEM and evidence
costs cannot distinguish direction. Do not add the two scores.

```r
direction_report <- rc_report_condition_directions(
  result,
  reaction_ids = reaction_id,
  medium_scenarios = medium_id,
  cell_types = cell_type,
  condition_col = "Group",
  conditions = c("Control", "JQ1", "MS177"),
  comparisons = list(
    c("Control", "JQ1"),
    c("Control", "MS177"),
    c("JQ1", "MS177")
  ),
  source_label = "original_layer2_core",
  outdir = "RegCompass_result/07_direction_report"
)
```

Inspect the primary direction-specific tables:

```r
direction_report$directional_pairwise
direction_report$directional_omnibus
```

Inspect the non-additive reaction-level tables:

```r
direction_report$reaction_pairwise
direction_report$reaction_omnibus
```

The reaction-level `report_metric` values are:

```text
any_direction_support = max(forward_support, reverse_support)
directional_balance  = forward_support - reverse_support
```

`any_direction_support` is the best-supported available direction and avoids
double counting identical forward/reverse rows. `directional_balance` describes
support asymmetry and is not net flux.

Inspect direction identifiability:

```r
direction_report$direction_diagnostics[
  ,
  c(
    "reaction_id",
    "condition",
    "direction_pair_status",
    "max_abs_forward_reverse_difference",
    "directionally_indistinguishable",
    "preferred_direction"
  )
]
```

See [Direction-aware final reporting](direction-aware-condition-reporting.md)
for target-union reporting, combined core/non-core testing families, and
interpretation rules.

## Plot one reaction direction

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
cache[, c(
  "medium_scenario",
  "file",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "build_strategy"
)]
```

The same shared structural GEM, bounds, medium, target direction, and target-flux
fraction must be used for every compared unit. Metacell-level P values quantify
within-dataset condition-associated separation and are not biological-replicate
level treatment inference.
