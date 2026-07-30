# Tutorial 5: compare reactions between conditions

Use this tutorial after Stage 6 or a one-shot run. Mathematical definitions are
in [Tutorial 3](tutorial-03-mathematical-model.md). Public API:
[functions.md](functions.md).

## Primary output

```r
ranking <- result$reaction_ranking
contrast <- result$condition_contrast
```

These tables use the primary `condition_full_oof` penalty. Lower normalized
penalty indicates stronger network-constrained support for a reaction direction;
it is not measured flux.

A valid comparison fixes the reaction, direction, medium and broad cell type:

```r
reaction_id <- "MAR04324"
direction <- "forward"
medium_id <- "normal_human_plasma"
cell_type <- "stem-cell_like"

one <- ranking[
  ranking$reaction_id == reaction_id &
    ranking$target_direction == direction &
    ranking$medium_scenario == medium_id &
    ranking$cell_type == cell_type,
  , drop = FALSE
]
```

Use the exact `medium_scenario` identifier stored in the result. Comparisons must
not mix plasma, culture-challenge, or user-defined backgrounds.

## Pairwise and omnibus tests

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

A positive `delta_median_score_b_minus_a` indicates stronger support in
`condition_b`. Pairwise tests use Wilcoxon statistics; the optional omnibus test
uses Kruskal-Wallis statistics.

## Direction-aware report

Forward and reverse directions are separate LP targets and must not be added.

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
  source_label = "condition_full_oof",
  outdir = "RegCompass_result/07_direction_report"
)
```

## Decomposition outputs

```r
result$common_support_component_summary
result$common_support_component_contrast
result$condition_unique_penalty_increment_summary
result$rna_only_control_summary
```

The common-support tables isolate jointly estimable edges. The condition-unique
increment is the condition-full LP penalty minus the common-support LP penalty;
it is a decomposition, not the primary ranking.

## Structural comparability

```r
result$microcompass$model_cache_summary[, c(
  "medium_scenario",
  "file",
  "file_checksum",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions"
)]
```

Compared metacells must share the same cached model, reaction order, bounds,
target direction, target-flux fraction and `vmax`.

Metacells are the statistical units. Reported P values describe within-dataset
condition-associated separation and are not sample- or donor-level inference.