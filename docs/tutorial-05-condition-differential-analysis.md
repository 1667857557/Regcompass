# Tutorial 5: compare reactions between conditions

Use these functions after Stage 6, a one-shot run, or targeted reaction scoring:

- `rc_test_condition_reactions()`;
- `rc_report_condition_directions()`;
- `rc_plot_condition_reaction()`;
- `rc_plot_condition_gene_reactions()`.

Mathematical definitions are in [Mathematical model](mathematical-model.md).

## Inspect available targets

```r
ranking <- result$reaction_ranking
contrast <- result$condition_contrast

head(ranking)
head(contrast)
```

Lower normalized penalty indicates stronger network-constrained support for the
specified reaction direction. It is not a measured flux.

## Select a fixed target

```r
reaction_id <- "MAR04324"
direction <- "forward"
medium_id <- "physiologic"
cell_type <- "stem-cell_like"

one <- ranking[
  ranking$reaction_id == reaction_id &
    ranking$target_direction == direction &
    ranking$medium_scenario == medium_id &
    ranking$cell_type == cell_type,
  , drop = FALSE
]
```

Condition comparisons must use the same reaction, direction, medium, and broad
cell type.

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
`condition_b`.

Pairwise tests use Wilcoxon statistics. With at least three conditions, the
optional omnibus test uses Kruskal-Wallis statistics. Multiple-testing scope is
controlled by `p_adjust_scope`.

## Direction-aware report

Forward and reverse directions are separate targets. Do not add their scores.

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

Inspect:

```r
direction_report$directional_pairwise
direction_report$directional_omnibus
direction_report$reaction_pairwise
direction_report$reaction_omnibus
direction_report$direction_diagnostics
```

`any_direction_support` reports the best-supported available direction.
`directional_balance` reports forward-versus-reverse support asymmetry; it is
not net flux.

## Plot one target

```r
rc_plot_condition_reaction(
  result,
  reaction_id = reaction_id,
  target_direction = direction,
  medium_scenario = medium_id,
  cell_type = cell_type,
  condition_col = "Group",
  annotation_p = "p_adj"
)
```

For gene-based selection and plotting:

```r
rc_plot_condition_gene_reactions(
  result,
  genes = c("SLC7A11", "GCLC"),
  condition_col = "Group",
  medium_scenario = medium_id
)
```

## Verify structural comparability

```r
result$microcompass$model_cache_summary[, c(
  "medium_scenario",
  "file",
  "file_checksum",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions"
)]
```

Compared units must share the same cached model, reaction order, bounds, target
direction, target-flux fraction, and `vmax`.

Metacells are the statistical units. Reported P values describe within-dataset
condition-associated separation and are not sample- or donor-level inference.

Public API: [functions.md](functions.md).
