# Tutorial 5: compare reactions between conditions

Use this tutorial after Stage 6 or a one-shot run. Mathematical definitions are
in [Tutorial 3](tutorial-03-mathematical-model.md); the stepwise workflow is in
[Tutorial 2](tutorial-02-stepwise-audit.md).

## Primary reaction output

```r
ranking <- result$reaction_ranking
contrast <- result$condition_contrast
```

These tables use the primary fixed-dictionary condition penalty. The historical
label `condition_full_oof` may still appear in compatibility fields and function
arguments, but the current calculation is not OOF. Lower normalized penalty
indicates stronger network-constrained support for a reaction direction; it is
not measured flux.

A valid comparison fixes reaction, direction, medium, and broad cell type:

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
`condition_b`. Pairwise tests use Wilcoxon rank-sum statistics; the optional
omnibus test uses Kruskal-Wallis statistics.

Metacells are the statistical units. Reported P values describe
within-dataset condition-associated separation and are not sample- or donor-level
biological-replicate inference.

## Direction-aware report

Forward and reverse directions are separate LP targets and must not be added or
subtracted as if they were one scalar flux.

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

`source_label = "condition_full_oof"` is retained for API compatibility and
selects the current primary fixed-dictionary penalty route.

## Inspect condition-GRN edges supporting a reaction

Condition effects are absolute fixed-dictionary GLM coefficients, not deviations
from a pooled coefficient. For a shared edge,

\[
\Delta\beta_{e,c_1,c_2}=
\widehat\beta_{e,c_1}-\widehat\beta_{e,c_2}.
\]

Start from the complete edge table rather than the significant-only table:

```r
edges_all <- result$grn$tf_peak_gene_condition_effect_all
edges_active <- result$grn$tf_peak_gene_condition_effect
```

For a valid edge-level comparison, verify:

- identical `edge_id`, target, TF, region, and ATAC feature mapping;
- estimability in both conditions;
- residual degrees of freedom and rank diagnostics;
- coefficient direction and standard error;
- raw and BH-adjusted P values.

Example:

```r
edge_id <- "TARGET||TF||REGION"
edge_compare <- edges_all[
  edges_all$edge_id == edge_id &
    edges_all$condition %in% c("Control", "JQ1"),
  , drop = FALSE
]

edge_compare[, c(
  "cell_type", "condition", "edge_id", "estimate", "std_err",
  "pval", "padj", "estimable", "direction"
)]
```

The penalty uses each condition's own estimable edges with `padj < 0.05`. A
non-significant edge remains in the complete table and must not be interpreted
as a biological zero. An unavailable edge has `estimate = NA` and is also not a
biological zero.

## Regulatory and RNA-only outputs

```r
result$layer1$gene_regulatory_modifier
result$rna_only_control_summary
```

Targets without a significant estimable condition edge use the neutral RNA-only
fallback. This is a missing-regulatory-evidence policy, not evidence that the
true regulatory effect equals zero.

Historical decomposition fields remain available:

```r
result$common_support_component_summary
result$common_support_component_contrast
result$condition_unique_penalty_increment_summary
```

In the current model, `common` is a compatibility alias of the primary
fixed-dictionary route and the condition-unique increment is a zero compatibility
matrix. These fields no longer represent a jointly estimable shared-slope or OOF
decomposition.

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
target direction, target-flux fraction, and `vmax`. A difference in any of these
quantities invalidates direct interpretation as a condition effect.

## Multiple media

Run condition tests separately for each medium. The medium identifier is part of
the comparison key:

```r
media <- unique(result$reaction_ranking$medium_scenario)

by_medium <- lapply(media, function(medium_id) {
  rc_test_condition_reactions(
    result,
    reaction_ids = reaction_id,
    target_directions = direction,
    medium_scenarios = medium_id,
    cell_types = cell_type,
    condition_col = "Group",
    comparisons = list(c("Control", "JQ1")),
    include_scores = TRUE
  )
})
names(by_medium) <- media
```

A condition effect that changes across media is a model-context interaction. It
should not be reported as a medium-independent change in metabolic flux.

## Reporting checklist

Report at minimum:

- reaction ID and direction;
- medium scenario and structural-model checksum;
- broad cell type and condition labels;
- number of metacells per condition;
- median score or normalized penalty per condition;
- effect direction and adjusted P value;
- whether the reaction's GPR genes had condition-GRN regulatory support or used
  RNA-only fallback;
- that metacells, not biological donors, are the statistical units.


## Valid comparison scope

Condition contrasts are performed for the same reaction, direction, medium,
and cell type. Each row is evaluated only in metacells whose cell type matches
the row's union-GEM scope. Cross-cell-type penalty or `vmax` comparisons are
not treated as condition contrasts.
