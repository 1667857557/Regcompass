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

A valid condition comparison fixes:

```text
cell type × reaction × direction × medium
```

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

Use the exact `medium_scenario` and `cell_type` identifiers stored in the result.
Comparisons must not mix plasma, culture-challenge, or user-defined backgrounds,
and must not combine reaction rows from different cell-type structural models.

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

The function filters reaction rows by the selected cell type before testing. A
row built from another cell type's union GEM is excluded rather than treated as
a missing observation in the selected cell type.

Metacells are the statistical units. Reported P values describe within-dataset
condition-associated separation and are not sample- or donor-level biological-
replicate inference.

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
from a pooled coefficient. Within one cell type, for a shared edge,

\[
\Delta\beta_{e,t,c_1,c_2}=
\widehat\beta_{e,t,c_1}-\widehat\beta_{e,t,c_2}.
\]

Start from the complete edge table rather than the significant-only table:

```r
edges_all <- result$grn$tf_peak_gene_condition_effect_all
edges_active <- result$grn$tf_peak_gene_condition_effect
```

For a valid edge-level comparison, verify:

- identical cell type, `edge_id`, target, TF, region, and ATAC feature mapping;
- estimability in both conditions;
- residual degrees of freedom and rank diagnostics;
- coefficient direction and standard error;
- raw and BH-adjusted P values.

Example:

```r
edge_id <- "TARGET||TF||REGION"
edge_compare <- edges_all[
  edges_all$cell_type == cell_type &
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
```

fixed-dictionary route and the condition-unique increment is a zero compatibility
matrix. These fields no longer represent a jointly estimable shared-slope or OOF
decomposition. Their summaries retain the same row-level cell-type filter as the
primary route.

## Structural comparability

For `model_mode = "meta_module_gem"`, inspect the cell-type-specific cache:

```r
cache <- result$microcompass$model_cache_summary

cache[, c(
  "cell_type",
  "medium_scenario",
  "file",
  "file_checksum",
  "n_celltype_biological_reactions",
  "n_celltype_fastcore_support_reactions",
  "build_strategy",
  "completion_stage"
)]
```

A valid condition contrast uses one cache row identified by the same
`cell_type × medium_scenario` key. Conditions of that cell type share:

- the same cell-type biological reaction union;
- the same independently completed cell-type FASTCORE model;
- the same model file and checksum;
- the same reaction order and bounds;
- the same target direction and target-flux fraction;
- the same directional `vmax`.

The required structural provenance is:

```text
structural_scope = cell_type_x_medium
shared_across_conditions = TRUE
shared_across_cell_types = FALSE
build_strategy = celltype_medium_union_gem
completion_stage = celltype_specific_fastcore_after_condition_module_union
```

A difference in model file, checksum, reaction order, bounds, direction, medium,
or `vmax` invalidates direct interpretation as a condition effect.

Different cell types deliberately use different structural models. A comparison
between their penalties or `vmax` values is a cross-cell-type model comparison,
not a condition contrast.

## Verify row-to-model assignment

Reaction rows contain cell-type scope in their row IDs and diagnostics:

```r
row_meta <- rc_parse_microcompass_row_id(
  rownames(result$microcompass$penalty)
)

head(row_meta[, c(
  "cell_type", "reaction_id", "target_direction", "medium_scenario"
)])
```

For a selected cell type, only rows with the same `row_meta$cell_type` are used.
The Layer 2 comparison table records the same relationship:

```r
result$reaction_comparison_by_metacell[, c(
  "row_id", "cell_type", "condition", "reaction_id", "direction", "medium"
)]
```

## Multiple media

Run condition tests separately for each medium. The medium identifier is part of
the comparison key, while cell type remains fixed:

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

Each medium uses a separate union GEM for the same cell type and its own
independent FASTCORE completion. A condition effect that changes across media is
a model-context interaction. It should not be reported as a medium-independent
change in metabolic flux.

## Optional targeted-reaction results

Targeted remapping retains cell-type scope in both the mapped relation table and
the scoring rows:

```r
targeted$expanded_scoring_targets[, c(
  "cell_type", "reaction_id", "anchor_core_reaction_ids"
)]

targeted$microcompass$model_cache_summary[, c(
  "cell_type", "medium_scenario", "file", "file_checksum"
)]
```

Targeted condition comparisons follow the same
`cell type × reaction × direction × medium` rule. Targeted models are exact
reuses of the corresponding Stage 5 cell-type caches; FASTCORE is not rerun.

## Reporting checklist

Report at minimum:

- broad cell type, reaction ID, and direction;
- medium scenario and cell-type structural-model checksum;
- condition labels and number of metacells per condition;
- median score or normalized penalty per condition;
- effect direction and adjusted P value;
- `structural_scope = cell_type_x_medium`;
- that FASTCORE was completed independently for the selected cell type and
  medium;
- whether the reaction's GPR genes had condition-GRN regulatory support or used
  RNA-only fallback;
- that metacells, not biological donors, are the statistical units;
- that no cross-cell-type penalty or `vmax` comparison was interpreted as a
  condition effect.
