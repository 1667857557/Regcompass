# Tutorial 4: post analysis

This tutorial starts from a completed RegCompass result. It covers targeted reaction rescoring, condition statistics, ranking, and plotting. Mathematical definitions of reaction support scores and penalties are in [mathematical-model.md](mathematical-model.md).

## 1. Targeted reaction rescoring

`rc_regcompass_step_target_union()` maps selected reaction or gene anchors to directly linked database reactions and reuses the existing Stage 5 structural models. It does not rerun CORDA2 or FASTCORE.

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "run/07_target_union",
  core_reaction_ids = c("MAR00123", "MAR00456"),
  core_genes = c("SLC7A11", "GCLC", "GCLM"),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  ),
  workers = 10L
)
```

Useful outputs include `selected_anchor_reactions`, `expanded_reaction_catalog`, `expanded_scoring_targets`, `summary`, and `microcompass`.

## 2. Condition statistics

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "condition",
  celltype_col = "cell_type",
  conditions = c("Control", "Treatment"),
  cell_types = "T_cell",
  comparisons = list(c("Control", "Treatment")),
  target_directions = "forward",
  medium_scenarios = "normal_human_plasma",
  min_units = 5L,
  p_adjust_method = "BH",
  p_adjust_scope = "celltype_contrast_medium"
)

condition_stats$pairwise
condition_stats$omnibus
```

Main filters are `reaction_ids`, `target_directions`, `medium_scenarios`, `conditions`, and `cell_types`. `include_omnibus = TRUE` enables Kruskal-Wallis testing when at least three conditions are retained. Metacells are within-dataset statistical units rather than donor-level biological replicates.

## 3. Plot one reaction across conditions

```r
p <- rc_plot_condition_reaction(
  result,
  reaction_id = "MAR06231",
  cell_type = "T_cell",
  target_direction = "forward",
  medium_scenario = "normal_human_plasma",
  condition_col = "condition",
  celltype_col = "cell_type",
  conditions = c("Control", "Treatment"),
  plot_type = "violin_boxplot"
)
print(p)
```

`plot_type` accepts `"violin_boxplot"`, `"violin"`, or `"boxplot"`. Pairwise annotation behavior is controlled by `comparisons`, `annotation_p`, `significance_threshold`, and the P-adjustment arguments.

## 4. Rank reactions within one cell type

```r
p_rank <- plot_top_celltype_reaction_rank(
  result,
  cell_type = "T_cell",
  target_direction = "forward",
  medium_scenario = "normal_human_plasma",
  conditions = c("Control", "Treatment"),
  top_n = 20L
)
print(p_rank)
```

The helper uses the score and evidence classes stored in the completed result; their quantitative definitions are maintained in [mathematical-model.md](mathematical-model.md).

## 5. Select reactions by metabolic genes

```r
selection <- rc_select_gene_reactions(
  result,
  genes = c("SLC7A11", "GCLC", "GCLM"),
  match = "any",
  cell_types = "T_cell"
)

plots <- rc_plot_condition_gene_reactions(
  result,
  genes = c("SLC7A11", "GCLC", "GCLM"),
  cell_type = "T_cell",
  conditions = c("Control", "Treatment"),
  medium_scenario = "normal_human_plasma",
  max_reactions = 12L
)
```

Use [functions.md](functions.md) for the public function index and Rd help for the complete parameter lists.
