# Condition-associated reaction statistics

`rc_test_condition_reactions()` compares fixed reaction targets between conditions within the same cell type, target direction, and medium. It operates on the completed Layer 2/Stage 6 result and does not rebuild the structural model.

## Basic use

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

Optional target filters are `reaction_ids`, `target_directions`, and `medium_scenarios`. `include_omnibus = TRUE` adds Kruskal-Wallis tests when at least three conditions are retained. `include_scores = TRUE` keeps the filtered unit-level score data in the returned object.

## Reaction annotations and evidence

Completed results can contain:

```r
result$reaction_catalog
result$reaction_evidence
```

`reaction_catalog` contains reaction names, formulas, GPRs, genes, subsystem/role fields, and available external identifiers. `reaction_evidence` records condition-by-cell-type reaction evidence classes and provenance.

For older completed results, annotations can be attached with:

```r
result <- rc_attach_reaction_annotations(
  result,
  gem,
  condition_col = "condition",
  celltype_col = "cell_type"
)
```

Use `rc_select_gene_reactions()` to select GPR-associated reactions and `rc_plot_condition_gene_reactions()` to test/plot selected gene-associated targets.

## Interpretation boundary

Pairwise tests use Wilcoxon rank-sum tests; optional omnibus tests use Kruskal-Wallis. Metacells are within-dataset statistical units, not donor/sample biological replicates. Forward and reverse reaction targets remain separate.

The exact support-score transformation, RNA-versus-multiome evidence definitions, GPR aggregation, structural-model sharing, and Layer 2 penalty normalization are maintained only in [mathematical-model.md](mathematical-model.md).
