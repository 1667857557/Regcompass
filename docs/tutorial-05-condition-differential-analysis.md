# Tutorial Level 5: compare reaction support between conditions

Use this tutorial after Stage 6 or a complete one-shot run to compare the same reaction target across conditions within a fixed cell type, target direction, and medium scenario.

The statistical unit is one metacell. Therefore, these tests describe within-dataset condition-associated separation; metacells are not independent biological replicates.

## Load a completed result

```r
library(RegCompassR)

result <- readRDS(
  "RegCompass_result/06_results/regcompass_result.rds"
)

stopifnot(
  identical(result$version, "1.8.3"),
  nrow(result$reaction_ranking) > 0,
  nrow(result$reaction_catalog) > 0,
  nrow(result$reaction_evidence) > 0
)
```

For older result objects without reaction annotations, attach the exact GEM used in the original analysis:

```r
gem <- rc_prepare_gem("human")

result <- rc_attach_reaction_annotations(
  result,
  gem,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem"
)
```

Older Layer 1 objects that lack enough information to reconstruct RNA-only and integrated GPR-aggregated reaction capacities are marked `evidence_resolution = "reaction_capacity_unavailable"`. Gene-level ATAC changes remain visible, but such reactions are not promoted to `RNA+ATAC` without a reaction-level capacity comparison.

## Run all requested condition comparisons

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c(
    "control_24hr",
    "JQ1_24hr",
    "MS177_24hr"
  ),
  cell_types = "stem-cell_like",
  min_units = 5,
  p_adjust_method = "BH",
  p_adjust_scope = "celltype_contrast_medium",
  outdir = "RegCompass_result/07_condition_statistics"
)
```

When three or more conditions are retained, RegCompass returns:

- a Kruskal-Wallis omnibus test for each fixed `cell type × reaction × direction × medium` target;
- pairwise Wilcoxon tests for every requested condition pair;
- adjusted P values, effect sizes, direction of change, reaction annotation, and RNA-versus-multiome evidence provenance.

```r
condition_stats$omnibus
condition_stats$pairwise
```

## Inspect the annotated result columns

```r
condition_stats$omnibus[
  ,
  c(
    "reaction_id",
    "reaction_name",
    "target_direction",
    "tested_formula",
    "genes",
    "gpr_rule",
    "evidence_by_condition",
    "evidence_resolution_by_condition",
    "p_adj"
  )
]

condition_stats$pairwise[
  ,
  c(
    "reaction_id",
    "reaction_name",
    "target_direction",
    "condition_a",
    "condition_b",
    "delta_median_score_b_minus_a",
    "rank_biserial_b_minus_a",
    "evidence_class_a",
    "evidence_class_b",
    "evidence_resolution_a",
    "evidence_resolution_b",
    "p_adj"
  )
]
```

Positive `delta_median_score_b_minus_a` and positive `rank_biserial_b_minus_a` indicate stronger reaction support in `condition_b` than in `condition_a`.

Forward and reverse targets are distinct LP objectives. Do not combine them into a single net-flux estimate.

## Select reactions associated with target genes

Gene selection uses the Boolean GEM GPR annotation rather than reaction-name text matching. Matching is case-insensitive, while reported symbols retain the case stored in the source GEM; this preserves standard mouse symbols such as `Slc22a17`.

```r
rela_metabolic_genes <- c(
  "SLC7A11",
  "GCLC",
  "GCLM",
  "GSS",
  "GSR",
  "G6PD",
  "PGD"
)

gene_reactions <- rc_select_gene_reactions(
  result,
  genes = rela_metabolic_genes,
  match = "any",
  conditions = c(
    "control_24hr",
    "JQ1_24hr",
    "MS177_24hr"
  ),
  cell_types = "stem-cell_like"
)

gene_reactions$reactions[
  ,
  c(
    "reaction_id",
    "reaction_name",
    "model_formula",
    "genes",
    "gpr_rule",
    "matched_genes"
  )
]
```

A matched gene participates in the reaction GPR; it is not necessarily sufficient for reaction activity. Inspect `gpr_rule` for multisubunit complexes and alternative isozymes.

To retain only reactions whose GPR-aggregated capacity is actively changed by ATAC integration in at least one selected group:

```r
multiome_gene_reactions <- rc_select_gene_reactions(
  result,
  genes = rela_metabolic_genes,
  match = "any",
  conditions = c(
    "control_24hr",
    "JQ1_24hr",
    "MS177_24hr"
  ),
  cell_types = "stem-cell_like",
  evidence_class = "RNA+ATAC"
)
```

## Plot one reaction across conditions

```r
p <- rc_plot_condition_reaction(
  result,
  reaction_id = "MAR06231",
  cell_type = "stem-cell_like",
  target_direction = "reverse",
  medium_scenario = "high_glucose",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c(
    "control_24hr",
    "JQ1_24hr",
    "MS177_24hr"
  ),
  min_units = 5,
  annotation_p = "p_adj"
)

print(p)
```

The plot shows one point per metacell and adjusted significance brackets. Multiplicity correction is computed over the full scored reaction family within the selected `p_adjust_scope`, not only over the displayed reaction. The full annotated result is passed into the statistics layer, so reaction names, formulas, GPRs, and evidence classes remain attached to the plot.

## Plot significant reactions for a gene set

```r
gene_plots <- rc_plot_condition_gene_reactions(
  result,
  genes = rela_metabolic_genes,
  cell_type = "stem-cell_like",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c(
    "control_24hr",
    "JQ1_24hr",
    "MS177_24hr"
  ),
  comparisons = list(
    c("control_24hr", "MS177_24hr"),
    c("JQ1_24hr", "MS177_24hr")
  ),
  min_units = 5,
  target_directions = c("forward", "reverse"),
  medium_scenario = "high_glucose",
  evidence_class = "RNA+ATAC",
  p_adj_max = 0.05,
  min_abs_rank_biserial = 0.30,
  max_reactions = 12,
  outdir = paste0(
    "RegCompass_result/07_condition_statistics/",
    "RELA_gene_plots"
  )
)

names(gene_plots$plots)
gene_plots$selected_targets
gene_plots$pairwise_hits
print(gene_plots$plots[[1]])
```

The `conditions` filter is applied both to evidence-class selection and to condition testing. `min_units` is forwarded to the same statistical engine used by `rc_test_condition_reactions()`; it is not a hidden fixed threshold.

## Define a biologically interpretable candidate table

Use adjusted significance together with effect size and evidence source:

```r
hits <- subset(
  condition_stats$pairwise,
  p_adj < 0.05 &
    abs(rank_biserial_b_minus_a) >= 0.30 &
    abs(delta_median_score_b_minus_a) >= 0.10 &
    (
      evidence_class_a == "RNA+ATAC" |
        evidence_class_b == "RNA+ATAC"
    )
)

hits <- hits[order(hits$p_adj, -abs(hits$rank_biserial_b_minus_a)), ]
```

The exported statistics explicitly retain:

```text
inference_level = metacell_within_dataset
descriptive_only = TRUE
biological_replicate_inference = FALSE
```

A significant result indicates differential support for sustaining a reaction direction under the shared model and medium. It is not direct observation of intracellular flux. Population-level treatment inference requires independent biological replicates, and flux validation requires targeted metabolomics or isotope tracing.
