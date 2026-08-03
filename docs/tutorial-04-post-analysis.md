# Tutorial 4: post analysis

This tutorial starts from a completed RegCompass workflow. It covers focused reaction expansion, condition-level differential analysis, and visualization of selected reactions. It does not refit Pando or rebuild metacells.

## 1. Targeted reaction expansion and rescoring

`rc_regcompass_step_target_union()` uses selected core reactions or metabolic genes as anchors, identifies directly linked non-core reactions through KEGG, Reactome, or master-Rhea annotations, and reuses the exact cell-type-by-medium union GEMs from the original Layer 2 run. FASTCORE is not rerun.

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "run/07_post_analysis/target_union",
  core_reaction_ids = c("MAR00123", "MAR00456"),
  core_genes = c("SLC7A11", "GCLC", "GCLM"),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  )
)
```

Inspect the mapped anchors, expanded targets, and reused-model scores:

```r
targeted$selected_anchor_reactions
targeted$expanded_reaction_catalog
targeted$expanded_scoring_targets
targeted$summary
targeted$microcompass$penalty
```

The expansion remains cell-type-specific. Reactions from different cell-type union GEMs are never merged into one structural model.

## 2. Condition-level reaction table

The complete workflow returns one row per scored reaction target and metacell:

```r
reaction_long <- result$reaction_comparison_by_metacell

reaction_long <- subset(
  reaction_long,
  cell_type == "T_cell" &
    medium == "normal_human_plasma" &
    direction == "forward" &
    condition %in% c("Control", "Treatment") &
    penalty_available
)

reaction_long$support_score <- -log(
  reaction_long$penalty_per_target_flux + 1e-8
)
```

For a fixed cell type, direction, and medium, larger `support_score` means that less multiome-derived discordance penalty is required to sustain the target reaction flux.

## 3. limma differential analysis between conditions

The following analysis treats metacells as the statistical units. It compares all reaction support scores between two conditions within one fixed cell type, direction, and medium.

```r
library(limma)

score_wide <- reshape(
  reaction_long[, c("row_id", "metacell_id", "support_score")],
  idvar = "row_id",
  timevar = "metacell_id",
  direction = "wide"
)

rownames(score_wide) <- score_wide$row_id
score_matrix <- as.matrix(score_wide[, -1, drop = FALSE])
colnames(score_matrix) <- sub("^support_score\\.", "", colnames(score_matrix))

unit_meta <- unique(
  reaction_long[, c("metacell_id", "condition"), drop = FALSE]
)
unit_meta <- unit_meta[match(colnames(score_matrix), unit_meta$metacell_id), ]
stopifnot(identical(colnames(score_matrix), unit_meta$metacell_id))

condition <- factor(
  unit_meta$condition,
  levels = c("Control", "Treatment")
)
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

min_metacells <- 5L
keep <- apply(score_matrix, 1, function(z) {
  counts <- tapply(is.finite(z), condition, sum)
  all(counts >= min_metacells)
})
score_matrix <- score_matrix[keep, , drop = FALSE]

fit <- lmFit(score_matrix, design)
contrast_matrix <- makeContrasts(
  Treatment_vs_Control = Treatment - Control,
  levels = design
)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

limma_result <- topTable(
  fit,
  coef = "Treatment_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
limma_result$row_id <- rownames(limma_result)

reaction_key <- unique(
  reaction_long[, c(
    "row_id", "reaction_id", "direction", "medium", "cell_type"
  )]
)
limma_result <- merge(
  limma_result,
  reaction_key,
  by = "row_id",
  all.x = TRUE,
  sort = FALSE
)

limma_hits <- subset(
  limma_result,
  adj.P.Val < 0.05 & abs(logFC) >= 0.25
)
```

With the contrast `Treatment - Control`:

- positive `logFC` indicates stronger reaction support in Treatment;
- negative `logFC` indicates stronger reaction support in Control;
- `adj.P.Val` is the BH-adjusted limma P value across tested reaction targets.

This is within-dataset metacell-level inference. It is not donor-level biological-replicate inference. When independent donors or samples are available, build a donor-aware design and avoid treating metacells from the same donor as independent biological replicates.

## 4. Package-native nonparametric condition analysis

For pairwise Wilcoxon tests, effect sizes, and optional Kruskal-Wallis tests across three or more conditions, use the package-native function:

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
  p_adjust_scope = "celltype_contrast_medium",
  outdir = "run/07_post_analysis/condition_statistics"
)

condition_stats$pairwise
condition_stats$omnibus
```

Interpret adjusted P values together with the median shift and rank-biserial effect size:

```r
condition_hits <- subset(
  condition_stats$pairwise,
  p_adj < 0.05 &
    abs(rank_biserial_b_minus_a) >= 0.30 &
    abs(delta_median_score_b_minus_a) >= 0.10
)
```

## 5. Violin plot for a selected metabolic reaction

The selected reaction, target direction, medium, and cell type must be explicit. The default geometry combines a violin distribution, a compact boxplot, and individual metacell points.

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
  comparisons = list(c("Control", "Treatment")),
  plot_type = "violin_boxplot",
  annotation_p = "p_adj",
  show_nonsignificant = TRUE
)

print(p)
ggplot2::ggsave(
  "run/07_post_analysis/MAR06231_T_cell_violin.pdf",
  p,
  width = 5.5,
  height = 5
)
```

Available geometries are:

```r
plot_type = "violin_boxplot"
plot_type = "violin"
plot_type = "boxplot"
```

The plotted metacell values and statistical results remain accessible:

```r
attr(p, "plot_data")
attr(p, "condition_statistics")$pairwise
attr(p, "annotation_data")
```

## 6. Reactions associated with selected metabolic genes

Select reactions through their Boolean GPR annotations, then generate one condition plot per significant reaction direction:

```r
gene_reactions <- rc_select_gene_reactions(
  result,
  genes = c("SLC7A11", "GCLC", "GCLM", "GSS", "GSR"),
  match = "any",
  conditions = c("Control", "Treatment"),
  cell_types = "T_cell"
)

gene_plots <- rc_plot_condition_gene_reactions(
  result,
  genes = c("SLC7A11", "GCLC", "GCLM", "GSS", "GSR"),
  cell_type = "T_cell",
  condition_col = "condition",
  celltype_col = "cell_type",
  conditions = c("Control", "Treatment"),
  comparisons = list(c("Control", "Treatment")),
  target_directions = c("forward", "reverse"),
  medium_scenario = "normal_human_plasma",
  p_adj_max = 0.05,
  min_abs_rank_biserial = 0.30,
  max_reactions = 12,
  outdir = "run/07_post_analysis/gene_reaction_plots"
)

names(gene_plots$plots)
gene_plots$selected_targets
```

A reaction difference represents differential model support for a directional reaction target. It is not a direct measurement of net metabolic flux; metabolomics or isotope-tracing validation remains necessary for flux-level conclusions.
