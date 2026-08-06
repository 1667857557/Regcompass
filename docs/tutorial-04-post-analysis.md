# Tutorial 4: post analysis

This tutorial starts from a completed RegCompass workflow. It covers focused reaction expansion, condition-level differential analysis, same-cell-type reaction ranking, and visualization of selected reactions. It does not refit Pando or rebuild metacells.

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

## 5. Top reactions within one cell type

The following helper ranks reactions within one fixed cell type, target direction, and medium. The ranking statistic is the median metacell support score

\[
-\log\left(\frac{\mathrm{penalty}}{\omega\,v_{\max}}+10^{-8}\right).
\]

Only reactions whose formal evidence class is `RNA-only` or `RNA+ATAC` are eligible. Structural reactions without a GPR and GPR reactions without observed RNA support are excluded rather than being mislabeled as RNA-only.

A reaction is labeled **Multiome-supported** when at least one selected metacell has an active ATAC contribution that changes the GPR-aggregated reaction capacity. The result object stores this criterion as `has_active_multiome_contribution`; the helper applies `any()` across the selected conditions for the requested cell type. Otherwise, an eligible reaction is labeled **RNA-only**.

Reaction labels use the formal independent reaction name from `result$reaction_catalog`. When that name is missing or only repeats the model reaction ID, the direction-specific chemical formula is used. Duplicate names receive the reaction ID as a suffix so every bar remains uniquely identifiable.

```r
# BEGIN top-celltype-reaction-rank
plot_top_celltype_reaction_rank <- function(
    result,
    cell_type,
    target_direction = c("forward", "reverse"),
    medium_scenario = NULL,
    conditions = NULL,
    top_n = 20L,
    label_width = 58L) {
  target_direction <- match.arg(target_direction)
  if (!is.character(cell_type) || length(cell_type) != 1L ||
      is.na(cell_type) || !nzchar(trimws(cell_type))) {
    stop("`cell_type` must be one non-empty cell-type label.", call. = FALSE)
  }
  cell_type <- trimws(cell_type)
  if (length(top_n) != 1L || is.na(top_n) || !is.finite(top_n) || top_n < 1) {
    stop("`top_n` must be one positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)
  if (length(label_width) != 1L || is.na(label_width) ||
      !is.finite(label_width) || label_width < 20) {
    stop("`label_width` must be one integer of at least 20.", call. = FALSE)
  }
  label_width <- as.integer(label_width)

  reaction_long <- result$reaction_comparison_by_metacell
  reaction_catalog <- result$reaction_catalog
  reaction_evidence <- result$reaction_evidence
  required_long <- c(
    "reaction_id", "direction", "medium", "cell_type", "condition",
    "metacell_id", "penalty_available", "penalty_per_target_flux"
  )
  required_catalog <- c(
    "reaction_id", "reaction_name", "forward_formula", "reverse_formula"
  )
  required_evidence <- c(
    "reaction_id", "condition", "cell_type", "evidence_class",
    "has_active_multiome_contribution"
  )
  if (!is.data.frame(reaction_long) ||
      !all(required_long %in% colnames(reaction_long))) {
    stop("The result lacks the annotated metacell reaction table.",
         call. = FALSE)
  }
  if (!is.numeric(reaction_long$penalty_per_target_flux)) {
    stop("`penalty_per_target_flux` must be numeric.", call. = FALSE)
  }
  if (!is.data.frame(reaction_catalog) ||
      !all(required_catalog %in% colnames(reaction_catalog))) {
    stop("The result lacks the formal reaction-name/formula catalog.",
         call. = FALSE)
  }
  if (!is.data.frame(reaction_evidence) ||
      !all(required_evidence %in% colnames(reaction_evidence))) {
    stop("The result lacks reaction-level RNA/ATAC evidence provenance.",
         call. = FALSE)
  }

  selected <- reaction_long[
    trimws(as.character(reaction_long$cell_type)) == cell_type &
      as.character(reaction_long$direction) == target_direction &
      reaction_long$penalty_available %in% TRUE &
      is.finite(reaction_long$penalty_per_target_flux),
    , drop = FALSE
  ]
  if (!is.null(conditions)) {
    conditions <- unique(trimws(as.character(conditions)))
    conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
    if (!length(conditions)) {
      stop("`conditions` did not contain a usable condition label.",
           call. = FALSE)
    }
    selected <- selected[
      as.character(selected$condition) %in% conditions,
      , drop = FALSE
    ]
  }
  if (!nrow(selected)) {
    stop("No scored metacells match the requested cell type and direction.",
         call. = FALSE)
  }

  available_media <- unique(as.character(selected$medium))
  available_media <- available_media[
    !is.na(available_media) & nzchar(available_media)
  ]
  if (is.null(medium_scenario)) {
    if (length(available_media) != 1L) {
      stop(
        "Specify `medium_scenario`; available media are: ",
        paste(available_media, collapse = ", "),
        call. = FALSE
      )
    }
    medium_scenario <- available_media[[1L]]
  }
  if (!is.character(medium_scenario) || length(medium_scenario) != 1L ||
      is.na(medium_scenario) || !nzchar(trimws(medium_scenario))) {
    stop("`medium_scenario` must be one non-empty medium label.",
         call. = FALSE)
  }
  medium_scenario <- trimws(medium_scenario)
  selected <- selected[
    as.character(selected$medium) == medium_scenario,
    , drop = FALSE
  ]
  if (!nrow(selected)) {
    stop("No scored metacells match the requested medium.", call. = FALSE)
  }

  selected$support_score <- -log(
    pmax(as.numeric(selected$penalty_per_target_flux), 0) + 1e-8
  )
  selected <- selected[is.finite(selected$support_score), , drop = FALSE]
  if (!nrow(selected)) {
    stop("No finite reaction support scores remain.", call. = FALSE)
  }

  score_rows <- lapply(
    split(selected, as.character(selected$reaction_id), drop = TRUE),
    function(one) {
      finite <- is.finite(one$support_score)
      data.frame(
        reaction_id = as.character(one$reaction_id[[1L]]),
        median_support_score = stats::median(one$support_score[finite]),
        mean_support_score = mean(one$support_score[finite]),
        median_penalty_per_target_flux = stats::median(
          one$penalty_per_target_flux[finite]
        ),
        n_metacells = length(unique(as.character(one$metacell_id[finite]))),
        n_conditions = length(unique(as.character(one$condition[finite]))),
        stringsAsFactors = FALSE
      )
    }
  )
  scores <- do.call(rbind, score_rows)
  rownames(scores) <- NULL

  selected_conditions <- unique(as.character(selected$condition))
  evidence_selected <- reaction_evidence[
    trimws(as.character(reaction_evidence$cell_type)) == cell_type &
      as.character(reaction_evidence$condition) %in% selected_conditions,
    , drop = FALSE
  ]
  evidence_rows <- lapply(
    split(
      evidence_selected,
      as.character(evidence_selected$reaction_id),
      drop = TRUE
    ),
    function(one) {
      evidence_class <- as.character(one$evidence_class)
      eligible <- any(evidence_class %in% c("RNA-only", "RNA+ATAC"))
      active_multiome <- any(
        one$has_active_multiome_contribution %in% TRUE
      )
      data.frame(
        reaction_id = as.character(one$reaction_id[[1L]]),
        evidence_eligible = eligible,
        support_class = if (active_multiome) {
          "Multiome-supported"
        } else {
          "RNA-only"
        },
        stringsAsFactors = FALSE
      )
    }
  )
  evidence_summary <- if (length(evidence_rows)) {
    do.call(rbind, evidence_rows)
  } else {
    data.frame(
      reaction_id = character(),
      evidence_eligible = logical(),
      support_class = character(),
      stringsAsFactors = FALSE
    )
  }
  rownames(evidence_summary) <- NULL

  ranked <- merge(
    scores,
    evidence_summary,
    by = "reaction_id",
    all.x = TRUE,
    sort = FALSE
  )
  ranked <- ranked[ranked$evidence_eligible %in% TRUE, , drop = FALSE]
  if (!nrow(ranked)) {
    stop(
      "No reactions have RNA-only or active multiome evidence in the selection.",
      call. = FALSE
    )
  }

  catalog_index <- match(
    ranked$reaction_id,
    as.character(reaction_catalog$reaction_id)
  )
  reaction_name <- trimws(as.character(
    reaction_catalog$reaction_name[catalog_index]
  ))
  formula_column <- paste0(target_direction, "_formula")
  reaction_formula <- trimws(as.character(
    reaction_catalog[[formula_column]][catalog_index]
  ))
  name_available <- !is.na(reaction_name) & nzchar(reaction_name) &
    reaction_name != ranked$reaction_id
  formula_available <- !is.na(reaction_formula) & nzchar(reaction_formula)
  reaction_label <- ifelse(
    name_available,
    reaction_name,
    ifelse(formula_available, reaction_formula, ranked$reaction_id)
  )
  duplicate_label <- duplicated(reaction_label) |
    duplicated(reaction_label, fromLast = TRUE)
  reaction_label[duplicate_label] <- paste0(
    reaction_label[duplicate_label], " [",
    ranked$reaction_id[duplicate_label], "]"
  )
  ranked$reaction_name <- reaction_name
  ranked$reaction_formula <- reaction_formula
  ranked$reaction_label_text <- reaction_label

  ranked <- ranked[
    order(
      -ranked$median_support_score,
      -ranked$mean_support_score,
      ranked$reaction_id
    ),
    , drop = FALSE
  ]
  ranked$rank <- seq_len(nrow(ranked))
  ranked <- utils::head(ranked, top_n)
  wrap_one <- function(x) {
    paste(strwrap(x, width = label_width), collapse = "\n")
  }
  ranked$reaction_label <- vapply(
    ranked$reaction_label_text,
    wrap_one,
    character(1)
  )
  ranked$reaction_label <- factor(
    ranked$reaction_label,
    levels = rev(ranked$reaction_label)
  )
  ranked$support_class <- factor(
    ranked$support_class,
    levels = c("RNA-only", "Multiome-supported")
  )

  condition_text <- paste(selected_conditions, collapse = ", ")
  p <- ggplot2::ggplot(
    ranked,
    ggplot2::aes(
      x = reaction_label,
      y = median_support_score,
      fill = support_class
    )
  ) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Median reaction support score",
      fill = "Evidence support",
      title = paste0(
        "Top ", nrow(ranked), " reactions in ", cell_type
      ),
      subtitle = paste0(
        target_direction, " | ", medium_scenario, " | ", condition_text,
        "\nMultiome-supported if any selected metacell changes GPR reaction capacity"
      )
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      plot.title.position = "plot",
      axis.text.y = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5)
    )

  attr(p, "rank_data") <- transform(
    ranked,
    reaction_label = as.character(reaction_label),
    support_class = as.character(support_class)
  )
  attr(p, "selection") <- list(
    cell_type = cell_type,
    target_direction = target_direction,
    medium_scenario = medium_scenario,
    conditions = selected_conditions,
    top_n = top_n,
    ranking_statistic = "median_support_score",
    multiome_rule =
      "any(has_active_multiome_contribution) across selected conditions"
  )
  p
}
# END top-celltype-reaction-rank
```

Generate the default Top 20 plot for one cell type:

```r
library(ggplot2)

rank_plot <- plot_top_celltype_reaction_rank(
  result = result,
  cell_type = "T_cell",
  target_direction = "forward",
  medium_scenario = "normal_human_plasma",
  conditions = c("Control", "Treatment")
)

print(rank_plot)
attr(rank_plot, "rank_data")
attr(rank_plot, "selection")

ggplot2::ggsave(
  "run/07_post_analysis/T_cell_top20_reaction_rank.pdf",
  rank_plot,
  width = 8.5,
  height = 7
)
```

Set `top_n` to display a different number of reactions:

```r
rank_plot_top30 <- plot_top_celltype_reaction_rank(
  result = result,
  cell_type = "T_cell",
  target_direction = "forward",
  medium_scenario = "normal_human_plasma",
  conditions = c("Control", "Treatment"),
  top_n = 30L
)
```

When `conditions = NULL`, all conditions available for the selected cell type, direction, and medium are pooled for ranking and evidence classification. If only one medium is available after filtering, `medium_scenario` may be omitted; otherwise it must be specified explicitly.

## 6. Violin plot for a selected metabolic reaction

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

## 7. Reactions associated with selected metabolic genes

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
