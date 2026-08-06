if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required for the reaction rank regression.", call. = FALSE)
}

source(file.path("R", "reaction_rank_plot.R"))
source(file.path("R", "workflow_utils.R"))

metacells <- c("T_1", "T_2", "T_3", "T_4")
conditions <- c("Control", "Control", "Treatment", "Treatment")
support <- list(
  R1 = c(4.0, 4.2, 3.8, 4.1),
  R2 = c(3.0, 3.1, 2.9, 3.2),
  R3 = c(2.0, 2.1, 1.9, 2.2),
  R4 = c(10.0, 10.1, 9.9, 10.2)
)
reaction_long <- do.call(rbind, lapply(names(support), function(reaction_id) {
  data.frame(
    reaction_id = reaction_id,
    direction = "forward",
    medium = "plasma",
    cell_type = "T_cell",
    condition = conditions,
    metacell_id = metacells,
    penalty_available = TRUE,
    penalty_per_target_flux = exp(-support[[reaction_id]]),
    stringsAsFactors = FALSE
  )
}))
rownames(reaction_long) <- NULL

reaction_catalog <- data.frame(
  reaction_id = c("R1", "R2", "R3", "R4"),
  reaction_name = c("Shared reaction", "Shared reaction", "R3", "Structural"),
  forward_formula = c("A -> B", "A -> C", "C -> D", "X -> Y"),
  reverse_formula = c("B -> A", "C -> A", "D -> C", "Y -> X"),
  stringsAsFactors = FALSE
)

reaction_evidence <- data.frame(
  reaction_id = rep(c("R1", "R2", "R3", "R4"), each = 2L),
  condition = rep(c("Control", "Treatment"), 4L),
  cell_type = "T_cell",
  evidence_class = c(
    "RNA-only", "RNA+ATAC",
    "RNA-only", "RNA-only",
    "RNA-only", "RNA-only",
    "structural/no-GPR", "structural/no-GPR"
  ),
  has_active_multiome_contribution = c(
    FALSE, TRUE,
    FALSE, FALSE,
    FALSE, FALSE,
    FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)

result <- list(
  reaction_comparison_by_metacell = reaction_long,
  reaction_catalog = reaction_catalog,
  reaction_evidence = reaction_evidence
)

plot <- plot_top_celltype_reaction_rank(
  result = result,
  cell_type = "T_cell",
  target_direction = "forward"
)
stopifnot(inherits(plot, "ggplot"))
rank_data <- attr(plot, "rank_data")
selection <- attr(plot, "selection")
stopifnot(
  is.data.frame(rank_data),
  nrow(rank_data) == 3L,
  identical(as.character(rank_data$reaction_id), c("R1", "R2", "R3")),
  identical(
    as.character(rank_data$support_class),
    c("Multiome-supported", "RNA-only", "RNA-only")
  ),
  grepl("R1", rank_data$reaction_label_text[rank_data$reaction_id == "R1"],
        fixed = TRUE),
  grepl("R2", rank_data$reaction_label_text[rank_data$reaction_id == "R2"],
        fixed = TRUE),
  identical(
    rank_data$reaction_label_text[rank_data$reaction_id == "R3"],
    "C -> D"
  ),
  !"R4" %in% rank_data$reaction_id,
  identical(selection$top_n, 20L),
  identical(selection$medium_scenario, "plasma"),
  identical(selection$conditions, c("Control", "Treatment")),
  isTRUE(selection$condition_available),
  identical(selection$ranking_statistic, "median_support_score")
)

control_plot <- plot_top_celltype_reaction_rank(
  result = result,
  cell_type = "T_cell",
  target_direction = "forward",
  medium_scenario = "plasma",
  conditions = "Control",
  top_n = 1L
)
control_data <- attr(control_plot, "rank_data")
stopifnot(
  nrow(control_data) == 1L,
  identical(as.character(control_data$reaction_id), "R1"),
  identical(as.character(control_data$support_class), "RNA-only"),
  identical(control_data$n_conditions, 1L)
)

no_condition_result <- result
no_condition_result$reaction_comparison_by_metacell$condition <- NULL
no_condition_result$reaction_evidence$condition <- NULL
no_condition_plot <- plot_top_celltype_reaction_rank(
  result = no_condition_result,
  cell_type = "T_cell",
  target_direction = "forward",
  medium_scenario = "plasma"
)
no_condition_data <- attr(no_condition_plot, "rank_data")
no_condition_selection <- attr(no_condition_plot, "selection")
stopifnot(
  inherits(no_condition_plot, "ggplot"),
  identical(no_condition_selection$conditions, NULL),
  identical(no_condition_selection$condition_available, FALSE),
  all(no_condition_data$n_conditions == 0L)
)
condition_error <- tryCatch(
  {
    plot_top_celltype_reaction_rank(
      result = no_condition_result,
      cell_type = "T_cell",
      target_direction = "forward",
      medium_scenario = "plasma",
      conditions = "Control"
    )
    NULL
  },
  error = identity
)
stopifnot(
  inherits(condition_error, "error"),
  grepl("no usable condition metadata", conditionMessage(condition_error),
        fixed = TRUE)
)

single_condition_result <- result
single_condition_result$reaction_comparison_by_metacell$condition <- "Control"
single_condition_result$reaction_evidence <-
  single_condition_result$reaction_evidence[
    single_condition_result$reaction_evidence$condition == "Control",
    , drop = FALSE
  ]
single_condition_plot <- plot_top_celltype_reaction_rank(
  result = single_condition_result,
  cell_type = "T_cell",
  target_direction = "forward",
  medium_scenario = "plasma"
)
single_condition_selection <- attr(single_condition_plot, "selection")
single_condition_data <- attr(single_condition_plot, "rank_data")
stopifnot(
  identical(single_condition_selection$conditions, "Control"),
  isTRUE(single_condition_selection$condition_available),
  all(single_condition_data$n_conditions == 1L)
)

if (!methods::isClass("Seurat")) {
  methods::setClass("Seurat", slots = c(meta.data = "data.frame"))
}
object_without_condition <- methods::new(
  "Seurat",
  meta.data = data.frame(
    cell_type = c("T_cell", "T_cell"),
    row.names = c("cell_1", "cell_2"),
    stringsAsFactors = FALSE
  )
)
design <- .rc_resolve_condition_design(object_without_condition)
stopifnot(
  identical(design$analysis_mode, "standard_pando"),
  identical(design$condition_levels, "all"),
  identical(design$condition_supplied, FALSE),
  identical(design$fallback_reason, "default_condition_col_absent"),
  design$condition_col %in% colnames(design$object@meta.data),
  all(design$object@meta.data[[design$condition_col]] == "all")
)

explicit_error <- tryCatch(
  {
    .rc_resolve_condition_design(
      object_without_condition,
      condition_col = "missing_group"
    )
    NULL
  },
  error = identity
)
stopifnot(
  inherits(explicit_error, "error"),
  grepl("Explicitly requested", conditionMessage(explicit_error), fixed = TRUE)
)

cat("exported post-analysis reaction rank regression passed\n")
