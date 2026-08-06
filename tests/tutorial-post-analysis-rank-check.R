document_path <- file.path("docs", "tutorial-04-post-analysis.md")
if (!file.exists(document_path)) {
  stop("tutorial-04-post-analysis.md is unavailable.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required for the tutorial rank regression.", call. = FALSE)
}

lines <- readLines(document_path, warn = FALSE)
begin <- which(lines == "# BEGIN top-celltype-reaction-rank")
end <- which(lines == "# END top-celltype-reaction-rank")
if (length(begin) != 1L || length(end) != 1L || begin >= end) {
  stop("The tutorial rank helper markers are malformed.", call. = FALSE)
}

environment <- new.env(parent = globalenv())
eval(parse(text = lines[begin:end]), envir = environment)
plot_rank <- get(
  "plot_top_celltype_reaction_rank",
  envir = environment,
  inherits = FALSE
)

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

plot <- plot_rank(
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
  identical(selection$ranking_statistic, "median_support_score")
)

control_plot <- plot_rank(
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
  identical(as.character(control_data$support_class), "RNA-only")
)

cat("tutorial post-analysis reaction rank regression passed\n")
