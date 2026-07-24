# Safe biological annotation layer for the existing condition-reaction plot.
# Raw microCOMPASS inputs remain supported when no reaction catalog is attached.

rc_plot_condition_reaction <- function(
    x,
    reaction_id,
    cell_type,
    target_direction = NULL,
    medium_scenario = NULL,
    condition_col = NULL,
    celltype_col = NULL,
    conditions = NULL,
    comparisons = NULL,
    min_units = 5L,
    p_adjust_method = "BH",
    p_adjust_scope = c(
      "celltype_contrast_medium", "celltype_contrast", "celltype", "global"
    ),
    annotation_p = c("p_adj", "p_value"),
    significance_threshold = 0.05,
    show_nonsignificant = FALSE,
    show_omnibus = TRUE,
    point_size = 1.8,
    point_alpha = 0.75,
    jitter_width = 0.12,
    box_width = 0.55,
    bracket_step = 0.12,
    title = NULL,
    y_label = "Reaction support score") {
  plot <- .rc_condition_plot_core(
    x = x,
    reaction_id = reaction_id,
    cell_type = cell_type,
    target_direction = target_direction,
    medium_scenario = medium_scenario,
    condition_col = condition_col,
    celltype_col = celltype_col,
    conditions = conditions,
    comparisons = comparisons,
    min_units = min_units,
    p_adjust_method = p_adjust_method,
    p_adjust_scope = p_adjust_scope,
    annotation_p = annotation_p,
    significance_threshold = significance_threshold,
    show_nonsignificant = show_nonsignificant,
    show_omnibus = show_omnibus,
    point_size = point_size,
    point_alpha = point_alpha,
    jitter_width = jitter_width,
    box_width = box_width,
    bracket_step = bracket_step,
    title = title,
    y_label = y_label
  )

  statistics <- attr(plot, "condition_statistics")
  annotation_row <- data.frame()
  if (!is.list(statistics) || !is.data.frame(statistics$pairwise)) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }
  required <- c(
    "reaction_id", "cell_type", "target_direction", "medium_scenario",
    "reaction_name", "tested_formula", "genes", "evidence_comparison"
  )
  if (!all(required %in% colnames(statistics$pairwise))) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }

  keep <- as.character(statistics$pairwise$reaction_id) == reaction_id &
    as.character(statistics$pairwise$cell_type) == cell_type
  if (!is.null(target_direction)) {
    keep <- keep &
      as.character(statistics$pairwise$target_direction) == target_direction
  }
  if (!is.null(medium_scenario)) {
    keep <- keep &
      as.character(statistics$pairwise$medium_scenario) == medium_scenario
  }
  annotation_row <- statistics$pairwise[keep, , drop = FALSE]
  if (!nrow(annotation_row)) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }
  annotation_row <- annotation_row[1L, , drop = FALSE]
  reaction_name <- as.character(annotation_row$reaction_name[[1L]])
  evidence_text <- as.character(annotation_row$evidence_comparison[[1L]])

  if (is.null(title) && length(reaction_name) == 1L &&
      .rc_ra_nonempty(reaction_name)) {
    plot <- plot + ggplot2::labs(
      title = paste0(
        reaction_name, " (", reaction_id, ", ",
        annotation_row$target_direction[[1L]], ") in ", cell_type
      )
    )
  }
  caption <- .rc_ra_plot_caption(annotation_row, evidence_text)
  if (!is.null(caption) && length(caption) == 1L && !is.na(caption)) {
    plot <- plot + ggplot2::labs(caption = caption) + ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  }
  attr(plot, "reaction_annotation") <- annotation_row
  plot
}
