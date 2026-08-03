# Violin geometry extension for condition-associated reaction plots.

.rc_plot_condition_reaction_boxplot_impl <- rc_plot_condition_reaction

#' Plot one reaction target across conditions
#'
#' Extends the canonical condition-reaction plot with violin geometries while
#' preserving the full reaction-wide statistics, significance brackets, formal
#' reaction annotations, and attached plot data produced by
#' [rc_plot_condition_reaction()].
#'
#' @param plot_type Geometry used for the condition distributions:
#'   `"violin_boxplot"` (default), `"violin"`, or `"boxplot"`.
#' @param violin_width Width of the violin geometry.
#' @param violin_alpha Violin fill transparency.
#' @param violin_trim Trim violin densities to the observed score range.
#' @inheritParams rc_plot_condition_reaction
#' @return A `ggplot` object with `condition_statistics`, `plot_data`,
#'   `annotation_data`, and `plot_type` attributes.
#' @export
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
    box_width = 0.32,
    plot_type = c("violin_boxplot", "violin", "boxplot"),
    violin_width = 0.82,
    violin_alpha = 0.45,
    violin_trim = FALSE,
    bracket_step = 0.12,
    title = NULL,
    y_label = "Reaction support score") {
  plot_type <- match.arg(plot_type)
  if (!is.logical(violin_trim) || length(violin_trim) != 1L ||
      is.na(violin_trim)) {
    stop("`violin_trim` must be one non-missing TRUE/FALSE value.",
         call. = FALSE)
  }
  if (!is.numeric(violin_width) || length(violin_width) != 1L ||
      !is.finite(violin_width) || violin_width <= 0 ||
      !is.numeric(violin_alpha) || length(violin_alpha) != 1L ||
      !is.finite(violin_alpha) || violin_alpha <= 0 || violin_alpha > 1) {
    stop("Violin geometry controls contain invalid values.", call. = FALSE)
  }

  plot <- .rc_plot_condition_reaction_boxplot_impl(
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

  if (!identical(plot_type, "boxplot")) {
    violin_layer <- ggplot2::geom_violin(
      width = violin_width,
      trim = violin_trim,
      scale = "width",
      alpha = violin_alpha,
      linewidth = 0.45,
      show.legend = FALSE
    )
    layers <- plot$layers
    if (identical(plot_type, "violin")) {
      keep <- !vapply(layers, function(layer) {
        inherits(layer$geom, "GeomBoxplot")
      }, logical(1))
      layers <- layers[keep]
    }
    plot$layers <- c(list(violin_layer), layers)
  }

  attr(plot, "plot_type") <- plot_type
  plot
}
