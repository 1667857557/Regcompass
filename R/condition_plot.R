# Plot condition-associated reaction support under the shared union-GEM.

.rc_plot_condition_single_string <- function(value, name, allow_null = FALSE) {
  if (is.null(value) && isTRUE(allow_null)) return(NULL)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop("`", name, "` must be one non-empty string.", call. = FALSE)
  }
  value
}

.rc_plot_condition_target <- function(
    row_meta, reaction_id, target_direction, medium_scenario) {
  selected <- row_meta$reaction_id == reaction_id
  if (!any(selected)) {
    stop("Reaction target not found: ", reaction_id, call. = FALSE)
  }

  available_directions <- unique(as.character(
    row_meta$target_direction[selected]
  ))
  if (is.null(target_direction)) {
    if (length(available_directions) != 1L) {
      stop(
        "`target_direction` is required because reaction ", reaction_id,
        " has multiple scored directions: ",
        paste(available_directions, collapse = ", "),
        call. = FALSE
      )
    }
    target_direction <- available_directions[[1L]]
  }
  selected <- selected & row_meta$target_direction == target_direction
  if (!any(selected)) {
    stop(
      "Reaction ", reaction_id, " was not scored in direction ",
      target_direction, ". Available directions: ",
      paste(available_directions, collapse = ", "),
      call. = FALSE
    )
  }

  available_media <- unique(as.character(row_meta$medium_scenario[selected]))
  if (is.null(medium_scenario)) {
    if (length(available_media) != 1L) {
      stop(
        "`medium_scenario` is required because the target has multiple media: ",
        paste(available_media, collapse = ", "),
        call. = FALSE
      )
    }
    medium_scenario <- available_media[[1L]]
  }
  selected <- selected & row_meta$medium_scenario == medium_scenario
  if (sum(selected) != 1L) {
    stop(
      "Expected exactly one reaction-direction-medium target after filtering.",
      call. = FALSE
    )
  }

  list(
    row_id = row_meta$row_id[selected][[1L]],
    reaction_id = reaction_id,
    target_direction = target_direction,
    medium_scenario = medium_scenario
  )
}

.rc_plot_condition_star <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  "ns"
}

.rc_plot_condition_annotations <- function(
    pairwise, condition_levels, p_column,
    significance_threshold, show_nonsignificant,
    score_min, score_max, bracket_step) {
  if (!nrow(pairwise)) return(data.frame())
  pairwise$annotation_p <- pairwise[[p_column]]
  pairwise$label <- vapply(
    pairwise$annotation_p, .rc_plot_condition_star, character(1)
  )
  pairwise$xmin <- match(pairwise$condition_a, condition_levels)
  pairwise$xmax <- match(pairwise$condition_b, condition_levels)
  pairwise <- pairwise[
    is.finite(pairwise$annotation_p) &
      is.finite(pairwise$xmin) & is.finite(pairwise$xmax),
    , drop = FALSE
  ]
  if (!isTRUE(show_nonsignificant)) {
    pairwise <- pairwise[
      pairwise$annotation_p < significance_threshold,
      , drop = FALSE
    ]
  }
  if (!nrow(pairwise)) return(data.frame())

  swap <- pairwise$xmin > pairwise$xmax
  if (any(swap)) {
    old <- pairwise$xmin[swap]
    pairwise$xmin[swap] <- pairwise$xmax[swap]
    pairwise$xmax[swap] <- old
  }
  pairwise$span <- pairwise$xmax - pairwise$xmin
  pairwise <- pairwise[order(
    pairwise$span, pairwise$xmin, pairwise$xmax, pairwise$annotation_p
  ), , drop = FALSE]

  score_range <- score_max - score_min
  if (!is.finite(score_range) || score_range <= 0) {
    score_range <- max(1, abs(score_max), abs(score_min))
  }
  step <- score_range * bracket_step
  if (!is.finite(step) || step <= 0) step <- 0.1
  pairwise$y <- score_max + step * seq_len(nrow(pairwise))
  pairwise$tip <- step * 0.20
  pairwise$text_y <- pairwise$y + step * 0.12
  pairwise
}

#' Plot one reaction target across conditions
#'
#' Draws a condition-level boxplot for one fixed reaction, direction, medium,
#' and cell type. Every finite metacell score is displayed as a jittered point.
#' Pairwise significance brackets are derived from
#' [rc_test_condition_reactions()] and can use raw or multiplicity-adjusted P
#' values. The statistics layer is run over the full scored reaction set before
#' the selected target is extracted, so adjusted P values retain the requested
#' reaction-wide multiplicity correction. When three or more conditions are
#' present, the Kruskal-Wallis omnibus result is shown in the subtitle by default.
#'
#' The plotted score is
#' `-log(penalty / (omega * vmax) + eps)`, with larger values indicating stronger
#' multiome support for the selected reaction direction under the shared
#' structural GEM.
#'
#' @param x A Layer 2 microCOMPASS result or a complete RegCompass result.
#' @param reaction_id One reaction identifier.
#' @param cell_type One cell-type value.
#' @param target_direction Optional `"forward"` or `"reverse"`. Required when
#'   both directions were scored.
#' @param medium_scenario Optional medium identifier. Required when the selected
#'   target was scored under multiple media.
#' @param condition_col,celltype_col Metadata columns used by the Layer 2 units.
#' @param conditions Optional ordered condition values. The default uses all
#'   conditions represented in the selected cell type.
#' @param comparisons Optional pairwise contrasts passed to
#'   [rc_test_condition_reactions()]. The default uses every condition pair.
#' @param min_units Minimum finite metacells per condition for testing.
#' @param p_adjust_method,p_adjust_scope Multiple-testing settings passed to
#'   [rc_test_condition_reactions()].
#' @param annotation_p Use `"p_adj"` or `"p_value"` for significance labels.
#' @param significance_threshold Maximum P value displayed as significant.
#' @param show_nonsignificant Display non-significant comparisons as `ns`.
#' @param show_omnibus Add the Kruskal-Wallis result to the subtitle when
#'   available.
#' @param point_size,point_alpha,jitter_width Metacell point controls.
#' @param box_width Boxplot width.
#' @param bracket_step Vertical spacing between significance brackets, expressed
#'   as a fraction of the observed score range.
#' @param title Optional plot title.
#' @param y_label Y-axis label.
#' @return A `ggplot` object. The underlying statistics, plotted metacells, and
#'   annotation table are attached as `condition_statistics`, `plot_data`, and
#'   `annotation_data` attributes.
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
    box_width = 0.55,
    bracket_step = 0.12,
    title = NULL,
    y_label = "Reaction support score") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package `ggplot2` is required for `rc_plot_condition_reaction()`.",
      call. = FALSE
    )
  }
  reaction_id <- .rc_plot_condition_single_string(
    reaction_id, "reaction_id"
  )
  cell_type <- .rc_plot_condition_single_string(cell_type, "cell_type")
  target_direction <- .rc_plot_condition_single_string(
    target_direction, "target_direction", allow_null = TRUE
  )
  medium_scenario <- .rc_plot_condition_single_string(
    medium_scenario, "medium_scenario", allow_null = TRUE
  )
  p_adjust_scope <- match.arg(p_adjust_scope)
  annotation_p <- match.arg(annotation_p)

  logical_controls <- c(
    show_nonsignificant = show_nonsignificant,
    show_omnibus = show_omnibus
  )
  if (anyNA(logical_controls) || !all(logical_controls %in% c(TRUE, FALSE))) {
    stop("Plot logical controls must be non-missing TRUE/FALSE values.",
         call. = FALSE)
  }
  numeric_controls <- c(
    significance_threshold = significance_threshold,
    point_size = point_size,
    point_alpha = point_alpha,
    jitter_width = jitter_width,
    box_width = box_width,
    bracket_step = bracket_step
  )
  if (any(!is.finite(numeric_controls)) ||
      significance_threshold <= 0 || significance_threshold > 1 ||
      point_size <= 0 || point_alpha <= 0 || point_alpha > 1 ||
      jitter_width < 0 || box_width <= 0 || bracket_step <= 0) {
    stop("Plot numeric controls contain invalid values.", call. = FALSE)
  }

  microcompass <- .rc_condition_stats_microcompass(x)
  penalty <- as.matrix(microcompass$penalty)
  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  row_meta$row_id <- rownames(penalty)
  target <- .rc_plot_condition_target(
    row_meta = row_meta,
    reaction_id = reaction_id,
    target_direction = target_direction,
    medium_scenario = medium_scenario
  )

  meta <- microcompass$unit_meta
  if (!is.data.frame(meta)) {
    stop("microCOMPASS `unit_meta` must be a data frame.", call. = FALSE)
  }
  condition_col <- .rc_condition_stats_column(
    meta,
    condition_col,
    c("condition", "dataset", "Group", "group", "treatment"),
    "condition_col"
  )
  celltype_col <- .rc_condition_stats_column(
    meta,
    celltype_col,
    c("cell_type", "celltype", "epithelial_or_stem", "CellType"),
    "celltype_col"
  )
  unit_id <- .rc_condition_stats_unit_ids(meta)
  if (!setequal(colnames(penalty), unit_id)) {
    stop("microCOMPASS matrices and unit metadata contain different units.",
         call. = FALSE)
  }
  meta$.unit_id <- unit_id
  meta <- meta[match(colnames(penalty), meta$.unit_id), , drop = FALSE]
  available_conditions <- unique(as.character(
    meta[[condition_col]][as.character(meta[[celltype_col]]) == cell_type]
  ))
  available_conditions <- available_conditions[
    !is.na(available_conditions) & nzchar(trimws(available_conditions))
  ]
  if (!length(available_conditions)) {
    stop("Cell type not found in Layer 2 metadata: ", cell_type,
         call. = FALSE)
  }
  conditions <- if (is.null(conditions)) {
    available_conditions
  } else {
    unique(as.character(conditions))
  }
  if (length(conditions) < 2L || anyNA(conditions) ||
      any(!nzchar(trimws(conditions)))) {
    stop("At least two non-empty conditions are required for plotting.",
         call. = FALSE)
  }
  unknown <- setdiff(conditions, available_conditions)
  if (length(unknown)) {
    stop(
      "Selected cell type lacks conditions: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  # Run the statistics over every scored reaction target in the selected cell
  # type and condition set. This keeps BH adjustment scoped across reactions;
  # filtering to one target before testing would collapse adjusted P to raw P.
  statistics <- rc_test_condition_reactions(
    x = microcompass,
    condition_col = condition_col,
    celltype_col = celltype_col,
    conditions = conditions,
    cell_types = cell_type,
    comparisons = comparisons,
    min_units = min_units,
    include_omnibus = TRUE,
    p_adjust_method = p_adjust_method,
    p_adjust_scope = p_adjust_scope,
    include_scores = TRUE
  )
  if (!target$row_id %in% rownames(statistics$score)) {
    stop("Selected target score was not returned by the statistics layer.",
         call. = FALSE)
  }

  score <- statistics$score[target$row_id, , drop = TRUE]
  plot_meta <- meta[match(names(score), meta$.unit_id), , drop = FALSE]
  plot_data <- data.frame(
    unit_id = names(score),
    condition = factor(
      as.character(plot_meta[[condition_col]]), levels = conditions
    ),
    cell_type = as.character(plot_meta[[celltype_col]]),
    score = as.numeric(score),
    stringsAsFactors = FALSE
  )
  plot_data <- plot_data[
    plot_data$cell_type == cell_type & is.finite(plot_data$score),
    , drop = FALSE
  ]
  if (!nrow(plot_data)) {
    stop("No finite metacell scores remain for the selected target.",
         call. = FALSE)
  }

  pairwise <- statistics$pairwise[
    statistics$pairwise$row_id == target$row_id &
      statistics$pairwise$cell_type == cell_type,
    , drop = FALSE
  ]
  score_min <- min(plot_data$score)
  score_max <- max(plot_data$score)
  annotation_data <- .rc_plot_condition_annotations(
    pairwise = pairwise,
    condition_levels = conditions,
    p_column = annotation_p,
    significance_threshold = significance_threshold,
    show_nonsignificant = show_nonsignificant,
    score_min = score_min,
    score_max = score_max,
    bracket_step = bracket_step
  )

  omnibus_subtitle <- NULL
  if (isTRUE(show_omnibus) && nrow(statistics$omnibus)) {
    omnibus <- statistics$omnibus[
      statistics$omnibus$row_id == target$row_id &
        statistics$omnibus$cell_type == cell_type,
      , drop = FALSE
    ]
    if (nrow(omnibus) == 1L && is.finite(omnibus$p_adj[[1L]])) {
      omnibus_subtitle <- paste0(
        "Kruskal-Wallis ", p_adjust_method, "-adjusted P = ",
        format.pval(omnibus$p_adj[[1L]], digits = 3, eps = 1e-4)
      )
    }
  }

  plot_title <- title %||% paste(
    target$reaction_id,
    paste0("(", target$target_direction, ")"),
    "in", cell_type
  )
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = condition, y = score, fill = condition)
  ) +
    ggplot2::geom_boxplot(
      width = box_width,
      outlier.shape = NA,
      alpha = 0.65,
      linewidth = 0.5
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = condition),
      width = jitter_width,
      height = 0,
      size = point_size,
      alpha = point_alpha,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = omnibus_subtitle,
      x = NULL,
      y = y_label,
      fill = "Condition"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 18, 8, 8)
    )

  if (nrow(annotation_data)) {
    for (i in seq_len(nrow(annotation_data))) {
      one <- annotation_data[i, , drop = FALSE]
      plot <- plot +
        ggplot2::annotate(
          "segment",
          x = one$xmin,
          xend = one$xmax,
          y = one$y,
          yend = one$y,
          linewidth = 0.45
        ) +
        ggplot2::annotate(
          "segment",
          x = one$xmin,
          xend = one$xmin,
          y = one$y - one$tip,
          yend = one$y,
          linewidth = 0.45
        ) +
        ggplot2::annotate(
          "segment",
          x = one$xmax,
          xend = one$xmax,
          y = one$y - one$tip,
          yend = one$y,
          linewidth = 0.45
        ) +
        ggplot2::annotate(
          "text",
          x = (one$xmin + one$xmax) / 2,
          y = one$text_y,
          label = one$label,
          vjust = 0,
          size = 4
        )
    }
    score_range <- score_max - score_min
    if (!is.finite(score_range) || score_range <= 0) {
      score_range <- max(1, abs(score_min), abs(score_max))
    }
    upper <- max(annotation_data$text_y) + score_range * 0.08
    lower <- score_min - score_range * 0.06
    plot <- plot + ggplot2::coord_cartesian(
      ylim = c(lower, upper), clip = "off"
    )
  } else {
    plot <- plot + ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.06, 0.10))
    )
  }

  attr(plot, "condition_statistics") <- statistics
  attr(plot, "plot_data") <- plot_data
  attr(plot, "annotation_data") <- annotation_data
  plot
}
