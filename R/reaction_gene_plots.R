# Gene-centered collections of condition-reaction plots.

.rc_ra_plot_one_precomputed <- function(
    x, statistics, row_id, cell_type, condition_col, celltype_col,
    conditions, annotation_p, significance_threshold,
    show_nonsignificant, show_omnibus, point_size, point_alpha,
    jitter_width, box_width, bracket_step) {
  microcompass <- .rc_condition_stats_microcompass(x)
  meta <- microcompass$unit_meta
  condition_col <- .rc_condition_stats_column(
    meta, condition_col,
    c("condition", "dataset", "Group", "group", "treatment"),
    "condition_col"
  )
  celltype_col <- .rc_condition_stats_column(
    meta, celltype_col,
    c("cell_type", "celltype", "epithelial_or_stem", "CellType"),
    "celltype_col"
  )
  unit_ids <- .rc_condition_stats_unit_ids(meta)
  meta$.unit_id <- unit_ids
  score <- statistics$score[row_id, , drop = TRUE]
  plot_meta <- meta[match(names(score), meta$.unit_id), , drop = FALSE]
  plot_data <- data.frame(
    unit_id = names(score),
    condition = factor(as.character(plot_meta[[condition_col]]), levels = conditions),
    cell_type = as.character(plot_meta[[celltype_col]]),
    score = as.numeric(score),
    stringsAsFactors = FALSE
  )
  plot_data <- plot_data[
    plot_data$cell_type == cell_type & is.finite(plot_data$score), , drop = FALSE
  ]
  pairwise <- statistics$pairwise[
    statistics$pairwise$row_id == row_id &
      statistics$pairwise$cell_type == cell_type, , drop = FALSE
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
  target <- pairwise[1L, , drop = FALSE]
  omnibus_subtitle <- NULL
  if (isTRUE(show_omnibus) && nrow(statistics$omnibus)) {
    omnibus <- statistics$omnibus[
      statistics$omnibus$row_id == row_id &
        statistics$omnibus$cell_type == cell_type, , drop = FALSE
    ]
    if (nrow(omnibus) == 1L && is.finite(omnibus$p_adj[[1L]])) {
      omnibus_subtitle <- paste0(
        "Kruskal-Wallis ", statistics$params$p_adjust_method,
        "-adjusted P = ",
        format.pval(omnibus$p_adj[[1L]], digits = 3, eps = 1e-4)
      )
    }
  }
  plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = condition, y = score, fill = condition)
  ) +
    ggplot2::geom_boxplot(
      width = box_width, outlier.shape = NA, alpha = 0.65, linewidth = 0.5
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = condition), width = jitter_width, height = 0,
      size = point_size, alpha = point_alpha, show.legend = FALSE
    ) +
    ggplot2::labs(
      title = paste0(
        target$reaction_name[[1L]], " (", target$reaction_id[[1L]], ", ",
        target$target_direction[[1L]], ") in ", cell_type
      ),
      subtitle = omnibus_subtitle,
      caption = .rc_ra_plot_caption(
        target, target$evidence_comparison[[1L]] %||% NULL
      ),
      x = NULL, y = "Reaction support score", fill = "Condition"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 18, 8, 8),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  if (nrow(annotation_data)) {
    for (i in seq_len(nrow(annotation_data))) {
      one <- annotation_data[i, , drop = FALSE]
      plot <- plot +
        ggplot2::annotate(
          "segment", x = one$xmin, xend = one$xmax,
          y = one$y, yend = one$y, linewidth = 0.45
        ) +
        ggplot2::annotate(
          "segment", x = one$xmin, xend = one$xmin,
          y = one$y - one$tip, yend = one$y, linewidth = 0.45
        ) +
        ggplot2::annotate(
          "segment", x = one$xmax, xend = one$xmax,
          y = one$y - one$tip, yend = one$y, linewidth = 0.45
        ) +
        ggplot2::annotate(
          "text", x = (one$xmin + one$xmax) / 2,
          y = one$text_y, label = one$label, vjust = 0, size = 4
        )
    }
    score_range <- score_max - score_min
    if (!is.finite(score_range) || score_range <= 0) {
      score_range <- max(1, abs(score_min), abs(score_max))
    }
    plot <- plot + ggplot2::coord_cartesian(
      ylim = c(
        score_min - score_range * 0.06,
        max(annotation_data$text_y) + score_range * 0.08
      ),
      clip = "off"
    )
  }
  attr(plot, "plot_data") <- plot_data
  attr(plot, "annotation_data") <- annotation_data
  attr(plot, "reaction_annotation") <- target
  plot
}

#' Plot significant condition responses for reactions selected by genes
#'
#' Runs the condition statistics once, selects scored reactions containing the
#' requested metabolic genes, ranks significant reaction-direction targets, and
#' returns a named collection of annotated boxplots.
#'
#' @param x Annotated RegCompass result.
#' @param genes Metabolic gene symbols used to select GPR reactions.
#' @param cell_type One cell type.
#' @param condition_col,celltype_col Metadata columns.
#' @param conditions Ordered conditions.
#' @param comparisons Optional condition pairs.
#' @param target_directions Optional scored directions.
#' @param medium_scenario Optional medium.
#' @param evidence_class Optional group evidence classes used for gene-reaction
#'   selection.
#' @param p_adj_max Maximum adjusted pairwise P value.
#' @param min_abs_rank_biserial Minimum absolute rank-biserial effect.
#' @param max_reactions Maximum number of reaction-direction plots.
#' @param annotation_p Significance label column.
#' @param outdir Optional directory for PDF plots and selection tables.
#' @return A `regcompass_gene_reaction_plots` list.
#' @export
rc_plot_condition_gene_reactions <- function(
    x, genes, cell_type,
    condition_col = NULL, celltype_col = NULL,
    conditions = NULL, comparisons = NULL,
    target_directions = NULL, medium_scenario = NULL,
    evidence_class = NULL,
    p_adj_max = 0.05, min_abs_rank_biserial = 0.30,
    max_reactions = 12L,
    annotation_p = c("p_adj", "p_value"),
    significance_threshold = 0.05,
    show_nonsignificant = FALSE, show_omnibus = TRUE,
    point_size = 1.8, point_alpha = 0.75,
    jitter_width = 0.12, box_width = 0.55, bracket_step = 0.12,
    outdir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required.", call. = FALSE)
  }
  annotation_p <- match.arg(annotation_p)
  if (!is.numeric(p_adj_max) || length(p_adj_max) != 1L ||
      !is.finite(p_adj_max) || p_adj_max <= 0 || p_adj_max > 1 ||
      !is.numeric(min_abs_rank_biserial) ||
      length(min_abs_rank_biserial) != 1L ||
      !is.finite(min_abs_rank_biserial) || min_abs_rank_biserial < 0 ||
      min_abs_rank_biserial > 1 ||
      !is.numeric(max_reactions) || length(max_reactions) != 1L ||
      !is.finite(max_reactions) || max_reactions < 1) {
    stop("Gene-reaction plot selection thresholds are invalid.", call. = FALSE)
  }
  max_reactions <- as.integer(max_reactions)
  selection <- rc_select_gene_reactions(
    x = x,
    genes = genes,
    cell_types = cell_type,
    evidence_class = evidence_class
  )
  if (!length(selection$reaction_ids)) {
    stop("No scored reactions matched the requested genes and evidence filters.",
         call. = FALSE)
  }
  statistics <- rc_test_condition_reactions(
    x = x,
    condition_col = condition_col,
    celltype_col = celltype_col,
    conditions = conditions,
    cell_types = cell_type,
    comparisons = comparisons,
    min_units = 5L,
    include_omnibus = TRUE,
    p_adjust_method = "BH",
    p_adjust_scope = "celltype_contrast_medium",
    include_scores = TRUE
  )
  pairwise <- statistics$pairwise[
    statistics$pairwise$reaction_id %in% selection$reaction_ids &
      statistics$pairwise$cell_type == cell_type &
      is.finite(statistics$pairwise$p_adj) &
      statistics$pairwise$p_adj <= p_adj_max &
      is.finite(statistics$pairwise$rank_biserial_b_minus_a) &
      abs(statistics$pairwise$rank_biserial_b_minus_a) >=
        min_abs_rank_biserial,
    , drop = FALSE
  ]
  if (!is.null(target_directions)) {
    pairwise <- pairwise[
      pairwise$target_direction %in% as.character(target_directions), , drop = FALSE
    ]
  }
  if (!is.null(medium_scenario)) {
    pairwise <- pairwise[
      pairwise$medium_scenario %in% as.character(medium_scenario), , drop = FALSE
    ]
  }
  if (!nrow(pairwise)) {
    stop("No gene-associated targets passed the significance and effect filters.",
         call. = FALSE)
  }
  split_rows <- split(seq_len(nrow(pairwise)), pairwise$row_id)
  ranked <- do.call(rbind, lapply(split_rows, function(rows) {
    one <- pairwise[rows, , drop = FALSE]
    data.frame(
      row_id = one$row_id[[1L]],
      reaction_id = one$reaction_id[[1L]],
      reaction_name = one$reaction_name[[1L]],
      target_direction = one$target_direction[[1L]],
      medium_scenario = one$medium_scenario[[1L]],
      tested_formula = one$tested_formula[[1L]],
      genes = one$genes[[1L]],
      min_p_adj = min(one$p_adj, na.rm = TRUE),
      max_abs_rank_biserial = max(
        abs(one$rank_biserial_b_minus_a), na.rm = TRUE
      ),
      max_abs_delta_median = max(
        abs(one$delta_median_score_b_minus_a), na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
  ranked <- ranked[order(
    ranked$min_p_adj,
    -ranked$max_abs_rank_biserial,
    -ranked$max_abs_delta_median,
    ranked$reaction_id,
    ranked$target_direction
  ), , drop = FALSE]
  ranked <- utils::head(ranked, max_reactions)
  if (is.null(conditions)) conditions <- statistics$params$conditions
  plots <- lapply(ranked$row_id, function(row_id) {
    .rc_ra_plot_one_precomputed(
      x = x,
      statistics = statistics,
      row_id = row_id,
      cell_type = cell_type,
      condition_col = condition_col,
      celltype_col = celltype_col,
      conditions = conditions,
      annotation_p = annotation_p,
      significance_threshold = significance_threshold,
      show_nonsignificant = show_nonsignificant,
      show_omnibus = show_omnibus,
      point_size = point_size,
      point_alpha = point_alpha,
      jitter_width = jitter_width,
      box_width = box_width,
      bracket_step = bracket_step
    )
  })
  names(plots) <- ranked$row_id
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    .rc_write_tsv_gz(
      ranked,
      file.path(outdir, "selected_gene_reaction_targets.tsv.gz")
    )
    .rc_write_tsv_gz(
      pairwise,
      file.path(outdir, "selected_gene_reaction_pairwise.tsv.gz")
    )
    for (i in seq_along(plots)) {
      safe <- gsub("[^A-Za-z0-9_.-]+", "_", ranked$row_id[[i]])
      ggplot2::ggsave(
        filename = file.path(outdir, paste0(safe, ".pdf")),
        plot = plots[[i]], width = 7, height = 5.5
      )
    }
  }
  answer <- list(
    plots = plots,
    selected_targets = ranked,
    pairwise_hits = pairwise,
    statistics = statistics,
    gene_selection = selection
  )
  class(answer) <- c("regcompass_gene_reaction_plots", "list")
  answer
}
