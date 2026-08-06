# Cell-type reaction ranking and visualization.

.rc_rank_nonempty_string <- function(value, name, allow_null = FALSE) {
  if (is.null(value) && isTRUE(allow_null)) return(NULL)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop("`", name, "` must be one non-empty string.", call. = FALSE)
  }
  trimws(value)
}

.rc_rank_optional_condition <- function(data, context) {
  if (!"condition" %in% colnames(data)) {
    return(list(
      available = FALSE,
      values = rep(NA_character_, nrow(data)),
      levels = character()
    ))
  }
  values <- trimws(as.character(data$condition))
  usable <- !is.na(values) & nzchar(values)
  if (any(usable) && any(!usable)) {
    stop(
      context,
      " contains a partially missing `condition` column. Supply complete labels or remove the column.",
      call. = FALSE
    )
  }
  if (!any(usable)) {
    return(list(
      available = FALSE,
      values = rep(NA_character_, nrow(data)),
      levels = character()
    ))
  }
  list(
    available = TRUE,
    values = values,
    levels = unique(values)
  )
}

.rc_rank_requested_conditions <- function(conditions, available) {
  if (is.null(conditions)) return(available)
  conditions <- unique(trimws(as.character(conditions)))
  if (!length(conditions) || anyNA(conditions) || any(!nzchar(conditions))) {
    stop("`conditions` must contain at least one non-empty label.",
         call. = FALSE)
  }
  unknown <- setdiff(conditions, available)
  if (length(unknown)) {
    stop(
      "Requested conditions were not found: ",
      paste(unknown, collapse = ", "),
      ". Available conditions: ", paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  conditions
}

.rc_rank_evidence_summary <- function(
    reaction_evidence, cell_type, selected_conditions) {
  evidence_selected <- reaction_evidence[
    trimws(as.character(reaction_evidence$cell_type)) == cell_type,
    , drop = FALSE
  ]
  condition <- .rc_rank_optional_condition(
    evidence_selected, "`reaction_evidence`"
  )
  if (condition$available && !is.null(selected_conditions)) {
    evidence_selected <- evidence_selected[
      condition$values %in% selected_conditions,
      , drop = FALSE
    ]
  }
  rows <- lapply(
    split(
      evidence_selected,
      as.character(evidence_selected$reaction_id),
      drop = TRUE
    ),
    function(one) {
      evidence_class <- as.character(one$evidence_class)
      active_multiome <- any(
        one$has_active_multiome_contribution %in% TRUE
      )
      data.frame(
        reaction_id = as.character(one$reaction_id[[1L]]),
        evidence_eligible = any(
          evidence_class %in% c("RNA-only", "RNA+ATAC")
        ),
        support_class = if (active_multiome) {
          "Multiome-supported"
        } else {
          "RNA-only"
        },
        stringsAsFactors = FALSE
      )
    }
  )
  if (!length(rows)) {
    return(data.frame(
      reaction_id = character(),
      evidence_eligible = logical(),
      support_class = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Plot the top reaction targets within one cell type
#'
#' Ranks directional reaction targets by the median metacell support score
#' `-log(penalty_per_target_flux + 1e-8)` within one cell type and medium.
#' Reactions are eligible only when their formal evidence class is `RNA-only`
#' or `RNA+ATAC`. A reaction is labelled `Multiome-supported` when at least one
#' selected evidence row has an active ATAC contribution that changes the
#' GPR-aggregated reaction capacity.
#'
#' Condition metadata are optional. When `condition` is absent or entirely
#' empty in the annotated result tables, all matching metacells are pooled and
#' `conditions` must remain `NULL`. A single condition is valid and does not
#' require a contrast. When condition labels are available, `conditions = NULL`
#' uses every represented condition.
#'
#' @param result A complete RegCompass result containing
#'   `reaction_comparison_by_metacell`, `reaction_catalog`, and
#'   `reaction_evidence`.
#' @param cell_type One cell-type label.
#' @param target_direction `"forward"` or `"reverse"`.
#' @param medium_scenario Optional medium label. It may be omitted when exactly
#'   one medium remains after cell-type, direction, and condition filtering.
#' @param conditions Optional condition labels. Leave `NULL` to pool all
#'   available conditions, or when condition metadata were not supplied.
#' @param top_n Maximum number of reaction targets to display.
#' @param label_width Approximate character width used to wrap reaction labels.
#' @return A `ggplot` object. Ranked data and resolved selection are attached as
#'   `rank_data` and `selection` attributes.
#' @export
plot_top_celltype_reaction_rank <- function(
    result,
    cell_type,
    target_direction = c("forward", "reverse"),
    medium_scenario = NULL,
    conditions = NULL,
    top_n = 20L,
    label_width = 58L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package `ggplot2` is required for `plot_top_celltype_reaction_rank()`.",
      call. = FALSE
    )
  }
  target_direction <- match.arg(target_direction)
  cell_type <- .rc_rank_nonempty_string(cell_type, "cell_type")
  medium_scenario <- .rc_rank_nonempty_string(
    medium_scenario, "medium_scenario", allow_null = TRUE
  )
  if (length(top_n) != 1L || is.na(top_n) || !is.finite(top_n) ||
      top_n < 1 || top_n != as.integer(top_n)) {
    stop("`top_n` must be one positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)
  if (length(label_width) != 1L || is.na(label_width) ||
      !is.finite(label_width) || label_width < 20 ||
      label_width != as.integer(label_width)) {
    stop("`label_width` must be one integer of at least 20.", call. = FALSE)
  }
  label_width <- as.integer(label_width)

  reaction_long <- result$reaction_comparison_by_metacell
  reaction_catalog <- result$reaction_catalog
  reaction_evidence <- result$reaction_evidence
  required_long <- c(
    "reaction_id", "direction", "medium", "cell_type", "metacell_id",
    "penalty_available", "penalty_per_target_flux"
  )
  required_catalog <- c(
    "reaction_id", "reaction_name", "forward_formula", "reverse_formula"
  )
  required_evidence <- c(
    "reaction_id", "cell_type", "evidence_class",
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
  if (!nrow(selected)) {
    stop("No scored metacells match the requested cell type and direction.",
         call. = FALSE)
  }

  condition <- .rc_rank_optional_condition(
    selected, "`reaction_comparison_by_metacell`"
  )
  if (!condition$available && !is.null(conditions)) {
    stop(
      "`conditions` cannot be supplied because the result has no usable condition metadata.",
      call. = FALSE
    )
  }
  selected_conditions <- if (condition$available) {
    .rc_rank_requested_conditions(conditions, condition$levels)
  } else {
    NULL
  }
  if (condition$available) {
    selected <- selected[
      condition$values %in% selected_conditions,
      , drop = FALSE
    ]
  }

  available_media <- unique(trimws(as.character(selected$medium)))
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
  selected <- selected[
    trimws(as.character(selected$medium)) == medium_scenario,
    , drop = FALSE
  ]
  if (!nrow(selected)) {
    stop("No scored metacells match the requested medium.", call. = FALSE)
  }

  if (condition$available) {
    selected_conditions <- unique(trimws(as.character(selected$condition)))
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
      condition_count <- if (condition$available) {
        length(unique(trimws(as.character(one$condition[finite]))))
      } else {
        0L
      }
      data.frame(
        reaction_id = as.character(one$reaction_id[[1L]]),
        median_support_score = stats::median(one$support_score[finite]),
        mean_support_score = mean(one$support_score[finite]),
        median_penalty_per_target_flux = stats::median(
          one$penalty_per_target_flux[finite]
        ),
        n_metacells = length(unique(as.character(one$metacell_id[finite]))),
        n_conditions = as.integer(condition_count),
        stringsAsFactors = FALSE
      )
    }
  )
  scores <- do.call(rbind, score_rows)
  rownames(scores) <- NULL

  evidence_summary <- .rc_rank_evidence_summary(
    reaction_evidence = reaction_evidence,
    cell_type = cell_type,
    selected_conditions = selected_conditions
  )
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
  ranked$reaction_label <- vapply(
    ranked$reaction_label_text,
    function(x) paste(strwrap(x, width = label_width), collapse = "\n"),
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

  condition_text <- if (condition$available) {
    paste(selected_conditions, collapse = ", ")
  } else {
    "all metacells; no condition grouping"
  }
  plot <- ggplot2::ggplot(
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
        "\nMultiome-supported if any selected evidence changes GPR reaction capacity"
      )
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      plot.title.position = "plot",
      axis.text.y = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5)
    )

  attr(plot, "rank_data") <- transform(
    ranked,
    reaction_label = as.character(reaction_label),
    support_class = as.character(support_class)
  )
  attr(plot, "selection") <- list(
    cell_type = cell_type,
    target_direction = target_direction,
    medium_scenario = medium_scenario,
    conditions = selected_conditions,
    condition_available = condition$available,
    top_n = top_n,
    ranking_statistic = "median_support_score",
    multiome_rule =
      "any(has_active_multiome_contribution) across selected evidence rows"
  )
  plot
}
