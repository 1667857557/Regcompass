# Direction-aware condition reporting for reversible reaction targets.

.rc_direction_report_source_label <- function(source_label) {
  if (is.null(source_label)) return(NA_character_)
  valid <- is.character(source_label) && length(source_label) == 1L &&
    !is.na(source_label) && nzchar(trimws(source_label))
  if (!valid) {
    stop("`source_label` must be NULL or one non-empty string.", call. = FALSE)
  }
  trimws(source_label)
}

.rc_direction_report_unit_table <- function(
    statistics, microcompass, condition_col, celltype_col,
    direction_tolerance, source_label) {
  score <- as.matrix(statistics$score)
  row_meta <- rc_parse_microcompass_row_id(rownames(score))
  row_meta$row_id <- rownames(score)

  meta <- microcompass$unit_meta
  unit_id <- .rc_condition_stats_unit_ids(meta)
  meta$.unit_id <- unit_id
  meta <- meta[match(colnames(score), meta$.unit_id), , drop = FALSE]
  if (anyNA(meta$.unit_id)) {
    stop(
      "Condition-report unit metadata do not align with score columns.",
      call. = FALSE
    )
  }

  condition <- as.character(meta[[condition_col]])
  cell_type <- as.character(meta[[celltype_col]])
  keep_unit <- condition %in% statistics$params$conditions &
    cell_type %in% statistics$params$cell_types
  score <- score[, keep_unit, drop = FALSE]
  condition <- condition[keep_unit]
  cell_type <- cell_type[keep_unit]

  key <- paste(row_meta$reaction_id, row_meta$medium_scenario, sep = "\001")
  groups <- split(seq_len(nrow(row_meta)), key)
  rows <- lapply(groups, function(index) {
    one <- row_meta[index, , drop = FALSE]
    forward <- index[one$target_direction == "forward"]
    reverse <- index[one$target_direction == "reverse"]
    if (length(forward) > 1L || length(reverse) > 1L) {
      stop(
        paste(
          "Direction report requires at most one forward and one reverse",
          "target per reaction and medium."
        ),
        call. = FALSE
      )
    }

    score_forward <- rep(NA_real_, ncol(score))
    score_reverse <- rep(NA_real_, ncol(score))
    if (length(forward)) score_forward <- score[forward, , drop = TRUE]
    if (length(reverse)) score_reverse <- score[reverse, , drop = TRUE]

    finite_forward <- is.finite(score_forward)
    finite_reverse <- is.finite(score_reverse)
    any_direction <- rep(NA_real_, ncol(score))
    both_finite <- finite_forward & finite_reverse
    any_direction[both_finite] <- pmax(
      score_forward[both_finite], score_reverse[both_finite]
    )
    only_forward <- finite_forward & !finite_reverse
    only_reverse <- finite_reverse & !finite_forward
    any_direction[only_forward] <- score_forward[only_forward]
    any_direction[only_reverse] <- score_reverse[only_reverse]

    directional_balance <- rep(NA_real_, ncol(score))
    directional_balance[both_finite] <-
      score_forward[both_finite] - score_reverse[both_finite]

    data.frame(
      unit_id = colnames(score),
      reaction_id = one$reaction_id[[1L]],
      medium_scenario = one$medium_scenario[[1L]],
      condition = condition,
      cell_type = cell_type,
      forward_target_available = length(forward) == 1L,
      reverse_target_available = length(reverse) == 1L,
      score_forward = score_forward,
      score_reverse = score_reverse,
      any_direction_support = any_direction,
      directional_balance = directional_balance,
      source_label = source_label,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "direction_tolerance") <- direction_tolerance
  out
}

.rc_direction_report_diagnostics <- function(unit_table, direction_tolerance) {
  key <- interaction(
    unit_table$reaction_id,
    unit_table$medium_scenario,
    unit_table$cell_type,
    unit_table$condition,
    drop = TRUE,
    lex.order = TRUE
  )
  groups <- split(seq_len(nrow(unit_table)), key)
  rows <- lapply(groups, function(index) {
    one <- unit_table[index, , drop = FALSE]
    forward_available <- isTRUE(one$forward_target_available[[1L]])
    reverse_available <- isTRUE(one$reverse_target_available[[1L]])
    finite_forward <- is.finite(one$score_forward)
    finite_reverse <- is.finite(one$score_reverse)
    paired <- finite_forward & finite_reverse
    same_missing_pattern <- if (forward_available && reverse_available) {
      identical(finite_forward, finite_reverse)
    } else {
      NA
    }
    max_abs_difference <- if (any(paired)) {
      max(abs(one$score_forward[paired] - one$score_reverse[paired]))
    } else {
      NA_real_
    }
    indistinguishable <- forward_available && reverse_available &&
      isTRUE(same_missing_pattern) && any(paired) &&
      is.finite(max_abs_difference) &&
      max_abs_difference <= direction_tolerance

    median_forward <- if (any(finite_forward)) {
      stats::median(one$score_forward[finite_forward])
    } else {
      NA_real_
    }
    median_reverse <- if (any(finite_reverse)) {
      stats::median(one$score_reverse[finite_reverse])
    } else {
      NA_real_
    }
    median_balance <- if (any(paired)) {
      stats::median(one$directional_balance[paired])
    } else {
      NA_real_
    }

    pair_status <- if (forward_available && reverse_available) {
      if (!any(paired)) {
        "bidirectional_no_paired_scores"
      } else if (indistinguishable) {
        "bidirectional_indistinguishable"
      } else {
        "bidirectional_distinguishable"
      }
    } else if (forward_available) {
      "forward_only"
    } else if (reverse_available) {
      "reverse_only"
    } else {
      "no_scored_direction"
    }

    preferred_direction <- if (identical(pair_status, "forward_only")) {
      "forward"
    } else if (identical(pair_status, "reverse_only")) {
      "reverse"
    } else if (identical(pair_status, "bidirectional_indistinguishable")) {
      "indistinguishable"
    } else if (!is.finite(median_balance)) {
      NA_character_
    } else if (median_balance > direction_tolerance) {
      "forward"
    } else if (median_balance < -direction_tolerance) {
      "reverse"
    } else {
      "balanced"
    }

    data.frame(
      reaction_id = one$reaction_id[[1L]],
      medium_scenario = one$medium_scenario[[1L]],
      cell_type = one$cell_type[[1L]],
      condition = one$condition[[1L]],
      n_units = nrow(one),
      n_forward_scores = sum(finite_forward),
      n_reverse_scores = sum(finite_reverse),
      n_paired_direction_scores = sum(paired),
      median_forward_support = median_forward,
      median_reverse_support = median_reverse,
      median_any_direction_support = if (
        any(is.finite(one$any_direction_support))
      ) {
        stats::median(one$any_direction_support, na.rm = TRUE)
      } else {
        NA_real_
      },
      median_directional_balance = median_balance,
      same_missing_pattern = same_missing_pattern,
      max_abs_forward_reverse_difference = max_abs_difference,
      directionally_indistinguishable = indistinguishable,
      direction_pair_status = pair_status,
      preferred_direction = preferred_direction,
      source_label = one$source_label[[1L]],
      interpretation = paste(
        "Forward and reverse are independent counterfactual LP targets;",
        "directional balance is not net flux."
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.rc_direction_report_metric_table <- function(unit_table) {
  any_direction <- data.frame(
    unit_id = unit_table$unit_id,
    reaction_id = unit_table$reaction_id,
    medium_scenario = unit_table$medium_scenario,
    condition = unit_table$condition,
    cell_type = unit_table$cell_type,
    report_metric = "any_direction_support",
    metric_value = unit_table$any_direction_support,
    source_label = unit_table$source_label,
    stringsAsFactors = FALSE
  )
  balance_keep <- unit_table$forward_target_available &
    unit_table$reverse_target_available
  directional_balance <- data.frame(
    unit_id = unit_table$unit_id[balance_keep],
    reaction_id = unit_table$reaction_id[balance_keep],
    medium_scenario = unit_table$medium_scenario[balance_keep],
    condition = unit_table$condition[balance_keep],
    cell_type = unit_table$cell_type[balance_keep],
    report_metric = "directional_balance",
    metric_value = unit_table$directional_balance[balance_keep],
    source_label = unit_table$source_label[balance_keep],
    stringsAsFactors = FALSE
  )
  rbind(any_direction, directional_balance)
}

.rc_direction_report_adjust_groups <- function(scope, omnibus = FALSE) {
  if (isTRUE(omnibus)) {
    return(switch(
      scope,
      celltype_contrast_medium = c(
        "report_metric", "cell_type", "medium_scenario"
      ),
      celltype_contrast = c("report_metric", "cell_type"),
      celltype = c("report_metric", "cell_type"),
      global = "report_metric"
    ))
  }
  switch(
    scope,
    celltype_contrast_medium = c(
      "report_metric", "cell_type", "condition_a", "condition_b",
      "medium_scenario"
    ),
    celltype_contrast = c(
      "report_metric", "cell_type", "condition_a", "condition_b"
    ),
    celltype = c("report_metric", "cell_type"),
    global = "report_metric"
  )
}

.rc_direction_report_pairwise <- function(
    metric_table, comparisons, min_units, wilcox_correct,
    p_adjust_method, p_adjust_scope, analysis_unit,
    biological_replicate_inference, direction_tolerance) {
  group_key <- interaction(
    metric_table$reaction_id,
    metric_table$medium_scenario,
    metric_table$cell_type,
    metric_table$report_metric,
    drop = TRUE,
    lex.order = TRUE
  )
  groups <- split(seq_len(nrow(metric_table)), group_key)
  rows <- list()
  index <- 0L
  for (group_rows in groups) {
    one <- metric_table[group_rows, , drop = FALSE]
    for (pair in comparisons) {
      condition_a <- pair[[1L]]
      condition_b <- pair[[2L]]
      a <- one$metric_value[one$condition == condition_a]
      b <- one$metric_value[one$condition == condition_b]
      effect <- .rc_condition_stats_effect(
        a,
        b,
        min_units = min_units,
        wilcox_correct = wilcox_correct
      )
      delta <- effect$delta_median_score_b_minus_a
      index <- index + 1L
      rows[[index]] <- data.frame(
        reaction_id = one$reaction_id[[1L]],
        medium_scenario = one$medium_scenario[[1L]],
        cell_type = one$cell_type[[1L]],
        report_metric = one$report_metric[[1L]],
        condition_a = condition_a,
        condition_b = condition_b,
        n_a = effect$n_a,
        n_b = effect$n_b,
        median_metric_a = effect$median_score_a,
        median_metric_b = effect$median_score_b,
        delta_median_metric_b_minus_a = delta,
        mean_metric_a = effect$mean_score_a,
        mean_metric_b = effect$mean_score_b,
        cohens_d_b_minus_a = effect$cohens_d_b_minus_a,
        rank_biserial_b_minus_a = effect$rank_biserial_b_minus_a,
        common_language_b_greater_a = effect$common_language_b_greater_a,
        p_value = effect$p_value,
        test_status = effect$test_status,
        higher_metric_condition = if (is.finite(delta)) {
          if (delta > 0) condition_b else if (delta < 0) condition_a else "tie"
        } else {
          NA_character_
        },
        direction_shift_b_minus_a = if (
          identical(one$report_metric[[1L]], "directional_balance") &&
          is.finite(delta)
        ) {
          if (delta > direction_tolerance) {
            "toward_forward"
          } else if (delta < -direction_tolerance) {
            "toward_reverse"
          } else {
            "no_shift"
          }
        } else {
          NA_character_
        },
        metric_definition = if (
          identical(one$report_metric[[1L]], "any_direction_support")
        ) {
          paste(
            "max(forward_support, reverse_support);",
            "one available direction is retained"
          )
        } else {
          "forward_support - reverse_support"
        },
        metric_semantics = if (
          identical(one$report_metric[[1L]], "any_direction_support")
        ) {
          "best-supported direction; non-additive reaction potential"
        } else {
          "directional support asymmetry; not net flux"
        },
        analysis_unit = analysis_unit,
        inference_level = if (biological_replicate_inference) {
          "biological_sample_celltype"
        } else if (identical(analysis_unit, "metacell")) {
          "metacell_within_dataset"
        } else {
          paste0(analysis_unit, "_within_dataset")
        },
        descriptive_only = !biological_replicate_inference,
        biological_replicate_inference = biological_replicate_inference,
        source_label = one$source_label[[1L]],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out <- .rc_condition_stats_adjust(
    out,
    .rc_direction_report_adjust_groups(p_adjust_scope, omnibus = FALSE),
    p_adjust_method
  )
  out$p_adjust_method <- p_adjust_method
  out$p_adjust_scope <- p_adjust_scope
  out <- out[order(
    out$cell_type,
    out$condition_a,
    out$condition_b,
    out$medium_scenario,
    out$report_metric,
    out$p_adj,
    -abs(out$rank_biserial_b_minus_a),
    out$reaction_id,
    na.last = TRUE
  ), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.rc_direction_report_omnibus <- function(
    metric_table, conditions, min_units, p_adjust_method, p_adjust_scope,
    analysis_unit, biological_replicate_inference) {
  group_key <- interaction(
    metric_table$reaction_id,
    metric_table$medium_scenario,
    metric_table$cell_type,
    metric_table$report_metric,
    drop = TRUE,
    lex.order = TRUE
  )
  groups_split <- split(seq_len(nrow(metric_table)), group_key)
  rows <- list()
  index <- 0L
  for (group_rows in groups_split) {
    one <- metric_table[group_rows, , drop = FALSE]
    groups <- lapply(conditions, function(condition) {
      values <- one$metric_value[one$condition == condition]
      values[is.finite(values)]
    })
    names(groups) <- conditions
    counts <- vapply(groups, length, integer(1))
    p_value <- NA_real_
    status <- "insufficient_units"
    if (all(counts >= min_units)) {
      values <- unlist(groups, use.names = FALSE)
      labels <- factor(
        rep(names(groups), lengths(groups)),
        levels = names(groups)
      )
      if (length(unique(values)) == 1L) {
        p_value <- 1
        status <- "constant_equal"
      } else {
        test <- suppressWarnings(stats::kruskal.test(x = values, g = labels))
        p_value <- unname(test$p.value)
        status <- if (is.finite(p_value)) "ok" else "test_failed"
      }
    }
    index <- index + 1L
    rows[[index]] <- data.frame(
      reaction_id = one$reaction_id[[1L]],
      medium_scenario = one$medium_scenario[[1L]],
      cell_type = one$cell_type[[1L]],
      report_metric = one$report_metric[[1L]],
      n_conditions = length(conditions),
      n_units_total = sum(counts),
      min_units_per_condition = min(counts),
      max_units_per_condition = max(counts),
      units_by_condition = paste(
        paste(names(counts), counts, sep = "="),
        collapse = ";"
      ),
      p_value = p_value,
      test_status = status,
      metric_definition = if (
        identical(one$report_metric[[1L]], "any_direction_support")
      ) {
        paste(
          "max(forward_support, reverse_support);",
          "one available direction is retained"
        )
      } else {
        "forward_support - reverse_support"
      },
      metric_semantics = if (
        identical(one$report_metric[[1L]], "any_direction_support")
      ) {
        "best-supported direction; non-additive reaction potential"
      } else {
        "directional support asymmetry; not net flux"
      },
      analysis_unit = analysis_unit,
      inference_level = if (biological_replicate_inference) {
        "biological_sample_celltype"
      } else if (identical(analysis_unit, "metacell")) {
        "metacell_within_dataset"
      } else {
        paste0(analysis_unit, "_within_dataset")
      },
      descriptive_only = !biological_replicate_inference,
      biological_replicate_inference = biological_replicate_inference,
      source_label = one$source_label[[1L]],
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out <- .rc_condition_stats_adjust(
    out,
    .rc_direction_report_adjust_groups(p_adjust_scope, omnibus = TRUE),
    p_adjust_method
  )
  out$p_adjust_method <- p_adjust_method
  out$p_adjust_scope <- p_adjust_scope
  out <- out[order(
    out$cell_type,
    out$medium_scenario,
    out$report_metric,
    out$p_adj,
    out$reaction_id,
    na.last = TRUE
  ), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.rc_direction_report_static_annotations <- function(statistics) {
  source <- if (nrow(statistics$pairwise)) {
    statistics$pairwise
  } else {
    statistics$omnibus
  }
  columns <- intersect(
    c(
      "reaction_id", "reaction_name", "subsystem", "reaction_role",
      "lower_bound", "upper_bound", "reversible", "model_formula",
      "forward_substrates", "forward_products", "forward_formula",
      "reverse_substrates", "reverse_products", "reverse_formula",
      "genes", "gpr_rule", "n_gpr_genes", "has_gpr",
      "kegg_reaction_id", "reactome_reaction_id", "rhea_reaction_id",
      "rhea_master_id"
    ),
    colnames(source)
  )
  if (!length(columns) || !"reaction_id" %in% columns) return(NULL)
  annotation <- unique(source[, columns, drop = FALSE])
  annotation <- annotation[
    !duplicated(annotation$reaction_id),
    ,
    drop = FALSE
  ]
  rownames(annotation) <- NULL
  annotation
}

.rc_direction_report_attach_annotations <- function(data, annotation) {
  if (is.null(annotation) || !nrow(data)) return(data)
  index <- match(
    as.character(data$reaction_id),
    as.character(annotation$reaction_id)
  )
  extra <- setdiff(colnames(annotation), "reaction_id")
  cbind(data, annotation[index, extra, drop = FALSE])
}

#' Build a direction-aware condition report
#'
#' Runs [rc_test_condition_reactions()] with both allowed LP directions and
#' preserves the COMPASS-style direction-specific results. It additionally
#' derives two non-additive reaction-level summaries: `any_direction_support`
#' is the larger available forward/reverse support score, and
#' `directional_balance` is forward support minus reverse support. The latter
#' describes support asymmetry and must not be interpreted as net flux.
#'
#' @param x A microCOMPASS result or RegCompass result containing one.
#' @param condition_col,celltype_col,conditions,cell_types,comparisons,
#'   reaction_ids,medium_scenarios,min_units,include_omnibus,p_adjust_method,
#'   p_adjust_scope,wilcox_correct,eps,vmax_tolerance Arguments passed to
#'   [rc_test_condition_reactions()]. Both allowed target directions are always
#'   requested.
#' @param direction_tolerance Non-negative tolerance used to classify paired
#'   forward/reverse unit scores as numerically indistinguishable.
#' @param include_unit_metrics Include the unit-level forward, reverse,
#'   any-direction, and balance table in the returned object.
#' @param source_label Optional provenance label, for example
#'   `"original_layer2_core"` or `"target_union_second_pass_noncore"`.
#' @param outdir Optional output directory.
#' @return A `regcompass_condition_direction_report` list containing the
#'   original `directional_pairwise` and `directional_omnibus` tables,
#'   direction-collapsed `reaction_pairwise` and `reaction_omnibus` tables,
#'   `direction_diagnostics`, the underlying `directional_statistics`, and
#'   optionally `unit_metrics`.
#' @export
rc_report_condition_directions <- function(
    x,
    condition_col = NULL,
    celltype_col = NULL,
    conditions = NULL,
    cell_types = NULL,
    comparisons = NULL,
    reaction_ids = NULL,
    medium_scenarios = NULL,
    min_units = 5L,
    include_omnibus = TRUE,
    p_adjust_method = "BH",
    p_adjust_scope = c(
      "celltype_contrast_medium", "celltype_contrast", "celltype", "global"
    ),
    wilcox_correct = FALSE,
    eps = 1e-8,
    vmax_tolerance = 1e-6,
    direction_tolerance = 1e-10,
    include_unit_metrics = FALSE,
    source_label = NULL,
    outdir = NULL) {
  p_adjust_scope <- match.arg(p_adjust_scope)
  if (!is.numeric(direction_tolerance) || length(direction_tolerance) != 1L ||
      !is.finite(direction_tolerance) || direction_tolerance < 0) {
    stop(
      "`direction_tolerance` must be one non-negative finite number.",
      call. = FALSE
    )
  }
  if (!is.logical(include_unit_metrics) || length(include_unit_metrics) != 1L ||
      is.na(include_unit_metrics)) {
    stop(
      "`include_unit_metrics` must be one non-missing TRUE/FALSE value.",
      call. = FALSE
    )
  }
  source_label <- .rc_direction_report_source_label(source_label)

  statistics <- rc_test_condition_reactions(
    x = x,
    condition_col = condition_col,
    celltype_col = celltype_col,
    conditions = conditions,
    cell_types = cell_types,
    comparisons = comparisons,
    reaction_ids = reaction_ids,
    target_directions = c("forward", "reverse"),
    medium_scenarios = medium_scenarios,
    min_units = min_units,
    include_omnibus = include_omnibus,
    p_adjust_method = p_adjust_method,
    p_adjust_scope = p_adjust_scope,
    wilcox_correct = wilcox_correct,
    eps = eps,
    vmax_tolerance = vmax_tolerance,
    include_scores = TRUE,
    outdir = NULL
  )
  microcompass <- .rc_condition_stats_microcompass(x)
  condition_col <- statistics$params$condition_col
  celltype_col <- statistics$params$celltype_col
  unit_table <- .rc_direction_report_unit_table(
    statistics = statistics,
    microcompass = microcompass,
    condition_col = condition_col,
    celltype_col = celltype_col,
    direction_tolerance = direction_tolerance,
    source_label = source_label
  )
  diagnostics <- .rc_direction_report_diagnostics(
    unit_table,
    direction_tolerance = direction_tolerance
  )
  metric_table <- .rc_direction_report_metric_table(unit_table)
  analysis_unit <- statistics$params$analysis_unit
  biological_replicate_inference <- identical(
    analysis_unit,
    "sample_celltype"
  )
  reaction_pairwise <- .rc_direction_report_pairwise(
    metric_table = metric_table,
    comparisons = statistics$params$comparisons,
    min_units = statistics$params$min_units,
    wilcox_correct = statistics$params$wilcox_correct,
    p_adjust_method = statistics$params$p_adjust_method,
    p_adjust_scope = statistics$params$p_adjust_scope,
    analysis_unit = analysis_unit,
    biological_replicate_inference = biological_replicate_inference,
    direction_tolerance = direction_tolerance
  )
  reaction_omnibus <- data.frame()
  if (isTRUE(include_omnibus) && length(statistics$params$conditions) >= 3L) {
    reaction_omnibus <- .rc_direction_report_omnibus(
      metric_table = metric_table,
      conditions = statistics$params$conditions,
      min_units = statistics$params$min_units,
      p_adjust_method = statistics$params$p_adjust_method,
      p_adjust_scope = statistics$params$p_adjust_scope,
      analysis_unit = analysis_unit,
      biological_replicate_inference = biological_replicate_inference
    )
  }

  annotation <- .rc_direction_report_static_annotations(statistics)
  reaction_pairwise <- .rc_direction_report_attach_annotations(
    reaction_pairwise,
    annotation
  )
  reaction_omnibus <- .rc_direction_report_attach_annotations(
    reaction_omnibus,
    annotation
  )
  diagnostics <- .rc_direction_report_attach_annotations(
    diagnostics,
    annotation
  )

  answer <- list(
    directional_pairwise = statistics$pairwise,
    directional_omnibus = statistics$omnibus,
    reaction_pairwise = reaction_pairwise,
    reaction_omnibus = reaction_omnibus,
    direction_diagnostics = diagnostics,
    directional_statistics = statistics,
    params = c(
      statistics$params,
      list(
        direction_tolerance = direction_tolerance,
        source_label = source_label,
        any_direction_formula = paste(
          "max(forward_support, reverse_support);",
          "retain one available direction"
        ),
        directional_balance_formula =
          "forward_support - reverse_support"
      )
    ),
    reporting_policy = paste(
      "Direction-specific LP targets remain the primary COMPASS-style report.",
      "any_direction_support is a non-additive best-direction summary;",
      "directional_balance is support asymmetry and is not net flux."
    )
  )
  if (isTRUE(include_unit_metrics)) answer$unit_metrics <- unit_table
  class(answer) <- c("regcompass_condition_direction_report", "list")

  if (!is.null(outdir)) {
    valid_outdir <- is.character(outdir) && length(outdir) == 1L &&
      !is.na(outdir) && nzchar(trimws(outdir))
    if (!valid_outdir) {
      stop("`outdir` must be one non-empty path.", call. = FALSE)
    }
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    .rc_write_tsv_gz(
      answer$directional_pairwise,
      file.path(outdir, "condition_directional_pairwise.tsv.gz")
    )
    if (nrow(answer$directional_omnibus)) {
      .rc_write_tsv_gz(
        answer$directional_omnibus,
        file.path(outdir, "condition_directional_omnibus.tsv.gz")
      )
    }
    .rc_write_tsv_gz(
      answer$reaction_pairwise,
      file.path(outdir, "condition_reaction_level_pairwise.tsv.gz")
    )
    if (nrow(answer$reaction_omnibus)) {
      .rc_write_tsv_gz(
        answer$reaction_omnibus,
        file.path(outdir, "condition_reaction_level_omnibus.tsv.gz")
      )
    }
    .rc_write_tsv_gz(
      answer$direction_diagnostics,
      file.path(outdir, "condition_direction_diagnostics.tsv.gz")
    )
    if (isTRUE(include_unit_metrics)) {
      .rc_write_tsv_gz(
        answer$unit_metrics,
        file.path(outdir, "condition_direction_unit_metrics.tsv.gz")
      )
    }
    saveRDS(answer, file.path(outdir, "condition_direction_report.rds"))
  }
  answer
}
