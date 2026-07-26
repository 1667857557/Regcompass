# Direction-aware condition reporting for reversible reaction targets.

.rc_direction_source_label <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop("`source_label` must be NULL or one non-empty string.", call. = FALSE)
  }
  trimws(x)
}

.rc_direction_unit_table <- function(statistics, microcompass, source_label) {
  score <- as.matrix(statistics$score)
  rows <- rc_parse_microcompass_row_id(rownames(score))
  meta <- microcompass$unit_meta
  unit_id <- .rc_condition_stats_unit_ids(meta)
  meta$.unit_id <- unit_id
  meta <- meta[match(colnames(score), meta$.unit_id), , drop = FALSE]
  if (anyNA(meta$.unit_id)) {
    stop("Condition-report metadata do not align with score columns.", call. = FALSE)
  }
  condition <- as.character(meta[[statistics$params$condition_col]])
  cell_type <- as.character(meta[[statistics$params$celltype_col]])
  keep <- condition %in% statistics$params$conditions &
    cell_type %in% statistics$params$cell_types
  score <- score[, keep, drop = FALSE]
  condition <- condition[keep]
  cell_type <- cell_type[keep]

  groups <- split(
    seq_len(nrow(rows)),
    paste(rows$reaction_id, rows$medium_scenario, sep = "\001")
  )
  answer <- lapply(groups, function(index) {
    one <- rows[index, , drop = FALSE]
    forward <- index[one$target_direction == "forward"]
    reverse <- index[one$target_direction == "reverse"]
    if (length(forward) > 1L || length(reverse) > 1L) {
      stop(
        "At most one forward and one reverse target are allowed per reaction and medium.",
        call. = FALSE
      )
    }
    sf <- sr <- rep(NA_real_, ncol(score))
    if (length(forward)) sf <- score[forward, , drop = TRUE]
    if (length(reverse)) sr <- score[reverse, , drop = TRUE]
    ff <- is.finite(sf)
    fr <- is.finite(sr)
    paired <- ff & fr
    any_support <- rep(NA_real_, ncol(score))
    any_support[paired] <- pmax(sf[paired], sr[paired])
    any_support[ff & !fr] <- sf[ff & !fr]
    any_support[fr & !ff] <- sr[fr & !ff]
    balance <- rep(NA_real_, ncol(score))
    balance[paired] <- sf[paired] - sr[paired]
    data.frame(
      unit_id = colnames(score),
      reaction_id = one$reaction_id[[1L]],
      medium_scenario = one$medium_scenario[[1L]],
      condition = condition,
      cell_type = cell_type,
      forward_target_available = length(forward) == 1L,
      reverse_target_available = length(reverse) == 1L,
      score_forward = sf,
      score_reverse = sr,
      any_direction_support = any_support,
      directional_balance = balance,
      source_label = source_label,
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, answer)
  rownames(answer) <- NULL
  answer
}

.rc_direction_diagnostics <- function(unit_table, tolerance) {
  groups <- split(
    seq_len(nrow(unit_table)),
    interaction(
      unit_table$reaction_id, unit_table$medium_scenario,
      unit_table$cell_type, unit_table$condition,
      drop = TRUE, lex.order = TRUE
    )
  )
  answer <- lapply(groups, function(index) {
    one <- unit_table[index, , drop = FALSE]
    has_forward <- isTRUE(one$forward_target_available[[1L]])
    has_reverse <- isTRUE(one$reverse_target_available[[1L]])
    ff <- is.finite(one$score_forward)
    fr <- is.finite(one$score_reverse)
    paired <- ff & fr
    same_missing <- if (has_forward && has_reverse) identical(ff, fr) else NA
    max_difference <- if (any(paired)) {
      max(abs(one$score_forward[paired] - one$score_reverse[paired]))
    } else NA_real_
    indistinguishable <- has_forward && has_reverse && isTRUE(same_missing) &&
      any(paired) && is.finite(max_difference) && max_difference <= tolerance
    median_value <- function(x) {
      x <- x[is.finite(x)]
      if (length(x)) stats::median(x) else NA_real_
    }
    median_forward <- median_value(one$score_forward)
    median_reverse <- median_value(one$score_reverse)
    median_balance <- median_value(one$directional_balance)
    status <- if (has_forward && has_reverse) {
      if (!any(paired)) "bidirectional_no_paired_scores" else
        if (indistinguishable) "bidirectional_indistinguishable" else
          "bidirectional_distinguishable"
    } else if (has_forward) "forward_only" else
      if (has_reverse) "reverse_only" else "no_scored_direction"
    preferred <- if (identical(status, "forward_only")) "forward" else
      if (identical(status, "reverse_only")) "reverse" else
        if (identical(status, "bidirectional_indistinguishable")) {
          "indistinguishable"
        } else if (!is.finite(median_balance)) NA_character_ else
          if (median_balance > tolerance) "forward" else
            if (median_balance < -tolerance) "reverse" else "balanced"
    data.frame(
      reaction_id = one$reaction_id[[1L]],
      medium_scenario = one$medium_scenario[[1L]],
      cell_type = one$cell_type[[1L]],
      condition = one$condition[[1L]],
      n_units = nrow(one),
      n_forward_scores = sum(ff),
      n_reverse_scores = sum(fr),
      n_paired_direction_scores = sum(paired),
      median_forward_support = median_forward,
      median_reverse_support = median_reverse,
      median_any_direction_support = median_value(one$any_direction_support),
      median_directional_balance = median_balance,
      same_missing_pattern = same_missing,
      max_abs_forward_reverse_difference = max_difference,
      directionally_indistinguishable = indistinguishable,
      direction_pair_status = status,
      preferred_direction = preferred,
      source_label = one$source_label[[1L]],
      interpretation = paste(
        "Forward and reverse are independent counterfactual LP targets;",
        "directional balance is not net flux."
      ),
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, answer)
  rownames(answer) <- NULL
  answer
}

.rc_direction_metric_table <- function(unit_table) {
  n <- nrow(unit_table)
  any_direction <- data.frame(
    unit_id = unit_table$unit_id,
    reaction_id = unit_table$reaction_id,
    medium_scenario = unit_table$medium_scenario,
    condition = unit_table$condition,
    cell_type = unit_table$cell_type,
    report_metric = rep("any_direction_support", n),
    metric_value = unit_table$any_direction_support,
    source_label = unit_table$source_label,
    stringsAsFactors = FALSE
  )
  keep <- unit_table$forward_target_available & unit_table$reverse_target_available
  balance <- data.frame(
    unit_id = unit_table$unit_id[keep],
    reaction_id = unit_table$reaction_id[keep],
    medium_scenario = unit_table$medium_scenario[keep],
    condition = unit_table$condition[keep],
    cell_type = unit_table$cell_type[keep],
    report_metric = rep("directional_balance", sum(keep)),
    metric_value = unit_table$directional_balance[keep],
    source_label = unit_table$source_label[keep],
    stringsAsFactors = FALSE
  )
  rbind(any_direction, balance)
}

.rc_direction_adjust_groups <- function(scope, omnibus = FALSE) {
  if (isTRUE(omnibus)) {
    return(switch(
      scope,
      celltype_contrast_medium = c("report_metric", "cell_type", "medium_scenario"),
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

.rc_direction_metric_description <- function(metric) {
  if (identical(metric, "any_direction_support")) {
    list(
      definition = paste(
        "max(forward_support, reverse_support);",
        "one available direction is retained"
      ),
      semantics = "best-supported direction; non-additive reaction potential"
    )
  } else {
    list(
      definition = "forward_support - reverse_support",
      semantics = "directional support asymmetry; not net flux"
    )
  }
}

.rc_direction_pairwise <- function(metric_table, statistics, tolerance) {
  groups <- split(
    seq_len(nrow(metric_table)),
    interaction(
      metric_table$reaction_id, metric_table$medium_scenario,
      metric_table$cell_type, metric_table$report_metric,
      drop = TRUE, lex.order = TRUE
    )
  )
  rows <- list()
  k <- 0L
  for (index in groups) {
    one <- metric_table[index, , drop = FALSE]
    metric <- one$report_metric[[1L]]
    description <- .rc_direction_metric_description(metric)
    for (pair in statistics$params$comparisons) {
      a <- one$metric_value[one$condition == pair[[1L]]]
      b <- one$metric_value[one$condition == pair[[2L]]]
      effect <- .rc_condition_stats_effect(
        a, b, statistics$params$min_units, statistics$params$wilcox_correct
      )
      delta <- effect$delta_median_score_b_minus_a
      k <- k + 1L
      rows[[k]] <- data.frame(
        reaction_id = one$reaction_id[[1L]],
        medium_scenario = one$medium_scenario[[1L]],
        cell_type = one$cell_type[[1L]],
        report_metric = metric,
        condition_a = pair[[1L]],
        condition_b = pair[[2L]],
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
          if (delta > 0) pair[[2L]] else if (delta < 0) pair[[1L]] else "tie"
        } else NA_character_,
        direction_shift_b_minus_a = if (
          identical(metric, "directional_balance") && is.finite(delta)
        ) {
          if (delta > tolerance) "toward_forward" else
            if (delta < -tolerance) "toward_reverse" else "no_shift"
        } else NA_character_,
        metric_definition = description$definition,
        metric_semantics = description$semantics,
        analysis_unit = statistics$params$analysis_unit,
        inference_level = statistics$pairwise$inference_level[[1L]],
        descriptive_only = statistics$pairwise$descriptive_only[[1L]],
        biological_replicate_inference =
          statistics$pairwise$biological_replicate_inference[[1L]],
        source_label = one$source_label[[1L]],
        stringsAsFactors = FALSE
      )
    }
  }
  answer <- do.call(rbind, rows)
  answer <- .rc_condition_stats_adjust(
    answer,
    .rc_direction_adjust_groups(statistics$params$p_adjust_scope),
    statistics$params$p_adjust_method
  )
  answer$p_adjust_method <- statistics$params$p_adjust_method
  answer$p_adjust_scope <- statistics$params$p_adjust_scope
  answer <- answer[order(
    answer$cell_type, answer$condition_a, answer$condition_b,
    answer$medium_scenario, answer$report_metric, answer$p_adj,
    -abs(answer$rank_biserial_b_minus_a), answer$reaction_id,
    na.last = TRUE
  ), , drop = FALSE]
  rownames(answer) <- NULL
  answer
}

.rc_direction_omnibus <- function(metric_table, statistics) {
  groups <- split(
    seq_len(nrow(metric_table)),
    interaction(
      metric_table$reaction_id, metric_table$medium_scenario,
      metric_table$cell_type, metric_table$report_metric,
      drop = TRUE, lex.order = TRUE
    )
  )
  rows <- lapply(groups, function(index) {
    one <- metric_table[index, , drop = FALSE]
    metric <- one$report_metric[[1L]]
    description <- .rc_direction_metric_description(metric)
    values <- lapply(statistics$params$conditions, function(condition) {
      x <- one$metric_value[one$condition == condition]
      x[is.finite(x)]
    })
    names(values) <- statistics$params$conditions
    counts <- vapply(values, length, integer(1))
    p_value <- NA_real_
    status <- "insufficient_units"
    if (all(counts >= statistics$params$min_units)) {
      x <- unlist(values, use.names = FALSE)
      g <- factor(rep(names(values), lengths(values)), levels = names(values))
      if (length(unique(x)) == 1L) {
        p_value <- 1
        status <- "constant_equal"
      } else {
        test <- suppressWarnings(stats::kruskal.test(x = x, g = g))
        p_value <- unname(test$p.value)
        status <- if (is.finite(p_value)) "ok" else "test_failed"
      }
    }
    data.frame(
      reaction_id = one$reaction_id[[1L]],
      medium_scenario = one$medium_scenario[[1L]],
      cell_type = one$cell_type[[1L]],
      report_metric = metric,
      n_conditions = length(values),
      n_units_total = sum(counts),
      min_units_per_condition = min(counts),
      max_units_per_condition = max(counts),
      units_by_condition = paste(paste(names(counts), counts, sep = "="), collapse = ";"),
      p_value = p_value,
      test_status = status,
      metric_definition = description$definition,
      metric_semantics = description$semantics,
      analysis_unit = statistics$params$analysis_unit,
      inference_level = statistics$pairwise$inference_level[[1L]],
      descriptive_only = statistics$pairwise$descriptive_only[[1L]],
      biological_replicate_inference =
        statistics$pairwise$biological_replicate_inference[[1L]],
      source_label = one$source_label[[1L]],
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  answer <- .rc_condition_stats_adjust(
    answer,
    .rc_direction_adjust_groups(statistics$params$p_adjust_scope, omnibus = TRUE),
    statistics$params$p_adjust_method
  )
  answer$p_adjust_method <- statistics$params$p_adjust_method
  answer$p_adjust_scope <- statistics$params$p_adjust_scope
  answer <- answer[order(
    answer$cell_type, answer$medium_scenario, answer$report_metric,
    answer$p_adj, answer$reaction_id, na.last = TRUE
  ), , drop = FALSE]
  rownames(answer) <- NULL
  answer
}

.rc_direction_static_annotations <- function(statistics) {
  source <- if (nrow(statistics$pairwise)) statistics$pairwise else statistics$omnibus
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
  if (!"reaction_id" %in% columns) return(NULL)
  answer <- unique(source[, columns, drop = FALSE])
  answer[!duplicated(answer$reaction_id), , drop = FALSE]
}

.rc_direction_attach_annotations <- function(data, annotation) {
  if (is.null(annotation) || !nrow(data)) return(data)
  index <- match(data$reaction_id, annotation$reaction_id)
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
    stop("`direction_tolerance` must be one non-negative finite number.", call. = FALSE)
  }
  if (!is.logical(include_unit_metrics) || length(include_unit_metrics) != 1L ||
      is.na(include_unit_metrics)) {
    stop("`include_unit_metrics` must be one non-missing TRUE/FALSE value.", call. = FALSE)
  }
  source_label <- .rc_direction_source_label(source_label)
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
  statistics$pairwise$source_label <- source_label
  if (nrow(statistics$omnibus)) statistics$omnibus$source_label <- source_label
  microcompass <- .rc_condition_stats_microcompass(x)
  unit_table <- .rc_direction_unit_table(statistics, microcompass, source_label)
  diagnostics <- .rc_direction_diagnostics(unit_table, direction_tolerance)
  metric_table <- .rc_direction_metric_table(unit_table)
  reaction_pairwise <- .rc_direction_pairwise(
    metric_table, statistics, direction_tolerance
  )
  reaction_omnibus <- data.frame()
  if (isTRUE(include_omnibus) && length(statistics$params$conditions) >= 3L) {
    reaction_omnibus <- .rc_direction_omnibus(metric_table, statistics)
  }
  annotation <- .rc_direction_static_annotations(statistics)
  reaction_pairwise <- .rc_direction_attach_annotations(reaction_pairwise, annotation)
  reaction_omnibus <- .rc_direction_attach_annotations(reaction_omnibus, annotation)
  diagnostics <- .rc_direction_attach_annotations(diagnostics, annotation)

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
        directional_balance_formula = "forward_support - reverse_support"
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
    if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) ||
        !nzchar(trimws(outdir))) {
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
