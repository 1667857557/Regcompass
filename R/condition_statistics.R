# Condition-associated reaction statistics on shared microCOMPASS scores.

.rc_condition_stats_microcompass <- function(x) {
  if (is.list(x) && !is.null(x$microcompass)) x <- x$microcompass
  if (!is.list(x) || is.null(x$penalty) || is.null(x$vmax) ||
      is.null(x$feasible) || is.null(x$unit_meta)) {
    stop(
      "`x` must be a microCOMPASS result or a RegCompass result containing `microcompass`.",
      call. = FALSE
    )
  }
  x
}

.rc_condition_stats_column <- function(meta, supplied, candidates, label) {
  if (!is.null(supplied)) {
    if (!is.character(supplied) || length(supplied) != 1L ||
        is.na(supplied) || !nzchar(supplied) || !supplied %in% colnames(meta)) {
      stop("`", label, "` must name one column in `unit_meta`.", call. = FALSE)
    }
    return(supplied)
  }
  hit <- candidates[candidates %in% colnames(meta)]
  if (!length(hit)) {
    stop(
      "Could not infer `", label, "`; supply it explicitly. Available columns: ",
      paste(colnames(meta), collapse = ", "),
      call. = FALSE
    )
  }
  hit[[1L]]
}

.rc_condition_stats_unit_ids <- function(meta) {
  unit_id <- if ("unit_id" %in% colnames(meta)) {
    as.character(meta$unit_id)
  } else if ("pool_id" %in% colnames(meta)) {
    as.character(meta$pool_id)
  } else if (!is.null(rownames(meta)) && !anyDuplicated(rownames(meta)) &&
             all(nzchar(rownames(meta)))) {
    rownames(meta)
  } else {
    stop("`unit_meta` lacks unique unit_id/pool_id values.", call. = FALSE)
  }
  if (anyNA(unit_id) || any(!nzchar(trimws(unit_id))) || anyDuplicated(unit_id)) {
    stop("microCOMPASS unit IDs must be unique and non-empty.", call. = FALSE)
  }
  unit_id
}

.rc_condition_stats_comparisons <- function(comparisons, conditions) {
  if (is.null(comparisons)) {
    return(utils::combn(conditions, 2L, simplify = FALSE))
  }
  if (is.data.frame(comparisons) || is.matrix(comparisons)) {
    if (ncol(comparisons) != 2L) {
      stop("`comparisons` must have exactly two columns.", call. = FALSE)
    }
    comparisons <- lapply(seq_len(nrow(comparisons)), function(i) {
      as.character(comparisons[i, , drop = TRUE])
    })
  }
  if (!is.list(comparisons) || !length(comparisons) ||
      any(!vapply(comparisons, function(z) {
        is.character(z) && length(z) == 2L && !anyNA(z) &&
          all(nzchar(trimws(z))) && z[[1L]] != z[[2L]]
      }, logical(1)))) {
    stop(
      "`comparisons` must be NULL, a two-column object, or a list of two-condition character vectors.",
      call. = FALSE
    )
  }
  comparisons <- lapply(comparisons, as.character)
  unknown <- setdiff(unique(unlist(comparisons, use.names = FALSE)), conditions)
  if (length(unknown)) {
    stop("Unknown conditions in `comparisons`: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  keys <- vapply(comparisons, paste, collapse = "\001", character(1))
  comparisons[!duplicated(keys)]
}

.rc_condition_stats_adjust <- function(data, group_cols, method) {
  data$p_adj <- NA_real_
  if (!nrow(data)) return(data)
  groups <- if (!length(group_cols)) {
    list(all = seq_len(nrow(data)))
  } else {
    split(
      seq_len(nrow(data)),
      interaction(data[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
    )
  }
  for (rows in groups) {
    valid <- rows[is.finite(data$p_value[rows])]
    if (length(valid)) {
      data$p_adj[valid] <- stats::p.adjust(data$p_value[valid], method = method)
    }
  }
  data
}

.rc_condition_stats_effect <- function(a, b, min_units, wilcox_correct) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  n_a <- length(a)
  n_b <- length(b)
  median_a <- if (n_a) stats::median(a) else NA_real_
  median_b <- if (n_b) stats::median(b) else NA_real_
  mean_a <- if (n_a) mean(a) else NA_real_
  mean_b <- if (n_b) mean(b) else NA_real_
  base <- list(
    n_a = n_a,
    n_b = n_b,
    median_score_a = median_a,
    median_score_b = median_b,
    delta_median_score_b_minus_a = median_b - median_a,
    mean_score_a = mean_a,
    mean_score_b = mean_b,
    cohens_d_b_minus_a = NA_real_,
    rank_biserial_b_minus_a = NA_real_,
    common_language_b_greater_a = NA_real_,
    p_value = NA_real_,
    test_status = "insufficient_units"
  )
  if (n_a < min_units || n_b < min_units) return(base)

  combined <- c(b, a)
  ranks <- rank(combined, ties.method = "average")
  u_b <- sum(ranks[seq_len(n_b)]) - n_b * (n_b + 1) / 2
  probability <- u_b / (n_a * n_b)
  base$common_language_b_greater_a <- probability
  base$rank_biserial_b_minus_a <- 2 * probability - 1

  pooled_variance <- (
    (n_a - 1) * stats::var(a) + (n_b - 1) * stats::var(b)
  ) / (n_a + n_b - 2)
  if (is.finite(pooled_variance) && pooled_variance > 0) {
    base$cohens_d_b_minus_a <- (mean_b - mean_a) / sqrt(pooled_variance)
  } else if (isTRUE(all.equal(mean_a, mean_b))) {
    base$cohens_d_b_minus_a <- 0
  }

  if (length(unique(combined)) == 1L) {
    base$p_value <- 1
    base$test_status <- "constant_equal"
    return(base)
  }
  test <- suppressWarnings(stats::wilcox.test(
    b, a,
    alternative = "two.sided",
    exact = FALSE,
    correct = wilcox_correct
  ))
  base$p_value <- unname(test$p.value)
  base$test_status <- if (is.finite(base$p_value)) "ok" else "test_failed"
  base
}

#' Test reaction support differences between conditions
#'
#' Compares the same reaction, direction, and medium between conditions within
#' each selected cell type. All units must have been scored against the shared
#' structural GEM stored in the microCOMPASS result. The tested score is
#' `-log(penalty / (omega * vmax) + eps)`, so larger values indicate stronger
#' multiome support for the target reaction direction.
#'
#' Pairwise differences use two-sided Wilcoxon rank-sum tests and report median
#' shifts, rank-biserial correlation, common-language effect size, and Cohen's
#' d. When at least three conditions are selected, an optional Kruskal-Wallis
#' omnibus test is also performed. P values are adjusted across reaction targets
#' using the requested scope.
#'
#' Metacell tests quantify within-dataset condition-associated separation. They
#' do not turn metacells into independent biological replicates. The returned
#' inference fields record this distinction explicitly.
#'
#' @param x A microCOMPASS result, such as `step5`, or a complete RegCompass
#'   result containing a `microcompass` element.
#' @param condition_col Metadata column defining conditions. When `NULL`, common
#'   condition-column names are searched in order.
#' @param celltype_col Metadata column defining cell types. When `NULL`, common
#'   cell-type-column names are searched in order.
#' @param conditions Optional conditions to retain. The default uses all.
#' @param cell_types Optional cell types to retain. The default uses all.
#' @param comparisons Optional pairwise contrasts. Supply a list of character
#'   vectors of length two or a two-column matrix/data frame. The default tests
#'   every pair of retained conditions.
#' @param reaction_ids,target_directions,medium_scenarios Optional target filters.
#' @param min_units Minimum finite units required in every tested condition.
#' @param include_omnibus Run Kruskal-Wallis tests when at least three conditions
#'   are retained.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param p_adjust_scope Multiple-testing scope: `"celltype_contrast_medium"`
#'   (default), `"celltype_contrast"`, `"celltype"`, or `"global"`.
#' @param wilcox_correct Apply continuity correction in Wilcoxon tests.
#' @param eps Positive offset used in the support-score transformation.
#' @param vmax_tolerance Relative tolerance for verifying that target vmax is
#'   invariant across units under the shared GEM.
#' @param include_scores Include the filtered unit-level score matrix in the
#'   returned object.
#' @param outdir Optional directory for pairwise/omnibus TSV files and an RDS.
#' @return A `regcompass_condition_statistics` list with `pairwise`, `omnibus`,
#'   `params`, `inference_policy`, and optionally `score`.
#' @export
rc_test_condition_reactions <- function(
    x,
    condition_col = NULL,
    celltype_col = NULL,
    conditions = NULL,
    cell_types = NULL,
    comparisons = NULL,
    reaction_ids = NULL,
    target_directions = NULL,
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
    include_scores = FALSE,
    outdir = NULL) {
  p_adjust_scope <- match.arg(p_adjust_scope)
  if (!is.numeric(min_units) || length(min_units) != 1L ||
      !is.finite(min_units) || min_units < 2 ||
      abs(min_units - round(min_units)) > sqrt(.Machine$double.eps)) {
    stop("`min_units` must be one integer of at least two.", call. = FALSE)
  }
  min_units <- as.integer(min_units)
  if (!is.logical(include_omnibus) || length(include_omnibus) != 1L ||
      is.na(include_omnibus) || !is.logical(wilcox_correct) ||
      length(wilcox_correct) != 1L || is.na(wilcox_correct) ||
      !is.logical(include_scores) || length(include_scores) != 1L ||
      is.na(include_scores)) {
    stop("Logical controls must be one non-missing TRUE/FALSE value.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0 ||
      !is.numeric(vmax_tolerance) || length(vmax_tolerance) != 1L ||
      !is.finite(vmax_tolerance) || vmax_tolerance < 0) {
    stop("`eps` must be positive and `vmax_tolerance` non-negative.", call. = FALSE)
  }
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) || !p_adjust_method %in% stats::p.adjust.methods) {
    stop("`p_adjust_method` is not supported by `stats::p.adjust()`.", call. = FALSE)
  }

  microcompass <- .rc_condition_stats_microcompass(x)
  penalty <- as.matrix(microcompass$penalty)
  vmax <- as.matrix(microcompass$vmax)
  feasible <- as.matrix(microcompass$feasible)
  valid_matrix <- function(z) {
    is.numeric(z) && !is.null(rownames(z)) && !is.null(colnames(z)) &&
      !anyDuplicated(rownames(z)) && !anyDuplicated(colnames(z))
  }
  if (!valid_matrix(penalty) || !valid_matrix(vmax) ||
      !is.logical(feasible) || is.null(dim(feasible)) ||
      !identical(dimnames(penalty), dimnames(vmax)) ||
      !identical(dimnames(penalty), dimnames(feasible))) {
    stop("microCOMPASS penalty, vmax, and feasible matrices must align exactly.",
         call. = FALSE)
  }

  meta <- microcompass$unit_meta
  if (!is.data.frame(meta)) stop("microCOMPASS `unit_meta` must be a data frame.",
                                 call. = FALSE)
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
  unit_id <- .rc_condition_stats_unit_ids(meta)
  if (!setequal(colnames(penalty), unit_id)) {
    stop("microCOMPASS matrices and unit metadata contain different units.",
         call. = FALSE)
  }
  meta$.unit_id <- unit_id
  meta <- meta[match(colnames(penalty), meta$.unit_id), , drop = FALSE]
  condition_values <- as.character(meta[[condition_col]])
  celltype_values <- as.character(meta[[celltype_col]])
  if (anyNA(condition_values) || anyNA(celltype_values) ||
      any(!nzchar(trimws(condition_values))) ||
      any(!nzchar(trimws(celltype_values)))) {
    stop("Condition and cell-type metadata must be complete.", call. = FALSE)
  }

  available_conditions <- unique(condition_values)
  available_cell_types <- unique(celltype_values)
  conditions <- if (is.null(conditions)) available_conditions else as.character(conditions)
  cell_types <- if (is.null(cell_types)) available_cell_types else as.character(cell_types)
  unknown_conditions <- setdiff(conditions, available_conditions)
  unknown_cell_types <- setdiff(cell_types, available_cell_types)
  if (length(unknown_conditions) || length(unknown_cell_types)) {
    stop(
      paste0(
        if (length(unknown_conditions)) {
          paste0("Unknown conditions: ", paste(unknown_conditions, collapse = ", "), ". ")
        } else "",
        if (length(unknown_cell_types)) {
          paste0("Unknown cell types: ", paste(unknown_cell_types, collapse = ", "), ".")
        } else ""
      ),
      call. = FALSE
    )
  }
  conditions <- unique(conditions)
  cell_types <- unique(cell_types)
  if (length(conditions) < 2L) {
    stop("At least two conditions are required.", call. = FALSE)
  }
  comparisons <- .rc_condition_stats_comparisons(comparisons, conditions)

  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  row_meta$row_id <- rownames(penalty)
  row_keep <- rep(TRUE, nrow(row_meta))
  if (!is.null(reaction_ids)) {
    row_keep <- row_keep & row_meta$reaction_id %in% as.character(reaction_ids)
  }
  if (!is.null(target_directions)) {
    row_keep <- row_keep &
      row_meta$target_direction %in% as.character(target_directions)
  }
  if (!is.null(medium_scenarios)) {
    row_keep <- row_keep &
      row_meta$medium_scenario %in% as.character(medium_scenarios)
  }
  if (!any(row_keep)) stop("No reaction targets remain after filtering.", call. = FALSE)
  penalty <- penalty[row_keep, , drop = FALSE]
  vmax <- vmax[row_keep, , drop = FALSE]
  feasible <- feasible[row_keep, , drop = FALSE]
  row_meta <- row_meta[row_keep, , drop = FALSE]

  vmax_invariant <- vapply(seq_len(nrow(vmax)), function(i) {
    values <- vmax[i, is.finite(vmax[i, ]), drop = TRUE]
    if (length(values) <= 1L) return(TRUE)
    diff(range(values)) <= vmax_tolerance * max(1, abs(stats::median(values)))
  }, logical(1))
  if (any(!vmax_invariant)) {
    stop(
      "Target vmax differs across units despite the shared structural GEM: ",
      paste(utils::head(rownames(vmax)[!vmax_invariant], 10L), collapse = ", "),
      call. = FALSE
    )
  }

  omega <- microcompass$params$omega %||% 0.95
  if (!is.numeric(omega) || length(omega) != 1L || !is.finite(omega) ||
      omega <= 0 || omega > 1) {
    stop("microCOMPASS `omega` must be in (0, 1].", call. = FALSE)
  }
  required_flux <- omega * vmax
  normalized <- penalty / required_flux
  score <- -log(normalized + eps)
  score[!feasible | !is.finite(normalized) | normalized < 0 |
          !is.finite(score) | required_flux <= 0] <- NA_real_

  selected_units <- condition_values %in% conditions & celltype_values %in% cell_types
  score <- score[, selected_units, drop = FALSE]
  meta <- meta[selected_units, , drop = FALSE]
  condition_values <- as.character(meta[[condition_col]])
  celltype_values <- as.character(meta[[celltype_col]])

  analysis_unit <- as.character(microcompass$params$unit %||% "unknown")
  biological_replicate_inference <- identical(analysis_unit, "sample_celltype")
  inference_level <- if (biological_replicate_inference) {
    "biological_sample_celltype"
  } else if (identical(analysis_unit, "metacell")) {
    "metacell_within_dataset"
  } else {
    paste0(analysis_unit, "_within_dataset")
  }

  pairwise_rows <- list()
  pairwise_index <- 0L
  for (cell_type in cell_types) {
    cell_keep <- celltype_values == cell_type
    for (pair in comparisons) {
      condition_a <- pair[[1L]]
      condition_b <- pair[[2L]]
      a_keep <- cell_keep & condition_values == condition_a
      b_keep <- cell_keep & condition_values == condition_b
      for (i in seq_len(nrow(score))) {
        effect <- .rc_condition_stats_effect(
          score[i, a_keep], score[i, b_keep], min_units, wilcox_correct
        )
        pairwise_index <- pairwise_index + 1L
        pairwise_rows[[pairwise_index]] <- data.frame(
          row_id = row_meta$row_id[[i]],
          reaction_id = row_meta$reaction_id[[i]],
          target_direction = row_meta$target_direction[[i]],
          medium_scenario = row_meta$medium_scenario[[i]],
          cell_type = cell_type,
          condition_a = condition_a,
          condition_b = condition_b,
          n_a = effect$n_a,
          n_b = effect$n_b,
          median_score_a = effect$median_score_a,
          median_score_b = effect$median_score_b,
          delta_median_score_b_minus_a =
            effect$delta_median_score_b_minus_a,
          mean_score_a = effect$mean_score_a,
          mean_score_b = effect$mean_score_b,
          cohens_d_b_minus_a = effect$cohens_d_b_minus_a,
          rank_biserial_b_minus_a = effect$rank_biserial_b_minus_a,
          common_language_b_greater_a = effect$common_language_b_greater_a,
          p_value = effect$p_value,
          test_status = effect$test_status,
          higher_supported_condition = if (
            is.finite(effect$delta_median_score_b_minus_a)
          ) {
            if (effect$delta_median_score_b_minus_a > 0) condition_b else
              if (effect$delta_median_score_b_minus_a < 0) condition_a else "tie"
          } else {
            NA_character_
          },
          analysis_unit = analysis_unit,
          inference_level = inference_level,
          descriptive_only = !biological_replicate_inference,
          biological_replicate_inference = biological_replicate_inference,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  pairwise <- do.call(rbind, pairwise_rows)
  pairwise_groups <- switch(
    p_adjust_scope,
    celltype_contrast_medium = c(
      "cell_type", "condition_a", "condition_b", "medium_scenario"
    ),
    celltype_contrast = c("cell_type", "condition_a", "condition_b"),
    celltype = "cell_type",
    global = character()
  )
  pairwise <- .rc_condition_stats_adjust(
    pairwise, pairwise_groups, p_adjust_method
  )
  pairwise$p_adjust_method <- p_adjust_method
  pairwise$p_adjust_scope <- p_adjust_scope
  pairwise <- pairwise[order(
    pairwise$cell_type, pairwise$condition_a, pairwise$condition_b,
    pairwise$medium_scenario, pairwise$p_adj,
    -abs(pairwise$rank_biserial_b_minus_a),
    pairwise$reaction_id, pairwise$target_direction,
    na.last = TRUE
  ), , drop = FALSE]
  rownames(pairwise) <- NULL

  omnibus <- data.frame()
  if (isTRUE(include_omnibus) && length(conditions) >= 3L) {
    omnibus_rows <- list()
    omnibus_index <- 0L
    for (cell_type in cell_types) {
      cell_keep <- celltype_values == cell_type
      for (i in seq_len(nrow(score))) {
        groups <- lapply(conditions, function(condition) {
          values <- score[i, cell_keep & condition_values == condition]
          values[is.finite(values)]
        })
        names(groups) <- conditions
        counts <- vapply(groups, length, integer(1))
        valid <- all(counts >= min_units)
        p_value <- NA_real_
        status <- "insufficient_units"
        if (valid) {
          combined <- unlist(groups, use.names = FALSE)
          if (length(unique(combined)) == 1L) {
            p_value <- 1
            status <- "constant_equal"
          } else {
            test <- suppressWarnings(stats::kruskal.test(groups))
            p_value <- unname(test$p.value)
            status <- if (is.finite(p_value)) "ok" else "test_failed"
          }
        }
        omnibus_index <- omnibus_index + 1L
        omnibus_rows[[omnibus_index]] <- data.frame(
          row_id = row_meta$row_id[[i]],
          reaction_id = row_meta$reaction_id[[i]],
          target_direction = row_meta$target_direction[[i]],
          medium_scenario = row_meta$medium_scenario[[i]],
          cell_type = cell_type,
          n_conditions = length(conditions),
          n_units_total = sum(counts),
          min_units_per_condition = min(counts),
          max_units_per_condition = max(counts),
          units_by_condition = paste(
            paste(names(counts), counts, sep = "="), collapse = ";"
          ),
          p_value = p_value,
          test_status = status,
          analysis_unit = analysis_unit,
          inference_level = inference_level,
          descriptive_only = !biological_replicate_inference,
          biological_replicate_inference = biological_replicate_inference,
          stringsAsFactors = FALSE
        )
      }
    }
    omnibus <- do.call(rbind, omnibus_rows)
    omnibus_groups <- switch(
      p_adjust_scope,
      celltype_contrast_medium = c("cell_type", "medium_scenario"),
      celltype_contrast = "cell_type",
      celltype = "cell_type",
      global = character()
    )
    omnibus <- .rc_condition_stats_adjust(
      omnibus, omnibus_groups, p_adjust_method
    )
    omnibus$p_adjust_method <- p_adjust_method
    omnibus$p_adjust_scope <- p_adjust_scope
    omnibus <- omnibus[order(
      omnibus$cell_type, omnibus$medium_scenario, omnibus$p_adj,
      omnibus$reaction_id, omnibus$target_direction,
      na.last = TRUE
    ), , drop = FALSE]
    rownames(omnibus) <- NULL
  }

  inference_policy <- if (biological_replicate_inference) {
    paste(
      "Scores were compared across sample-by-cell-type units under one shared",
      "structural GEM. P values require independent biological samples within",
      "each condition."
    )
  } else {
    paste(
      "Scores were compared across metacells under one shared structural GEM.",
      "P values quantify within-dataset condition-associated metacell",
      "separation and are not biological-replicate-level treatment inference."
    )
  }
  answer <- list(
    pairwise = pairwise,
    omnibus = omnibus,
    params = list(
      condition_col = condition_col,
      celltype_col = celltype_col,
      conditions = conditions,
      cell_types = cell_types,
      comparisons = comparisons,
      min_units = min_units,
      score_formula = "-log(penalty / (omega * vmax) + eps)",
      omega = omega,
      eps = eps,
      p_adjust_method = p_adjust_method,
      p_adjust_scope = p_adjust_scope,
      wilcox_correct = wilcox_correct,
      analysis_unit = analysis_unit
    ),
    inference_policy = inference_policy
  )
  if (isTRUE(include_scores)) answer$score <- score
  class(answer) <- c("regcompass_condition_statistics", "list")

  if (!is.null(outdir)) {
    if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) ||
        !nzchar(outdir)) {
      stop("`outdir` must be one non-empty path.", call. = FALSE)
    }
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    .rc_write_tsv_gz(
      pairwise,
      file.path(outdir, "condition_reaction_pairwise.tsv.gz")
    )
    if (nrow(omnibus)) {
      .rc_write_tsv_gz(
        omnibus,
        file.path(outdir, "condition_reaction_omnibus.tsv.gz")
      )
    }
    saveRDS(answer, file.path(outdir, "condition_reaction_statistics.rds"))
  }
  answer
}
