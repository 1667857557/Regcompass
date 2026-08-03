.rc_condition_penalty_comparison <- function(
    microcompass, condition_col = "condition", celltype_col = "cell_type",
    eps = 1e-8, vmax_tolerance = 1e-6) {
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(vmax_tolerance) || length(vmax_tolerance) != 1L ||
      !is.finite(vmax_tolerance) || vmax_tolerance < 0) {
    stop("`vmax_tolerance` must be one finite non-negative number.",
         call. = FALSE)
  }
  penalty <- as.matrix(microcompass$penalty)
  vmax <- as.matrix(microcompass$vmax)
  meta <- microcompass$unit_meta
  required <- c(condition_col, celltype_col)
  valid_matrix <- function(x) {
    is.numeric(x) && !is.null(rownames(x)) && !is.null(colnames(x)) &&
      !anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x))
  }
  if (!valid_matrix(penalty) || !valid_matrix(vmax) ||
      !setequal(rownames(penalty), rownames(vmax)) ||
      !setequal(colnames(penalty), colnames(vmax))) {
    stop("microCOMPASS penalty and vmax matrices must align.", call. = FALSE)
  }
  vmax <- vmax[rownames(penalty), colnames(penalty), drop = FALSE]
  if (!is.data.frame(meta) || !all(required %in% colnames(meta))) {
    stop("microCOMPASS unit metadata lack condition/cell-type columns.",
         call. = FALSE)
  }
  unit_id <- if ("unit_id" %in% colnames(meta)) {
    as.character(meta$unit_id)
  } else if ("pool_id" %in% colnames(meta)) {
    as.character(meta$pool_id)
  } else {
    stop("microCOMPASS unit metadata lack unit_id/pool_id.", call. = FALSE)
  }
  if (anyNA(unit_id) || any(!nzchar(trimws(unit_id))) || anyDuplicated(unit_id) ||
      !setequal(colnames(penalty), unit_id)) {
    stop("microCOMPASS units and metadata are invalid or different.",
         call. = FALSE)
  }
  meta$unit_id <- unit_id
  meta <- meta[match(colnames(penalty), meta$unit_id), , drop = FALSE]
  condition_value <- trimws(as.character(meta[[condition_col]]))
  celltype_value <- trimws(as.character(meta[[celltype_col]]))
  if (anyNA(condition_value) || any(!nzchar(condition_value)) ||
      anyNA(celltype_value) || any(!nzchar(celltype_value))) {
    stop("microCOMPASS condition/cell-type metadata are incomplete.",
         call. = FALSE)
  }

  omega <- microcompass$params$omega %||% 0.95
  if (!is.numeric(omega) || length(omega) != 1L ||
      !is.finite(omega) || omega <= 0 || omega > 1) {
    stop("microCOMPASS `omega` must be in (0, 1].", call. = FALSE)
  }
  vmax_invariant <- vapply(seq_len(nrow(vmax)), function(i) {
    values <- vmax[i, is.finite(vmax[i, ]), drop = TRUE]
    if (length(values) <= 1L) return(TRUE)
    diff(range(values)) <= vmax_tolerance *
      max(1, abs(stats::median(values)))
  }, logical(1))
  if (any(!vmax_invariant)) {
    stop(
      "Target vmax differs among units assigned to the same structural row: ",
      paste(utils::head(rownames(vmax)[!vmax_invariant], 10L), collapse = ", "),
      call. = FALSE
    )
  }
  required_target_flux <- omega * vmax
  normalized <- matrix(
    NA_real_, nrow(penalty), ncol(penalty), dimnames = dimnames(penalty)
  )
  valid <- is.finite(penalty) & is.finite(required_target_flux) &
    required_target_flux > 0
  normalized[valid] <- penalty[valid] / required_target_flux[valid]

  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  row_meta$row_id <- rownames(penalty)
  strata <- unique(meta[, c(condition_col, celltype_col), drop = FALSE])
  summary_rows <- lapply(seq_len(nrow(strata)), function(i) {
    condition <- as.character(strata[[condition_col]][[i]])
    cell_type <- as.character(strata[[celltype_col]][[i]])
    column_keep <- condition_value == condition & celltype_value == cell_type
    row_keep <- is.na(row_meta$cell_type) | row_meta$cell_type == cell_type
    if (!any(row_keep)) return(NULL)
    row_index <- which(row_keep)
    median_matrix <- function(x) {
      value <- matrixStats::rowMedians(
        x[row_index, column_keep, drop = FALSE], na.rm = TRUE
      )
      value[is.nan(value)] <- NA_real_
      value
    }
    median_penalty <- median_matrix(penalty)
    median_vmax <- median_matrix(vmax)
    median_required_flux <- median_matrix(required_target_flux)
    median_normalized <- median_matrix(normalized)
    data.frame(
      row_id = row_meta$row_id[row_index],
      reaction_id = row_meta$reaction_id[row_index],
      target_direction = row_meta$target_direction[row_index],
      medium_scenario = row_meta$medium_scenario[row_index],
      condition = condition,
      cell_type = cell_type,
      median_penalty = median_penalty,
      median_vmax = median_vmax,
      median_required_target_flux = median_required_flux,
      median_penalty_per_target_flux = median_normalized,
      support_score = -log(median_normalized + eps),
      priority_rank = NA_integer_,
      ranking_metric = "minimum_penalty_per_required_target_flux",
      ranking_scope = "condition_x_celltype_x_medium",
      n_metacells = sum(column_keep),
      stringsAsFactors = FALSE
    )
  })
  summary_rows <- Filter(Negate(is.null), summary_rows)
  if (!length(summary_rows)) {
    stop("No condition rows match their cell-type structural models.",
         call. = FALSE)
  }
  ranking <- do.call(rbind, summary_rows)
  rank_group <- interaction(
    ranking$cell_type, ranking$condition, ranking$medium_scenario,
    drop = TRUE, lex.order = TRUE
  )
  for (rows in split(seq_len(nrow(ranking)), rank_group)) {
    ranking$priority_rank[rows] <- as.integer(rank(
      ranking$median_penalty_per_target_flux[rows],
      ties.method = "min", na.last = "keep"
    ))
  }
  ranking <- ranking[order(
    ranking$cell_type, ranking$condition, ranking$medium_scenario,
    ranking$priority_rank, ranking$reaction_id,
    ranking$target_direction, na.last = TRUE
  ), , drop = FALSE]
  rownames(ranking) <- NULL

  contrast_rows <- list()
  index <- 0L
  for (cell_type in unique(ranking$cell_type)) {
    one <- ranking[ranking$cell_type == cell_type, , drop = FALSE]
    conditions <- unique(as.character(one$condition))
    if (length(conditions) < 2L) next
    for (pair in utils::combn(conditions, 2L, simplify = FALSE)) {
      a <- one[one$condition == pair[[1L]], , drop = FALSE]
      b <- one[one$condition == pair[[2L]], , drop = FALSE]
      b <- b[match(a$row_id, b$row_id), , drop = FALSE]
      if (anyNA(b$row_id)) {
        stop("Conditions within a cell type contain different reaction rows.",
             call. = FALSE)
      }
      index <- index + 1L
      contrast_rows[[index]] <- data.frame(
        row_id = a$row_id,
        reaction_id = a$reaction_id,
        target_direction = a$target_direction,
        medium_scenario = a$medium_scenario,
        cell_type = cell_type,
        condition_a = pair[[1L]],
        condition_b = pair[[2L]],
        median_penalty_a = a$median_penalty,
        median_penalty_b = b$median_penalty,
        median_penalty_per_target_flux_a =
          a$median_penalty_per_target_flux,
        median_penalty_per_target_flux_b =
          b$median_penalty_per_target_flux,
        priority_rank_a = a$priority_rank,
        priority_rank_b = b$priority_rank,
        delta_support_b_minus_a = b$support_score - a$support_score,
        higher_supported_condition = ifelse(
          b$support_score > a$support_score, pair[[2L]],
          ifelse(a$support_score > b$support_score, pair[[1L]], "tie")
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  contrast <- if (length(contrast_rows)) do.call(rbind, contrast_rows) else
    data.frame()
  if (nrow(contrast)) rownames(contrast) <- NULL
  list(
    summary = ranking,
    ranking = ranking,
    contrast = contrast,
    analysis_mode = if (length(unique(ranking$condition)) == 1L) {
      "single_condition_reaction_ranking"
    } else {
      "multi_condition_reaction_ranking_and_pairwise_comparison"
    },
    ranking_formula = "penalty / (omega * vmax)",
    ranking_scope = "condition_x_celltype_x_medium",
    structural_scope = "cell_type_x_medium",
    comparison_workflow = paste(
      "conditions are compared only within one cell type on that cell type's",
      "medium-specific union GEM"
    )
  )
}

.rc_condition_penalty_route <- function(
    microcompass, penalty, condition_col, celltype_col) {
  route <- microcompass
  route$penalty <- penalty
  .rc_condition_penalty_comparison(
    route,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
}

