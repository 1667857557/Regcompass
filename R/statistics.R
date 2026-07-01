#' Aggregate pool-level scores to sample-by-cell-type summaries
#'
#' This is the default v0.7 statistical unit: downstream differential tests should
#' use biological sample-level aggregates rather than treating pools as
#' independent replicates.
#'
#' @param score_mat Reaction-by-pool numeric matrix.
#' @param pool_meta Data frame with one row per pool and at least `pool_id`,
#' `sample_col`, and `celltype_col` columns.
#' @param sample_col Column in `pool_meta` containing biological sample IDs.
#' @param celltype_col Column in `pool_meta` containing annotated cell types.
#'
#' @return Reaction-by-sample-celltype matrix of row medians across pools.
#' @export
rc_sample_aggregate <- function(score_mat,
                                pool_meta,
                                sample_col = "sample_id",
                                celltype_col = "cell_type") {
  score_mat <- as.matrix(score_mat)
  storage.mode(score_mat) <- "numeric"
  if (is.null(rownames(score_mat)) || is.null(colnames(score_mat))) {
    stop("`score_mat` must have reaction row names and pool column names.", call. = FALSE)
  }
  if (!is.data.frame(pool_meta)) stop("`pool_meta` must be a data.frame.", call. = FALSE)
  required <- c("pool_id", sample_col, celltype_col)
  missing <- setdiff(required, colnames(pool_meta))
  if (length(missing) > 0L) stop("`pool_meta` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyNA(pool_meta$pool_id) || anyDuplicated(pool_meta$pool_id)) {
    stop("`pool_meta$pool_id` must be non-missing and unique.", call. = FALSE)
  }
  if (anyNA(pool_meta[[sample_col]]) || anyNA(pool_meta[[celltype_col]])) {
    stop("Sample and cell-type columns in `pool_meta` must not contain missing values.", call. = FALSE)
  }
  missing_pools <- setdiff(pool_meta$pool_id, colnames(score_mat))
  if (length(missing_pools) > 0L) {
    stop("`score_mat` is missing pools from `pool_meta`: ", paste(utils::head(missing_pools, 5L), collapse = ", "), call. = FALSE)
  }

  groups <- interaction(pool_meta[[sample_col]], pool_meta[[celltype_col]], drop = TRUE)
  split_pools <- split(as.character(pool_meta$pool_id), groups)
  res <- lapply(split_pools, function(pools) {
    matrixStats::rowMedians(score_mat[, pools, drop = FALSE], na.rm = TRUE)
  })
  out <- do.call(cbind, res)
  rownames(out) <- rownames(score_mat)
  colnames(out) <- names(split_pools)
  out
}

#' Fit simple sample-level linear models for each reaction
#'
#' v0.7 intentionally uses ordinary sample-level linear models. Pool-level tests
#' and mixed models are not used here; mixed models are reserved for later
#' milestones when sample size and repeated-pool structure justify them.
#'
#' @param Y Reaction-by-sample aggregate matrix. Columns must be named by sample
#' aggregate IDs present in `rownames(sample_meta)`.
#' @param sample_meta Data frame with one row per sample aggregate.
#' @param formula_str Right-hand-side model formula string beginning with `~`, for
#' example `"~ condition + batch"`.
#' @param BPPARAM Optional `BiocParallelParam` for reaction-level parallelism.
#'
#' @return Data frame with reaction ID, model term, estimate, standard error,
#' statistic, p-value, and BH q-value computed within each term.
#' @export
rc_lm_by_reaction <- function(Y,
                              sample_meta,
                              formula_str,
                              BPPARAM = NULL) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "numeric"
  if (is.null(rownames(Y)) || is.null(colnames(Y))) stop("`Y` must have reaction row names and sample columns.", call. = FALSE)
  if (!is.data.frame(sample_meta)) stop("`sample_meta` must be a data.frame.", call. = FALSE)
  if (is.null(rownames(sample_meta))) stop("`sample_meta` must have row names matching `colnames(Y)`.", call. = FALSE)
  missing_samples <- setdiff(colnames(Y), rownames(sample_meta))
  if (length(missing_samples) > 0L) {
    stop("`sample_meta` is missing rows for samples: ", paste(utils::head(missing_samples, 5L), collapse = ", "), call. = FALSE)
  }
  if (!is.character(formula_str) || length(formula_str) != 1L || !grepl("^\\s*~", formula_str)) {
    stop("`formula_str` must be a single right-hand-side formula string beginning with '~'.", call. = FALSE)
  }
  f <- stats::as.formula(paste("y", formula_str))
  sample_meta <- sample_meta[colnames(Y), , drop = FALSE]

  pieces <- rc_parallel_lapply(rownames(Y), function(r) {
    df <- sample_meta
    df$y <- as.numeric(Y[r, rownames(sample_meta)])
    fit <- stats::lm(f, data = df, na.action = stats::na.exclude)
    co <- summary(fit)$coefficients
    data.frame(
      reaction_id = r,
      term = rownames(co),
      estimate = co[, 1],
      std_error = co[, 2],
      statistic = co[, 3],
      p_value = co[, 4],
      n = stats::nobs(fit),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, BPPARAM = BPPARAM)
  out <- do.call(rbind, pieces)
  out$q_value <- ave(out$p_value, out$term, FUN = function(p) stats::p.adjust(p, method = "BH"))
  rownames(out) <- NULL
  out
}
