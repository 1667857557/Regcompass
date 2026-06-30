#' Create sample-aware micropools from cell metadata
#'
#' `rc_make_pools()` assigns cells to micropools within sample-aware strata. Cells
#' from different samples are never mixed. If `condition_col` is supplied, pools
#' are formed within `sample_id x condition x cell_type`; if `state_col` is also
#' supplied, the state/cluster column is added as a local-state stratum. If no
#' state column is supplied, v0.2 uses `sample_id x cell_type` (plus condition
#' when requested), matching the simplified first implementation in the
#' development specification.
#'
#' @param meta A cell metadata `data.frame` with cell IDs in `rownames(meta)`.
#' @param sample_col Metadata column containing sample identifiers.
#' @param celltype_col Metadata column containing cell type annotations.
#' @param condition_col Optional metadata column containing biological condition.
#' @param state_col Optional metadata column containing an existing cluster/local
#' state, for example `"seurat_clusters"`.
#' @param target_size Target number of cells per pool.
#' @param min_size Minimum informative pool size. Groups smaller than this value
#' are retained as one low-power pool; split pools below this value are flagged.
#' @param seed Random seed used before shuffling cells within each stratum.
#' @param BPPARAM Optional `BiocParallelParam` object. When provided and
#' BiocParallel is installed, independent strata are processed with
#' `BiocParallel::bplapply()`; otherwise base `lapply()` is used.
#'
#' @return A `data.frame` with one row per cell and columns `pool_id`, `cell_id`,
#' `low_power_pool`, `pool_size`, `group_key`, and the grouping metadata columns.
#' @export
rc_make_pools <- function(meta,
                          sample_col = "sample_id",
                          celltype_col = "cell_type",
                          condition_col = NULL,
                          state_col = NULL,
                          target_size = 80,
                          min_size = 30,
                          seed = 1,
                          BPPARAM = NULL) {
  if (!is.data.frame(meta)) {
    stop("`meta` must be a data.frame.", call. = FALSE)
  }
  if (is.null(rownames(meta)) || anyNA(rownames(meta)) || any(!nzchar(rownames(meta)))) {
    stop("`meta` must have non-empty cell IDs in rownames(meta).", call. = FALSE)
  }
  if (!is.numeric(target_size) || length(target_size) != 1L || is.na(target_size) || target_size < 1) {
    stop("`target_size` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(min_size) || length(min_size) != 1L || is.na(min_size) || min_size < 1) {
    stop("`min_size` must be a single positive number.", call. = FALSE)
  }

  group_cols <- c(sample_col, condition_col, celltype_col, state_col)
  group_cols <- group_cols[!is.na(group_cols) & nzchar(group_cols)]
  missing_cols <- setdiff(group_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  df <- meta
  df$cell_id <- rownames(df)
  split_key <- interaction(df[, group_cols, drop = FALSE], drop = TRUE, sep = "|")
  groups <- split(df$cell_id, split_key)
  group_values <- lapply(groups, function(cells) df[match(cells[[1]], df$cell_id), group_cols, drop = FALSE])

  set.seed(seed)
  group_seeds <- sample.int(.Machine$integer.max, length(groups))

  make_group <- function(i) {
    cells <- groups[[i]]
    n <- length(cells)
    values <- group_values[[i]]
    group_key <- names(groups)[[i]]

    if (n < min_size) {
      pool_index <- rep(1L, n)
    } else {
      set.seed(group_seeds[[i]])
      cells <- sample(cells)
      n_pool <- ceiling(n / target_size)
      pool_index <- rep(seq_len(n_pool), length.out = n)
    }

    pieces <- lapply(sort(unique(pool_index)), function(j) {
      cc <- cells[pool_index == j]
      out <- data.frame(
        pool_index = j,
        cell_id = cc,
        low_power_pool = length(cc) < min_size,
        pool_size = length(cc),
        group_key = group_key,
        stringsAsFactors = FALSE
      )
      for (col in group_cols) {
        out[[col]] <- values[[col]][[1]]
      }
      out
    })
    do.call(rbind, pieces)
  }

  pool_pieces <- if (!is.null(BPPARAM)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("BiocParallel must be installed when `BPPARAM` is provided.", call. = FALSE)
    }
    BiocParallel::bplapply(seq_along(groups), make_group, BPPARAM = BPPARAM)
  } else {
    lapply(seq_along(groups), make_group)
  }

  pool_map <- do.call(rbind, pool_pieces)
  pool_map$pool_id <- paste0("pool_", match(
    paste(pool_map$group_key, pool_map$pool_index, sep = "::"),
    unique(paste(pool_map$group_key, pool_map$pool_index, sep = "::"))
  ))
  pool_map$pool_index <- NULL
  pool_map <- pool_map[, c("pool_id", "cell_id", "low_power_pool", "pool_size", "group_key", group_cols), drop = FALSE]
  rownames(pool_map) <- NULL
  pool_map
}
