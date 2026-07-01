#' Drop cells with missing pooling labels
#' @export
rc_drop_na_grouping <- function(meta, group_cols) {
  bad <- rowSums(is.na(meta[, group_cols, drop = FALSE])) > 0
  if (any(bad)) warning(sum(bad), " cells removed due to NA in pooling columns", call. = FALSE)
  meta[!bad, , drop = FALSE]
}

#' Create sample-aware micropools
#'
#' Pools are formed within sample x cell_type, optionally further within condition.
#' This avoids treating mixed biological samples or mixed cell types as one
#' expression unit. The default target size 80 and minimum group/pool size 30 are
#' a practical compromise between single-cell sparsity and retaining rare annotated
#' cell types.
#' @export
rc_make_pools <- function(meta,
                          sample_col = "sample_id",
                          celltype_col = "cell_type",
                          condition_col = NULL,
                          target_size = 80,
                          min_pool_size = 30,
                          min_group_size = 30,
                          seed = 1) {
  if (!is.data.frame(meta)) stop("`meta` must be a data.frame.", call. = FALSE)
  if (is.null(rownames(meta)) || anyNA(rownames(meta)) || any(!nzchar(rownames(meta)))) stop("`meta` must have cell IDs in rownames(meta).", call. = FALSE)
  for (nm in c("target_size", "min_pool_size", "min_group_size")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val < 1) stop("`", nm, "` must be a single positive number.", call. = FALSE)
  }

  group_cols <- c(sample_col, celltype_col, condition_col)
  group_cols <- group_cols[!is.null(group_cols) & !is.na(group_cols) & nzchar(group_cols)]
  missing_cols <- setdiff(group_cols, colnames(meta))
  if (length(missing_cols) > 0L) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

  meta <- rc_drop_na_grouping(meta, group_cols)
  if (nrow(meta) == 0L) stop("No cells remain after dropping missing pooling labels.", call. = FALSE)
  meta$cell_id <- rownames(meta)
  keys <- interaction(meta[, group_cols, drop = FALSE], drop = TRUE, sep = "|")
  groups <- split(meta$cell_id, keys)
  set.seed(seed)

  out <- list()
  k <- 1L
  for (nm in names(groups)) {
    cells <- groups[[nm]]
    n <- length(cells)
    values <- meta[match(cells[[1]], meta$cell_id), group_cols, drop = FALSE]

    if (n < min_group_size) {
      one <- data.frame(group_id = nm, pool_id = NA_character_, cell_id = cells,
                        skipped = TRUE, skip_reason = "group_below_min_group_size",
                        low_power_pool = TRUE, pool_size = n,
                        no_within_group_pool_replicate = TRUE,
                        stringsAsFactors = FALSE)
      for (col in group_cols) one[[col]] <- values[[col]][[1]]
      out[[k]] <- one
      k <- k + 1L
      next
    }

    n_pool <- max(1L, floor(n / target_size))
    cells <- sample(cells)
    pool_assign <- rep(seq_len(n_pool), length.out = n)
    for (j in seq_len(n_pool)) {
      cc <- cells[pool_assign == j]
      one <- data.frame(group_id = nm, pool_id = paste0("pool_", k), cell_id = cc,
                        skipped = FALSE, skip_reason = NA_character_,
                        low_power_pool = length(cc) < min_pool_size,
                        pool_size = length(cc),
                        no_within_group_pool_replicate = n_pool == 1L,
                        stringsAsFactors = FALSE)
      for (col in group_cols) one[[col]] <- values[[col]][[1]]
      out[[k]] <- one
      k <- k + 1L
    }
  }

  pool_map <- do.call(rbind, out)
  rownames(pool_map) <- NULL
  pool_map
}

#' Optional seed-sensitivity helper; not used in the default workflow
#' @export
rc_make_pool_seed_replicates <- function(meta, seeds = seq_len(3), ...) {
  reps <- lapply(seeds, function(seed) {
    out <- rc_make_pools(meta = meta, seed = seed, ...)
    out$pool_seed_replicate <- seed
    idx <- !is.na(out$pool_id)
    out$pool_id[idx] <- paste0("seed", seed, "_", out$pool_id[idx])
    out
  })
  do.call(rbind, reps)
}
