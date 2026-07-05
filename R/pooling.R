#' Drop cells with NA grouping values
#'
#' @param meta Cell metadata.
#' @param group_cols Grouping columns that define pool strata.
#' @return Metadata with rows containing NA in grouping columns removed.
#' @export
rc_drop_na_grouping <- function(meta, group_cols) {
  bad <- rowSums(is.na(meta[, group_cols, drop = FALSE])) > 0
  if (any(bad)) warning(sum(bad), " cells removed due to NA in grouping columns", call. = FALSE)
  meta[!bad, , drop = FALSE]
}

#' Extract the reduction used for embedding-based pooling
#'
#' @param seu Seurat object containing PCA and/or Harmony reductions.
#' @param reduction Reduction name to use. Defaults to `harmony`, then `pca`.
#' @param dims Reduction dimensions to use.
#' @param scale_embedding Whether to scale embedding columns before pooling.
#' @return A list containing the embedding matrix and selected reduction metadata.
.rc_get_pool_embedding <- function(seu,
                                   reduction = NULL,
                                   dims = 1:30,
                                   scale_embedding = TRUE) {
  if (is.null(seu)) {
    stop("`seu` is required for embedding-based pooling.", call. = FALSE)
  }

  red_names <- names(seu@reductions)
  if (is.null(reduction)) {
    if ("harmony" %in% red_names) {
      reduction <- "harmony"
    } else if ("pca" %in% red_names) {
      reduction <- "pca"
    } else {
      stop("No `harmony` or `pca` reduction found. Run the single-cell multiomics workflow until one of these reductions is available, or set `pooling_method = 'random'` explicitly.", call. = FALSE)
    }
  }

  if (!reduction %in% red_names) {
    stop("Reduction not found: ", reduction, call. = FALSE)
  }
  if (!is.numeric(dims) || length(dims) == 0L || anyNA(dims) || any(dims < 1)) {
    stop("`dims` must contain positive numeric dimension indices.", call. = FALSE)
  }

  emb <- SeuratObject::Embeddings(seu, reduction = reduction)
  dims <- dims[dims <= ncol(emb)]
  if (length(dims) == 0L) {
    stop("No requested `dims` are available in reduction: ", reduction, call. = FALSE)
  }
  emb <- emb[, dims, drop = FALSE]

  if (scale_embedding) {
    emb <- scale(emb)
    emb[!is.finite(emb)] <- 0
  }

  list(
    emb = emb,
    reduction = reduction,
    dims = dims
  )
}

#' Randomly assign cells to pools
#'
#' @param cells Cell identifiers in one pooling stratum.
#' @param n_pool Number of pools to create.
#' @return A list of cell identifier vectors, one per pool.
.rc_assign_random_pools <- function(cells, n_pool) {
  cells <- sample(cells)
  split(cells, rep(seq_len(n_pool), length.out = length(cells)))
}

#' Assign cells to pools by embedding-space proximity
#'
#' @param cells Cell identifiers in one pooling stratum.
#' @param emb Embedding matrix with cells in row names.
#' @param n_pool Number of pools to create.
#' @param nstart Number of random starts passed to `stats::kmeans()`.
#' @return A list of cell identifier vectors, one per pool.
.rc_assign_embedding_pools <- function(cells,
                                       emb,
                                       n_pool,
                                       nstart = 20) {
  x <- emb[cells, , drop = FALSE]

  if (n_pool == 1L) {
    return(list(cells))
  }

  km <- stats::kmeans(
    x = x,
    centers = n_pool,
    nstart = nstart,
    iter.max = 100
  )

  split(cells, km$cluster)
}

#' Create sample-aware micropools from cell metadata
#'
#' Pools are formed only within the requested sample/condition/cell-type/state
#' strata. If `target_celltype` is supplied without a contrast, only that
#' cell type is pooled by default; set `include_other_celltypes_as_control = TRUE`
#' to keep all other cell types as an explicit control group.
#' Groups below `min_group_size` are represented as skipped rows and are
#' excluded from pseudobulk by default; groups with one valid pool are retained
#' and flagged as lacking within-group pool replication. With `pooling_method =
#' "auto"`, state-aware stratification is preferred when `state_col` is supplied,
#' otherwise Harmony/PCA embeddings from `seu` are used when available, and the
#' original random pooling logic is used only when no embedding input is supplied.
#' @export
rc_make_pools <- function(meta,
                          seu = NULL,
                          sample_col = "sample_id",
                          celltype_col = "cell_type",
                          condition_col = NULL,
                          state_col = NULL,
                          target_size = 80,
                          min_pool_size = 30,
                          min_group_size = 30,
                          min_size = NULL,
                          seed = 1,
                          pooling_method = c("auto", "random", "state", "embedding"),
                          target_celltype = NULL,
                          include_other_celltypes_as_control = FALSE,
                          contrast_col = NULL,
                          target_contrast_label = NULL,
                          other_contrast_label = "other",
                          reduction = NULL,
                          dims = 1:30,
                          nstart = 20,
                          scale_embedding = TRUE,
                          state_source = NA_character_,
                          state_resolution = NA_character_,
                          BPPARAM = NULL) {
  if (!is.null(min_size)) {
    min_pool_size <- min_size
    min_group_size <- min_size
  }
  pooling_method <- match.arg(pooling_method)
  if (pooling_method == "auto") {
    if (!is.null(state_col)) {
      pooling_method <- "state"
    } else if (!is.null(seu)) {
      pooling_method <- "embedding"
    } else {
      pooling_method <- "random"
    }
  }
  if (!is.data.frame(meta)) stop("`meta` must be a data.frame.", call. = FALSE)
  if (!is.logical(include_other_celltypes_as_control) || length(include_other_celltypes_as_control) != 1L || is.na(include_other_celltypes_as_control)) {
    stop("`include_other_celltypes_as_control` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(contrast_col)) contrast_col <- "celltype_contrast"
  if (!is.character(contrast_col) || length(contrast_col) != 1L || is.na(contrast_col) || !nzchar(contrast_col)) {
    stop("`contrast_col` must be a non-empty column name.", call. = FALSE)
  }
  if (!is.character(other_contrast_label) || length(other_contrast_label) != 1L || is.na(other_contrast_label) || !nzchar(other_contrast_label)) {
    stop("`other_contrast_label` must be a non-empty label.", call. = FALSE)
  }
  if (is.null(rownames(meta)) || anyNA(rownames(meta)) || any(!nzchar(rownames(meta)))) {
    stop("`meta` must have non-empty cell IDs in rownames(meta).", call. = FALSE)
  }
  for (nm in c("target_size", "min_pool_size", "min_group_size")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val < 1) stop("`", nm, "` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(nstart) || length(nstart) != 1L || is.na(nstart) || nstart < 1) {
    stop("`nstart` must be a single positive number.", call. = FALSE)
  }

  required_cols <- c(sample_col, condition_col, celltype_col, state_col)
  required_cols <- required_cols[!is.null(required_cols) & !is.na(required_cols) & nzchar(required_cols)]
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

  if (!is.null(target_celltype)) {
    target_celltype <- unique(as.character(target_celltype))
    target_celltype <- target_celltype[!is.na(target_celltype) & nzchar(target_celltype)]
    if (length(target_celltype) == 0L) stop("`target_celltype` must contain at least one non-empty cell-type label.", call. = FALSE)
    celltype_values <- as.character(meta[[celltype_col]])
    target_hit <- celltype_values %in% target_celltype
    if (!any(target_hit, na.rm = TRUE)) {
      stop("`target_celltype` did not match any values in `", celltype_col, "`.", call. = FALSE)
    }
    target_label <- target_contrast_label
    if (is.null(target_label)) target_label <- paste(target_celltype, collapse = "+")
    if (!is.character(target_label) || length(target_label) != 1L || is.na(target_label) || !nzchar(target_label)) {
      stop("`target_contrast_label` must be NULL or a non-empty label.", call. = FALSE)
    }

    if (isTRUE(include_other_celltypes_as_control)) {
      original_col <- paste0("original_", celltype_col)
      if (!original_col %in% colnames(meta)) meta[[original_col]] <- celltype_values
      meta[[celltype_col]] <- ifelse(target_hit, target_label, other_contrast_label)
      if (is.null(condition_col)) {
        meta[[contrast_col]] <- ifelse(target_hit, target_label, other_contrast_label)
        condition_col <- contrast_col
      }
    } else {
      meta <- meta[target_hit, , drop = FALSE]
      if (nrow(meta) == 0L) stop("No cells remain after selecting `target_celltype`.", call. = FALSE)
    }
  }

  group_cols <- c(sample_col, condition_col, celltype_col, state_col)
  group_cols <- group_cols[!is.null(group_cols) & !is.na(group_cols) & nzchar(group_cols)]

  emb <- NULL
  used_reduction <- NA_character_
  used_dims <- NA_character_
  if (pooling_method == "embedding") {
    emb_info <- .rc_get_pool_embedding(
      seu = seu,
      reduction = reduction,
      dims = dims,
      scale_embedding = scale_embedding
    )
    emb <- emb_info$emb
    used_reduction <- emb_info$reduction
    used_dims <- paste(emb_info$dims, collapse = ",")

    common <- intersect(rownames(meta), rownames(emb))
    if (length(common) == 0L) {
      stop("No cells overlap between `meta` rownames and the selected embedding rownames.", call. = FALSE)
    }
    meta <- meta[common, , drop = FALSE]
    emb <- emb[common, , drop = FALSE]
  }

  meta <- rc_drop_na_grouping(meta, group_cols)
  meta$cell_id <- rownames(meta)
  if (nrow(meta) == 0L) stop("No cells remain after dropping NA grouping values.", call. = FALSE)
  keys <- interaction(meta[, group_cols, drop = FALSE], drop = TRUE, sep = "|")
  groups <- split(meta$cell_id, keys)
  set.seed(seed)

  out <- list(); k <- 1L
  for (nm in names(groups)) {
    cells <- groups[[nm]]; n <- length(cells)
    values <- meta[match(cells[[1]], meta$cell_id), group_cols, drop = FALSE]
    if (n < min_group_size) {
      out[[k]] <- data.frame(group_id = nm, group_key = nm, pool_id = NA_character_, cell_id = cells,
                             skipped = TRUE, skip_reason = "group_below_min_group_size",
                             low_power_pool = TRUE, pool_size = n,
                             no_within_group_pool_replicate = TRUE,
                             pooling_method = pooling_method, pool_reduction = used_reduction,
                             pool_dims = used_dims, pool_seed = seed,
                             state_source = state_source, state_resolution = state_resolution, stringsAsFactors = FALSE)
      for (col in group_cols) out[[k]][[col]] <- values[[col]][[1]]
      k <- k + 1L; next
    }
    n_pool <- max(1L, floor(n / target_size))
    if (pooling_method %in% c("random", "state")) {
      pool_list <- .rc_assign_random_pools(cells, n_pool)
    } else if (pooling_method == "embedding") {
      pool_list <- .rc_assign_embedding_pools(cells = cells, emb = emb, n_pool = n_pool, nstart = nstart)
    }
    for (j in seq_along(pool_list)) {
      cc <- pool_list[[j]]
      out[[k]] <- data.frame(group_id = nm, group_key = nm, pool_id = paste0("pool_", k), cell_id = cc,
                             skipped = FALSE, skip_reason = NA_character_,
                             low_power_pool = length(cc) < min_pool_size, pool_size = length(cc),
                             no_within_group_pool_replicate = length(pool_list) == 1L,
                             pooling_method = pooling_method, pool_reduction = used_reduction,
                             pool_dims = used_dims, pool_seed = seed,
                             state_source = state_source, state_resolution = state_resolution, stringsAsFactors = FALSE)
      for (col in group_cols) out[[k]][[col]] <- values[[col]][[1]]
      k <- k + 1L
    }
  }
  pool_map <- do.call(rbind, out)
  rownames(pool_map) <- NULL
  pool_map
}


#' Create multiple pool maps across random seeds for seed sensitivity
#' @export
rc_make_pool_seed_replicates <- function(meta, seeds = seq_len(5), ...) {
  reps <- lapply(seeds, function(seed) {
    out <- rc_make_pools(meta = meta, seed = seed, ...)
    non_na <- !is.na(out$pool_id)
    out$pool_id[non_na] <- paste0("seed", seed, "_", out$pool_id[non_na])
    out$pool_seed_replicate <- seed
    out
  })
  do.call(rbind, reps)
}
