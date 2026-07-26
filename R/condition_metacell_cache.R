.rc_condition_metacell_md5 <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2, compress = FALSE)
  unname(as.character(tools::md5sum(path)))
}

.rc_condition_metacell_projection_weights <- function(n, a, b) {
  index <- as.numeric(seq_len(n))
  ((a * index + b * index^2) %% 1000003) + 1
}

.rc_condition_metacell_matrix_fingerprint <- function(x) {
  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop("Metacell cache fingerprinting requires a two-dimensional matrix.",
         call. = FALSE)
  }
  row_weight_a <- .rc_condition_metacell_projection_weights(nrow(x), 104729, 37)
  row_weight_b <- .rc_condition_metacell_projection_weights(nrow(x), 130363, 53)
  col_weight_a <- .rc_condition_metacell_projection_weights(ncol(x), 155921, 71)
  col_weight_b <- .rc_condition_metacell_projection_weights(ncol(x), 196613, 89)
  value_projection <- list(
    rows_a = as.numeric(x %*% col_weight_a),
    rows_b = as.numeric(x %*% col_weight_b),
    columns_a = as.numeric(Matrix::crossprod(x, row_weight_a)),
    columns_b = as.numeric(Matrix::crossprod(x, row_weight_b))
  )
  list(
    dim = as.integer(dim(x)),
    row_ids_md5 = .rc_condition_metacell_md5(as.character(rownames(x))),
    col_ids_md5 = .rc_condition_metacell_md5(as.character(colnames(x))),
    nnzero = as.numeric(Matrix::nnzero(x)),
    row_sums_md5 = .rc_condition_metacell_md5(as.numeric(Matrix::rowSums(x))),
    col_sums_md5 = .rc_condition_metacell_md5(as.numeric(Matrix::colSums(x))),
    values_md5 = .rc_condition_metacell_md5(value_projection)
  )
}

.rc_condition_metacell_reduction_fingerprint <- function(
    object, reduction, dims, cells) {
  reduction <- trimws(as.character(reduction[[1L]]))
  dims <- as.integer(dims)
  if (!nzchar(reduction) || !length(dims) || anyNA(dims) || any(dims < 1L)) {
    stop("Metacell reductions require a non-empty name and positive dimensions.",
         call. = FALSE)
  }
  if (!reduction %in% names(object@reductions)) {
    stop("Required metacell reduction is absent: `", reduction, "`.",
         call. = FALSE)
  }
  embeddings <- SeuratObject::Embeddings(object[[reduction]])
  if (max(dims) > ncol(embeddings)) {
    stop(
      "Reduction `", reduction, "` contains only ", ncol(embeddings),
      " dimensions but the metacell contract requests dimension ", max(dims),
      ".", call. = FALSE
    )
  }
  index <- match(cells, rownames(embeddings))
  if (anyNA(index)) {
    stop("Reduction `", reduction, "` lacks input cells required by Stage 2.",
         call. = FALSE)
  }
  selected <- embeddings[index, dims, drop = FALSE]
  rownames(selected) <- cells
  list(
    reduction = reduction,
    dims = dims,
    embedding = .rc_condition_metacell_matrix_fingerprint(selected)
  )
}

.rc_condition_metacell_cache_contract <- function(
    object, condition_col, celltype_col, rna_assay, atac_assay,
    metacell_args) {
  cells <- as.character(colnames(object))
  meta_index <- match(cells, rownames(object@meta.data))
  if (anyNA(meta_index)) {
    stop("Seurat cell order cannot be aligned to metadata for cache validation.",
         call. = FALSE)
  }
  defaults <- list(
    rna_reduction = "pca", atac_reduction = "lsi",
    rna_dims = 1:30, atac_dims = 2:30, gamma = 30L, seed = 12345L,
    min_cells_per_stratum = 100L, min_metacell_size = 20L,
    min_metacells_per_stratum = 2L
  )
  supplied <- metacell_args[intersect(names(metacell_args), names(defaults))]
  analysis_args <- modifyList(defaults, supplied)
  for (field in c(
    "gamma", "seed", "min_cells_per_stratum", "min_metacell_size",
    "min_metacells_per_stratum"
  )) analysis_args[[field]] <- as.integer(analysis_args[[field]])
  analysis_args$rna_dims <- as.integer(analysis_args$rna_dims)
  analysis_args$atac_dims <- as.integer(analysis_args$atac_dims)
  meta_signature <- data.frame(
    cell_id = cells,
    condition = as.character(object@meta.data[[condition_col]][meta_index]),
    cell_type = as.character(object@meta.data[[celltype_col]][meta_index]),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = "regcompass_condition_metacell_cache_v2_supercell_label",
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    strata_cols = condition_col,
    label_col = celltype_col,
    ordered_cell_metadata_md5 = .rc_condition_metacell_md5(meta_signature),
    analysis_args = analysis_args,
    rna_counts = .rc_condition_metacell_matrix_fingerprint(
      .rc_get_assay_counts(object, rna_assay)
    ),
    atac_counts = .rc_condition_metacell_matrix_fingerprint(
      .rc_get_assay_counts(object, atac_assay)
    ),
    rna_reduction = .rc_condition_metacell_reduction_fingerprint(
      object, analysis_args$rna_reduction, analysis_args$rna_dims, cells
    ),
    atac_reduction = .rc_condition_metacell_reduction_fingerprint(
      object, analysis_args$atac_reduction, analysis_args$atac_dims, cells
    )
  )
}

.rc_condition_metacell_has_checkpoints <- function(outdir) {
  if (!dir.exists(outdir)) return(FALSE)
  files <- list.files(
    outdir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE
  )
  any(basename(files) %in% c(
    "metacell_metadata.tsv.gz", "rna_counts.rds", "atac_counts.rds"
  ))
}

.rc_validate_condition_metacell_cache <- function(
    outdir, contract, overwrite = FALSE) {
  if (!.rc_condition_metacell_has_checkpoints(outdir) || isTRUE(overwrite)) {
    return(invisible(FALSE))
  }
  contract_file <- file.path(outdir, "condition_metacell_cache_contract.rds")
  if (!file.exists(contract_file)) {
    stop(
      "Existing condition-metacell checkpoints predate the current cache ",
      "contract. Set `metacell_args$overwrite = TRUE` to rebuild them.",
      call. = FALSE
    )
  }
  observed <- tryCatch(readRDS(contract_file), error = function(error) {
    stop(
      "The existing condition-metacell cache contract cannot be read: ",
      conditionMessage(error), ". Rebuild with `overwrite = TRUE`.",
      call. = FALSE
    )
  })
  if (!identical(observed, contract)) {
    stop(
      "Existing metacell checkpoints use different cells, labels, assay ",
      "contents, reductions, or construction parameters. Rebuild with ",
      "`metacell_args$overwrite = TRUE`.", call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_assign_metacell_dominant_celltype <- function(
    pooled, object, celltype_col) {
  membership <- pooled$membership
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership))) {
    stop("Condition-only metacells require cell-to-metacell membership.",
         call. = FALSE)
  }
  cell_index <- match(as.character(membership$cell_id), rownames(object@meta.data))
  if (anyNA(cell_index)) {
    stop("Metacell membership contains cells absent from the input object.",
         call. = FALSE)
  }
  cell_type <- trimws(as.character(object@meta.data[[celltype_col]][cell_index]))
  if (anyNA(cell_type) || any(!nzchar(cell_type))) {
    stop("Cell-type metadata are incomplete in metacell membership.",
         call. = FALSE)
  }
  membership[[celltype_col]] <- cell_type
  composition <- stats::aggregate(
    rep.int(1L, nrow(membership)),
    by = list(
      metacell_id = as.character(membership$metacell_id),
      cell_type = cell_type
    ), FUN = sum
  )
  colnames(composition) <- c("metacell_id", celltype_col, "n_cells")
  split_rows <- split(seq_len(nrow(composition)), composition$metacell_id)
  summary <- do.call(rbind, lapply(names(split_rows), function(id) {
    z <- composition[split_rows[[id]], , drop = FALSE]
    z <- z[order(-z$n_cells, as.character(z[[celltype_col]])), , drop = FALSE]
    total <- sum(z$n_cells)
    top_count <- max(z$n_cells)
    tied <- sum(z$n_cells == top_count) > 1L
    data.frame(
      metacell_id = id,
      dominant_celltype = if (tied) NA_character_ else
        as.character(z[[celltype_col]][[1L]]),
      dominant_celltype_fraction = top_count / total,
      n_celltypes = nrow(z),
      mixed_celltype_metacell = nrow(z) > 1L,
      dominant_celltype_tied = tied,
      stringsAsFactors = FALSE
    )
  }))
  tied_ids <- as.character(summary$metacell_id[
    summary$dominant_celltype_tied %in% TRUE
  ])
  if (length(tied_ids)) {
    stop(
      "SuperCell2 produced metacells with tied cell-type labels: ",
      paste(utils::head(tied_ids, 10L), collapse = ", "), call. = FALSE
    )
  }
  meta <- pooled$metacell_meta
  index <- match(as.character(meta$metacell_id), summary$metacell_id)
  if (anyNA(index)) {
    stop("Cell-type summaries do not align with metacell metadata.",
         call. = FALSE)
  }
  meta[[celltype_col]] <- summary$dominant_celltype[index]
  meta$dominant_celltype_fraction <- summary$dominant_celltype_fraction[index]
  meta$n_celltypes <- summary$n_celltypes[index]
  meta$mixed_celltype_metacell <- summary$mixed_celltype_metacell[index]
  meta$dominant_celltype_tied <- summary$dominant_celltype_tied[index]
  pooled$membership <- membership
  pooled$metacell_meta <- meta
  pooled$celltype_composition <- composition
  pooled$celltype_composition_summary <- summary
  pooled
}
