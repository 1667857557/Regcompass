.rc_condition_celltype_pool_col <- function(meta) {
  candidate <- ".rc_condition_pool_id"
  while (candidate %in% colnames(meta)) candidate <- paste0(candidate, "_")
  candidate
}

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
    row_sums_md5 = .rc_condition_metacell_md5(
      as.numeric(Matrix::rowSums(x))
    ),
    col_sums_md5 = .rc_condition_metacell_md5(
      as.numeric(Matrix::colSums(x))
    ),
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
  analysis_defaults <- list(
    rna_reduction = "pca",
    atac_reduction = "lsi",
    rna_dims = 1:30,
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 100L,
    min_metacell_size = 20L,
    min_metacells_per_stratum = 2L
  )
  supplied <- metacell_args[
    intersect(names(metacell_args), names(analysis_defaults))
  ]
  analysis_args <- modifyList(analysis_defaults, supplied)
  integer_fields <- c(
    "gamma", "seed", "min_cells_per_stratum", "min_metacell_size",
    "min_metacells_per_stratum"
  )
  for (field in integer_fields) {
    analysis_args[[field]] <- as.integer(analysis_args[[field]])
  }
  analysis_args$rna_dims <- as.integer(analysis_args$rna_dims)
  analysis_args$atac_dims <- as.integer(analysis_args$atac_dims)
  meta_signature <- data.frame(
    cell_id = cells,
    condition = as.character(
      object@meta.data[[condition_col]][meta_index]
    ),
    cell_type = as.character(
      object@meta.data[[celltype_col]][meta_index]
    ),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = "regcompass_condition_celltype_metacell_cache_v2",
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    label_col = NULL,
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
    outdir, recursive = TRUE, full.names = TRUE, all.files = FALSE,
    no.. = TRUE
  )
  any(basename(files) %in% c(
    "metacell_metadata.tsv.gz", "rna_counts.rds", "atac_counts.rds"
  ))
}

.rc_validate_condition_metacell_cache <- function(
    outdir, contract, overwrite = FALSE) {
  has_checkpoints <- .rc_condition_metacell_has_checkpoints(outdir)
  if (!has_checkpoints || isTRUE(overwrite)) return(invisible(FALSE))
  contract_file <- file.path(outdir, "condition_metacell_cache_contract.rds")
  if (!file.exists(contract_file)) {
    stop(
      "Existing condition-metacell checkpoints predate the audited cache contract. ",
      "Set `metacell_args$overwrite = TRUE` to rebuild them before reuse.",
      call. = FALSE
    )
  }
  observed <- tryCatch(
    readRDS(contract_file),
    error = function(e) {
      stop(
        "The existing condition-metacell cache contract cannot be read: ",
        conditionMessage(e),
        ". Set `metacell_args$overwrite = TRUE` to rebuild the checkpoints.",
        call. = FALSE
      )
    }
  )
  if (!identical(observed, contract)) {
    stop(
      "Existing condition-metacell checkpoints were created with different ",
      "cells, labels, assay contents, reduction embeddings, or construction ",
      "parameters. Set `metacell_args$overwrite = TRUE` to rebuild them.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_make_condition_pooled_metacells <- function(
    object, outdir,
    sample_col = NULL,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list()) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
  }
  required <- unique(c(condition_col, celltype_col, sample_col))
  required <- required[!is.na(required) & nzchar(required)]
  missing <- setdiff(required, colnames(object@meta.data))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invalid <- vapply(
    object@meta.data[, required, drop = FALSE],
    function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))),
    logical(1)
  )
  if (any(invalid)) {
    stop("Condition and cell-type metadata must be complete.",
         call. = FALSE)
  }
  if (!identical(fragment_files, FALSE) && !is.null(fragment_files)) {
    stop(
      paste(
        "The canonical condition-by-cell-type path requires",
        "`fragment_files = FALSE` and aggregates the existing ATAC",
        "peak-count assay."
      ),
      call. = FALSE
    )
  }
  unsupported <- intersect(
    names(metacell_args), c("sample_balance", "sample_balance_seed")
  )
  if (length(unsupported)) {
    stop(
      "Sample balancing is not part of the canonical workflow: ",
      paste(unsupported, collapse = ", "), call. = FALSE
    )
  }
  reserved <- intersect(names(metacell_args), c(
    "object", "outdir", "sample_col", "condition_col", "celltype_col",
    "label_col", "rna_assay", "atac_assay", "fragment_files",
    "save_metacell_object", "save_counts", "save_fragments",
    "require_fragment_aggregation", "fragment_aggregation_backend",
    "on_stratum_error"
  ))
  if (length(reserved)) {
    stop(
      "`metacell_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "), call. = FALSE
    )
  }
  if (is.null(metacell_args$gamma)) metacell_args$gamma <- 30L
  cache_contract <- .rc_condition_metacell_cache_contract(
    object = object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    metacell_args = metacell_args
  )
  .rc_validate_condition_metacell_cache(
    outdir = outdir,
    contract = cache_contract,
    overwrite = isTRUE(metacell_args$overwrite %||% FALSE)
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    cache_contract,
    file.path(outdir, "condition_metacell_cache_contract.rds")
  )
  internal_pool_col <- .rc_condition_celltype_pool_col(object@meta.data)
  object@meta.data[[internal_pool_col]] <- paste0(
    as.character(object@meta.data[[condition_col]]), "__",
    as.character(object@meta.data[[celltype_col]]),
    "__condition_celltype_pool"
  )
  defaults <- list(
    object = object,
    outdir = outdir,
    sample_col = internal_pool_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    label_col = NULL,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = FALSE,
    save_metacell_object = TRUE,
    save_counts = TRUE,
    save_fragments = FALSE,
    require_fragment_aggregation = FALSE,
    fragment_aggregation_backend = "none",
    on_stratum_error = "stop"
  )
  defaults[names(metacell_args)] <- NULL
  pooled <- do.call(
    rc_make_supercell2_metacells,
    c(defaults, metacell_args)
  )
  meta <- pooled$metacell_meta
  if (!is.data.frame(meta) || !nrow(meta)) {
    stop(
      "Condition-by-cell-type SuperCell2 produced no metacells.",
      call. = FALSE
    )
  }
  meta$condition_celltype_pool_id <- meta[[internal_pool_col]]
  cell_index <- match(
    as.character(pooled$membership$cell_id), rownames(object@meta.data)
  )
  if (anyNA(cell_index)) {
    stop("Metacell membership cannot be aligned to input cell metadata.",
         call. = FALSE)
  }
  pooled$membership[[condition_col]] <- as.character(
    object@meta.data[[condition_col]][cell_index]
  )
  pooled$membership[[celltype_col]] <- as.character(
    object@meta.data[[celltype_col]][cell_index]
  )
  membership_groups <- split(
    seq_len(nrow(pooled$membership)),
    as.character(pooled$membership$metacell_id)
  )
  impure <- names(membership_groups)[vapply(
    membership_groups, function(rows) {
      length(unique(pooled$membership[[condition_col]][rows])) != 1L ||
        length(unique(pooled$membership[[celltype_col]][rows])) != 1L
    }, logical(1)
  )]
  if (length(impure)) {
    stop(
      "SuperCell returned metacells mixing condition or broad cell type: ",
      paste(utils::head(impure, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  meta_index <- match(
    as.character(meta$metacell_id), names(membership_groups)
  )
  if (anyNA(meta_index)) {
    stop("SuperCell metadata do not cover every membership group.",
         call. = FALSE)
  }
  observed_condition <- vapply(
    membership_groups, function(rows) {
      pooled$membership[[condition_col]][rows[[1L]]]
    }, character(1)
  )
  observed_celltype <- vapply(
    membership_groups, function(rows) {
      pooled$membership[[celltype_col]][rows[[1L]]]
    }, character(1)
  )
  if (any(
    as.character(meta[[condition_col]]) !=
      observed_condition[as.character(meta$metacell_id)] |
    as.character(meta[[celltype_col]]) !=
      observed_celltype[as.character(meta$metacell_id)]
  )) {
    stop("SuperCell metadata disagree with membership hard strata.",
         call. = FALSE)
  }
  if (!is.null(sample_col) && nzchar(sample_col)) {
    pooled$membership[[sample_col]] <- as.character(
      object@meta.data[[sample_col]][cell_index]
    )
    pooled$sample_composition <- stats::aggregate(
      rep.int(1L, nrow(pooled$membership)),
      by = list(
        metacell_id = as.character(pooled$membership$metacell_id),
        sample_id = as.character(pooled$membership[[sample_col]])
      ),
      FUN = sum
    )
    colnames(pooled$sample_composition)[[3L]] <- "n_cells"
  } else {
    pooled$sample_composition <- data.frame()
  }
  pooled$celltype_composition <- data.frame(
    metacell_id = as.character(meta$metacell_id),
    value = as.character(meta[[celltype_col]]),
    n_cells = as.integer(meta$n_cells),
    stringsAsFactors = FALSE
  )
  colnames(pooled$celltype_composition)[[2L]] <- celltype_col
  pooled$celltype_composition_summary <- data.frame(
    metacell_id = as.character(meta$metacell_id),
    dominant_celltype = as.character(meta[[celltype_col]]),
    dominant_celltype_fraction = 1,
    n_celltypes = 1L,
    mixed_celltype_metacell = FALSE,
    dominant_celltype_tied = FALSE,
    stringsAsFactors = FALSE
  )
  meta$pooling_scope <- "condition_by_cell_type"
  meta$sample_weighting <- "none"
  meta$sample_col_role <- "internal_condition_celltype_pool_id"
  meta$celltype_role <- "hard_stratum"
  pooled$metacell_meta <- meta
  pooled$input_sample_col <- sample_col
  pooled$analysis_sample_col <- internal_pool_col
  pooled$condition_col <- condition_col
  pooled$celltype_col <- celltype_col
  pooled$pooling_scope <- "condition_by_cell_type"
  pooled$sample_weighting <- "none"
  pooled$cache_contract <- cache_contract
  pooled$input_design <- list(
    metacell_grouping = c(condition_col, celltype_col),
    condition_celltype_stratification = TRUE,
    condition_only_stratification = FALSE,
    supercell_label_col = NULL,
    celltype_assignment = "hard condition-by-cell-type stratum",
    ambiguous_celltype_policy = "not_applicable_strata_are_pure",
    gamma = metacell_args$gamma,
    cache_contract_schema = cache_contract$schema_version,
    inference_policy = paste(
      "cells are stratified by condition and broad cell type; sample metadata",
      "are retained for composition diagnostics but are not used for",
      "selection, weighting, or metacell grouping"
    )
  )
  pooled
}
