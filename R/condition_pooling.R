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
      " dimensions but dimension ", max(dims), " was requested.",
      call. = FALSE
    )
  }
  index <- match(cells, rownames(embeddings))
  if (anyNA(index)) {
    stop("Reduction `", reduction, "` lacks required input cells.",
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
    rna_reduction = "pca",
    atac_reduction = "lsi",
    rna_dims = 1:30,
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 20L,
    min_metacell_size = 1L,
    min_metacells_per_stratum = 1L,
    k.knn = 5L,
    do.approx = FALSE,
    approx.N = 20000L,
    block.size = 10000L,
    igraph.clustering = "walktrap"
  )
  analysis_args <- modifyList(defaults, metacell_args[
    intersect(names(metacell_args), names(defaults))
  ])
  meta_signature <- data.frame(
    cell_id = cells,
    condition = as.character(object@meta.data[[condition_col]][meta_index]),
    cell_type = as.character(object@meta.data[[celltype_col]][meta_index]),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = "regcompass_native_supercell_metacell_cache_v1",
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    native_supercell_api = "SCimplify_from_embedding",
    condition_argument = "cell.split.condition",
    celltype_argument = "cell.annotation",
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
  all(file.exists(c(
    file.path(outdir, "metacell_metadata.tsv.gz"),
    file.path(outdir, "membership.tsv.gz"),
    file.path(outdir, "rna_counts.rds"),
    file.path(outdir, "atac_counts.rds"),
    file.path(outdir, "metacell_object.rds")
  )))
}

.rc_validate_condition_metacell_cache <- function(
    outdir, contract, overwrite = FALSE) {
  if (!.rc_condition_metacell_has_checkpoints(outdir) || isTRUE(overwrite)) {
    return(invisible(FALSE))
  }
  contract_file <- file.path(outdir, "condition_metacell_cache_contract.rds")
  if (!file.exists(contract_file) ||
      !identical(readRDS(contract_file), contract)) {
    stop(
      "Existing metacell checkpoints use a different native SuperCell contract; set `metacell_args$overwrite = TRUE`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_scale_embedding_block <- function(x) {
  x <- as.matrix(x)
  x <- scale(x)
  x[!is.finite(x)] <- 0
  if (ncol(x)) x / sqrt(ncol(x)) else x
}

.rc_native_supercell_membership <- function(
    object, condition_col, celltype_col, rna_reduction, atac_reduction,
    rna_dims, atac_dims, gamma, seed, k.knn, do.approx, approx.N,
    block.size, igraph.clustering) {
  cells <- colnames(object)
  rna <- SeuratObject::Embeddings(object[[rna_reduction]])[
    cells, as.integer(rna_dims), drop = FALSE
  ]
  atac <- SeuratObject::Embeddings(object[[atac_reduction]])[
    cells, as.integer(atac_dims), drop = FALSE
  ]
  embedding <- cbind(
    .rc_scale_embedding_block(rna),
    .rc_scale_embedding_block(atac)
  )
  meta <- object@meta.data[cells, , drop = FALSE]
  condition <- as.character(meta[[condition_col]])
  condition_input <- if (grepl("^\\.regcompass_condition", condition_col) &&
      length(unique(condition)) == 1L) {
    NULL
  } else {
    condition
  }
  result <- SuperCell::SCimplify_from_embedding(
    X = embedding,
    cell.annotation = as.character(meta[[celltype_col]]),
    cell.split.condition = condition_input,
    gamma = as.integer(gamma),
    k.knn = as.integer(k.knn),
    n.pc = ncol(embedding),
    do.approx = isTRUE(do.approx),
    approx.N = as.integer(approx.N),
    block.size = as.integer(block.size),
    seed = as.integer(seed),
    igraph.clustering = igraph.clustering,
    return.singlecell.NW = FALSE,
    return.hierarchical.structure = TRUE
  )
  raw <- as.character(result$membership)
  if (length(raw) != length(cells)) {
    stop("SuperCell membership length differs from the input cell count.",
         call. = FALSE)
  }
  if (!is.null(names(result$membership))) {
    index <- match(cells, names(result$membership))
    if (!anyNA(index)) raw <- raw[index]
  }
  levels <- unique(raw)
  width <- max(3L, nchar(length(levels)))
  ids <- paste0("MC", sprintf(paste0("%0", width, "d"), seq_along(levels)))
  map <- stats::setNames(ids, levels)
  membership <- data.frame(
    cell_id = cells,
    metacell_id = unname(map[raw]),
    stringsAsFactors = FALSE
  )
  list(membership = membership, supercell = result)
}

.rc_aggregate_native_metacell_counts <- function(
    object, membership, rna_assay, atac_assay) {
  cells <- colnames(object)
  membership <- membership[match(cells, membership$cell_id), , drop = FALSE]
  groups <- unique(membership$metacell_id)
  design <- Matrix::sparseMatrix(
    i = seq_along(cells),
    j = match(membership$metacell_id, groups),
    x = 1,
    dims = c(length(cells), length(groups)),
    dimnames = list(cells, groups)
  )
  rna_counts <- .rc_as_sparse(
    .rc_get_assay_counts(object, rna_assay)[, cells, drop = FALSE] %*% design
  )
  atac_counts <- .rc_as_sparse(
    .rc_get_assay_counts(object, atac_assay)[, cells, drop = FALSE] %*% design
  )
  mc <- SeuratObject::CreateSeuratObject(
    counts = rna_counts,
    assay = rna_assay,
    project = "RegCompassNativeSuperCell"
  )
  mc[[atac_assay]] <- SeuratObject::CreateAssayObject(counts = atac_counts)
  list(object = mc, rna_counts = rna_counts, atac_counts = atac_counts)
}

.rc_make_condition_celltype_metacells <- function(
    object, outdir,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list()) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!is.list(metacell_args)) stop("`metacell_args` must be a list.", call. = FALSE)
  if (!identical(fragment_files, FALSE) && !is.null(fragment_files)) {
    stop("Native SuperCell metacells aggregate the existing ATAC count assay; `fragment_files` must be FALSE.", call. = FALSE)
  }
  .rc_validate_condition_celltype_metadata(
    object@meta.data, condition_col, celltype_col
  )
  if (!is.null(cell_type)) {
    requested <- unique(trimws(as.character(cell_type)))
    missing <- setdiff(requested, unique(as.character(object@meta.data[[celltype_col]])))
    if (length(missing)) {
      stop("Requested cell types were not found: ", paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    cells <- rownames(object@meta.data)[
      as.character(object@meta.data[[celltype_col]]) %in% requested
    ]
    object <- subset(object, cells = cells)
  }
  reserved <- intersect(names(metacell_args), c(
    "object", "outdir", "condition_col", "celltype_col", "cell_type",
    "rna_assay", "atac_assay", "fragment_files", "cell.annotation",
    "cell.split.condition"
  ))
  if (length(reserved)) {
    stop("`metacell_args` cannot override managed fields: ",
         paste(reserved, collapse = ", "), call. = FALSE)
  }
  defaults <- list(
    rna_reduction = "pca",
    atac_reduction = "lsi",
    rna_dims = 1:30,
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 20L,
    min_metacell_size = 1L,
    min_metacells_per_stratum = 1L,
    k.knn = 5L,
    do.approx = FALSE,
    approx.N = 20000L,
    block.size = 10000L,
    igraph.clustering = "walktrap",
    overwrite = FALSE
  )
  args <- modifyList(defaults, metacell_args)
  numeric_controls <- c(
    "gamma", "seed", "min_cells_per_stratum", "min_metacell_size",
    "min_metacells_per_stratum", "k.knn", "approx.N", "block.size"
  )
  for (field in numeric_controls) args[[field]] <- as.integer(args[[field]])
  if (any(vapply(args[numeric_controls], function(x) {
    length(x) != 1L || is.na(x) || x < 1L
  }, logical(1)))) {
    stop("SuperCell numeric controls must be positive integers.", call. = FALSE)
  }
  cells <- colnames(object)
  meta <- object@meta.data[cells, , drop = FALSE]
  stratum <- interaction(
    meta[, c(condition_col, celltype_col), drop = FALSE],
    drop = TRUE, lex.order = TRUE
  )
  stratum_size <- table(stratum)
  if (any(stratum_size < args$min_cells_per_stratum)) {
    stop(
      "Condition/cell-type strata below `min_cells_per_stratum`: ",
      paste(names(stratum_size)[stratum_size < args$min_cells_per_stratum], collapse = ", "),
      call. = FALSE
    )
  }
  contract <- .rc_condition_metacell_cache_contract(
    object, condition_col, celltype_col, rna_assay, atac_assay, args
  )
  .rc_validate_condition_metacell_cache(
    outdir, contract, overwrite = isTRUE(args$overwrite)
  )
  if (.rc_condition_metacell_has_checkpoints(outdir) && !isTRUE(args$overwrite)) {
    mc <- readRDS(file.path(outdir, "metacell_object.rds"))
    membership <- utils::read.delim(
      gzfile(file.path(outdir, "membership.tsv.gz")),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    mc_meta <- utils::read.delim(
      gzfile(file.path(outdir, "metacell_metadata.tsv.gz")),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    aggregated <- list(
      object = mc,
      rna_counts = readRDS(file.path(outdir, "rna_counts.rds")),
      atac_counts = readRDS(file.path(outdir, "atac_counts.rds"))
    )
    native <- list(supercell = NULL)
  } else {
    native <- .rc_native_supercell_membership(
      object = object,
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_reduction = args$rna_reduction,
      atac_reduction = args$atac_reduction,
      rna_dims = args$rna_dims,
      atac_dims = args$atac_dims,
      gamma = args$gamma,
      seed = args$seed,
      k.knn = args$k.knn,
      do.approx = args$do.approx,
      approx.N = args$approx.N,
      block.size = args$block.size,
      igraph.clustering = args$igraph.clustering
    )
    membership <- native$membership
    source_index <- match(membership$cell_id, rownames(object@meta.data))
    membership[[condition_col]] <- as.character(
      object@meta.data[[condition_col]][source_index]
    )
    membership[[celltype_col]] <- as.character(
      object@meta.data[[celltype_col]][source_index]
    )
    mc_meta <- rc_build_metacell_metadata(membership)
    groups <- split(seq_len(nrow(membership)), membership$metacell_id)
    impure <- names(groups)[vapply(groups, function(rows) {
      length(unique(membership[[condition_col]][rows])) != 1L ||
        length(unique(membership[[celltype_col]][rows])) != 1L
    }, logical(1))]
    if (length(impure)) {
      stop("Native SuperCell returned mixed condition/cell-type metacells: ",
           paste(utils::head(impure, 10L), collapse = ", "), call. = FALSE)
    }
    stratum_mc <- table(interaction(
      mc_meta[, c(condition_col, celltype_col), drop = FALSE],
      drop = TRUE, lex.order = TRUE
    ))
    if (any(stratum_mc < args$min_metacells_per_stratum)) {
      stop("Native SuperCell produced too few metacells in strata: ",
           paste(names(stratum_mc)[stratum_mc < args$min_metacells_per_stratum], collapse = ", "), call. = FALSE)
    }
    mc_meta$low_power_metacell <- mc_meta$n_cells < args$min_metacell_size
    mc_meta$effective_gamma <- args$gamma
    mc_meta$requested_gamma <- args$gamma
    mc_meta$fixed_gamma <- TRUE
    mc_meta$pooling_scope <- "native_supercell_condition_and_celltype"
    mc_meta$celltype_role <- "SuperCell_cell.annotation"
    aggregated <- .rc_aggregate_native_metacell_counts(
      object, membership, rna_assay, atac_assay
    )
    rownames(mc_meta) <- mc_meta$metacell_id
    aggregated$object@meta.data <- mc_meta[
      colnames(aggregated$object), , drop = FALSE
    ]
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    .rc_write_tsv_gz(membership, file.path(outdir, "membership.tsv.gz"))
    .rc_write_tsv_gz(mc_meta, file.path(outdir, "metacell_metadata.tsv.gz"))
    saveRDS(aggregated$rna_counts, file.path(outdir, "rna_counts.rds"))
    saveRDS(aggregated$atac_counts, file.path(outdir, "atac_counts.rds"))
    saveRDS(aggregated$object, file.path(outdir, "metacell_object.rds"))
    saveRDS(contract, file.path(outdir, "condition_metacell_cache_contract.rds"))
  }
  celltype_composition <- data.frame(
    metacell_id = as.character(mc_meta$metacell_id),
    value = as.character(mc_meta[[celltype_col]]),
    n_cells = as.integer(mc_meta$n_cells),
    stringsAsFactors = FALSE
  )
  colnames(celltype_composition)[[2L]] <- celltype_col
  celltype_summary <- data.frame(
    metacell_id = as.character(mc_meta$metacell_id),
    dominant_celltype = as.character(mc_meta[[celltype_col]]),
    dominant_celltype_fraction = 1,
    n_celltypes = 1L,
    mixed_celltype_metacell = FALSE,
    dominant_celltype_tied = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    metacell_objects = list(native_supercell = aggregated$object),
    metacell_object = aggregated$object,
    metacell_meta = mc_meta,
    membership = membership,
    rna_counts = aggregated$rna_counts,
    atac_counts = aggregated$atac_counts,
    celltype_composition = celltype_composition,
    celltype_composition_summary = celltype_summary,
    condition_col = condition_col,
    celltype_col = celltype_col,
    selected_cell_types = unique(as.character(mc_meta[[celltype_col]])),
    pooling_scope = "native_supercell_condition_and_celltype",
    cache_contract = contract,
    input_design = list(
      metacell_grouping = c(condition_col, celltype_col),
      native_supercell_api = "SCimplify_from_embedding",
      condition_argument = "cell.split.condition",
      celltype_argument = "cell.annotation",
      temporary_combined_stratum = FALSE,
      gamma = args$gamma,
      inference_policy = paste(
        "SuperCell receives condition and cell type as separate native inputs;",
        "no concatenated condition__cell_type field is created"
      ),
      sample_metadata = "not_used_or_retained"
    )
  )
}
