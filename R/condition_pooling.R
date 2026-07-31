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

.rc_condition_metacell_defaults <- function() {
  list(
    rna_reduction = "pca",
    atac_reduction = "lsi",
    rna_dims = 1:30,
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 20L,
    min_metacell_size = 1L,
    min_metacells_per_stratum = 1L,
    k.knn = 30L,
    kith = NULL,
    kernel = TRUE,
    graph.name = NULL,
    metacellNormalization = FALSE,
    avg.in.data = FALSE,
    verbose = FALSE,
    overwrite = FALSE
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
  defaults <- .rc_condition_metacell_defaults()
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
    schema_version = "regcompass_celltype_wnn_condition_joint_cache",
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    native_supercell_api = "SCimplify_by_graph_group",
    graph_group_argument = "cell.graph.group",
    condition_argument = "cell.split.condition",
    graph_scope = "one_independent_WNN_graph_per_cell_type",
    condition_scope = "all_conditions_joint_within_cell_type_graph",
    membership_split_timing = "after_joint_WNN_graph_clustering",
    graph_method = "SuperCell_multimodal_WNN_then_walktrap",
    modality_weighting = "adaptive_WNN_within_cell_type",
    aggregation_method = "SCimplify_for_Seurat_membership_mode",
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
      "Existing metacell checkpoints use a different grouped-WNN contract; set `metacell_args$overwrite = TRUE`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_require_supercell_api <- function() {
  if (!requireNamespace("SuperCell", quietly = TRUE)) {
    stop(
      "Package 'SuperCell' is required. Install the current default branch of 1667857557/SuperCell_Seurat_V4.",
      call. = FALSE
    )
  }
  exports <- getNamespaceExports("SuperCell")
  required_api <- c("SCimplify_by_graph_group", "SCimplify_for_Seurat")
  missing_api <- setdiff(required_api, exports)
  if (length(missing_api)) {
    stop("Installed SuperCell lacks required API(s): ",
         paste(missing_api, collapse = ", "), ".", call. = FALSE)
  }
  grouped <- getExportedValue("SuperCell", "SCimplify_by_graph_group")
  aggregate <- getExportedValue("SuperCell", "SCimplify_for_Seurat")
  grouped_required <- c(
    "seurat", "cell.graph.group", "cell.split.condition", "k.knn", "kith",
    "kernel", "gamma", "graph.name", "assay", "reduction", "dims", "seed"
  )
  aggregate_required <- c(
    "seurat", "assay", "reduction", "dims", "membership", "return.seurat"
  )
  missing_grouped <- setdiff(grouped_required, names(formals(grouped)))
  missing_aggregate <- setdiff(aggregate_required, names(formals(aggregate)))
  if (length(missing_grouped) || length(missing_aggregate)) {
    stop(
      "Installed SuperCell API formals are incompatible. Missing grouped formals: ",
      paste(missing_grouped, collapse = ", "),
      "; missing aggregation formals: ",
      paste(missing_aggregate, collapse = ", "), ".",
      call. = FALSE
    )
  }
  list(grouped = grouped, aggregate = aggregate)
}

.rc_call_supercell <- function(fun, args) {
  supported <- intersect(names(args), names(formals(fun)))
  do.call(fun, args[supported])
}

.rc_build_grouped_wnn_membership <- function(
    object, condition_col, celltype_col, rna_assay, atac_assay,
    rna_reduction, atac_reduction, rna_dims, atac_dims,
    gamma, seed, k.knn, kith, kernel, graph.name, verbose) {
  api <- .rc_require_supercell_api()
  cells <- as.character(colnames(object))
  meta <- object@meta.data[cells, , drop = FALSE]
  result <- .rc_call_supercell(api$grouped, list(
    seurat = object,
    cell.graph.group = stats::setNames(
      as.character(meta[[celltype_col]]), cells
    ),
    cell.split.condition = stats::setNames(
      as.character(meta[[condition_col]]), cells
    ),
    k.knn = as.integer(k.knn),
    kith = kith,
    kernel = isTRUE(kernel),
    gamma = as.numeric(gamma),
    graph.name = graph.name,
    assay = c(rna_assay, atac_assay),
    reduction = list(rna_reduction, atac_reduction),
    dims = list(as.integer(rna_dims), as.integer(atac_dims)),
    seed = as.integer(seed),
    return.group.results = FALSE,
    verbose = isTRUE(verbose)
  ))
  membership <- result$membership_table
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id", "parent_metacell_id") %in%
           colnames(membership))) {
    stop("SCimplify_by_graph_group() returned an incompatible membership table.",
         call. = FALSE)
  }
  index <- match(cells, membership$cell_id)
  if (anyNA(index)) {
    stop("Grouped SuperCell membership does not cover every input cell.",
         call. = FALSE)
  }
  membership <- membership[index, , drop = FALSE]
  if (anyDuplicated(membership$cell_id) ||
      anyNA(membership$metacell_id) || any(!nzchar(membership$metacell_id))) {
    stop("Grouped SuperCell membership IDs are invalid.", call. = FALSE)
  }
  list(
    membership = membership[, c(
      "cell_id", "metacell_id", "parent_metacell_id"
    ), drop = FALSE],
    parent_hierarchies = result$h_membership %||% list(),
    upstream_api = "SCimplify_by_graph_group"
  )
}

.rc_aggregate_metacell_counts <- function(
    object, membership, rna_assay, atac_assay,
    rna_reduction, atac_reduction, rna_dims, atac_dims,
    seed, metacellNormalization, avg.in.data, verbose,
    parent_hierarchies = list()) {
  api <- .rc_require_supercell_api()
  cells <- as.character(colnames(object))
  index <- match(cells, membership$cell_id)
  if (anyNA(index)) {
    stop("Metacell membership does not cover every input cell.", call. = FALSE)
  }
  membership <- membership[index, , drop = FALSE]
  groups <- unique(as.character(membership$metacell_id))
  numeric_membership <- match(membership$metacell_id, groups)
  names(numeric_membership) <- cells
  set.seed(as.integer(seed))
  mc <- .rc_call_supercell(api$aggregate, list(
    seurat = object,
    assay = c(rna_assay, atac_assay),
    reduction = list(rna_reduction, atac_reduction),
    dims = list(as.integer(rna_dims), as.integer(atac_dims)),
    membership = numeric_membership,
    metacellNormalization = isTRUE(metacellNormalization),
    avg.in.data = isTRUE(avg.in.data),
    fragmentFiles = NULL,
    return.seurat = TRUE,
    verbose = isTRUE(verbose)
  ))
  if (!inherits(mc, "Seurat") || ncol(mc) != length(groups)) {
    stop("SuperCell membership aggregation returned an incompatible object.",
         call. = FALSE)
  }
  current <- colnames(mc)
  parsed <- suppressWarnings(as.integer(sub("^.*_", "", current)))
  if (length(parsed) == length(groups) &&
      !anyNA(parsed) && setequal(parsed, seq_along(groups))) {
    mc <- mc[, order(parsed)]
  }
  if (ncol(mc) != length(groups)) {
    stop("SuperCell aggregate columns cannot be aligned to memberships.",
         call. = FALSE)
  }
  mc <- SeuratObject::RenameCells(mc, new.names = groups)
  rna_counts <- .rc_as_sparse(.rc_get_assay_counts(mc, rna_assay))
  atac_counts <- .rc_as_sparse(.rc_get_assay_counts(mc, atac_assay))
  mc@misc$membership_table <- membership
  mc@misc$regcompass_supercell_parent_hierarchies <- parent_hierarchies
  mc@misc$regcompass_supercell_aggregation <- list(
    api = "SCimplify_for_Seurat",
    mode = "provided_condition_pure_membership",
    metacellNormalization = isTRUE(metacellNormalization),
    avg.in.data = isTRUE(avg.in.data)
  )
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
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
  }
  if (!identical(fragment_files, FALSE) && !is.null(fragment_files)) {
    stop(
      "This condition-joint path aggregates the existing ATAC count assay; `fragment_files` must be FALSE.",
      call. = FALSE
    )
  }
  .rc_validate_condition_celltype_metadata(
    object@meta.data, condition_col, celltype_col
  )
  if (!is.null(cell_type)) {
    requested <- unique(trimws(as.character(cell_type)))
    available <- unique(as.character(object@meta.data[[celltype_col]]))
    missing <- setdiff(requested, available)
    if (length(missing)) {
      stop("Requested cell types were not found: ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    cells <- rownames(object@meta.data)[
      as.character(object@meta.data[[celltype_col]]) %in% requested
    ]
    object <- subset(object, cells = cells)
  }
  reserved <- intersect(names(metacell_args), c(
    "object", "outdir", "condition_col", "celltype_col", "cell_type",
    "rna_assay", "atac_assay", "fragment_files", "membership", "assay",
    "reduction", "dims", "label", "return.seurat",
    "cell.graph.group", "cell.split.condition"
  ))
  if (length(reserved)) {
    stop("`metacell_args` cannot override managed fields: ",
         paste(reserved, collapse = ", "), call. = FALSE)
  }
  retired <- intersect(names(metacell_args), c(
    "do.approx", "approx.N", "block.size", "igraph.clustering",
    "embedding_scaling"
  ))
  if (length(retired)) {
    stop(
      "Embedding-only graph controls are no longer accepted: ",
      paste(retired, collapse = ", "),
      ". Use grouped multimodal WNN controls.",
      call. = FALSE
    )
  }
  args <- modifyList(.rc_condition_metacell_defaults(), metacell_args)
  integer_controls <- c(
    "gamma", "seed", "min_cells_per_stratum", "min_metacell_size",
    "min_metacells_per_stratum", "k.knn"
  )
  for (field in integer_controls) args[[field]] <- as.integer(args[[field]])
  if (any(vapply(args[integer_controls], function(x) {
    length(x) != 1L || is.na(x) || x < 1L
  }, logical(1)))) {
    stop("SuperCell numeric controls must be positive integers.", call. = FALSE)
  }
  if (!is.null(args$kith)) {
    args$kith <- as.integer(args$kith)
    if (length(args$kith) != 1L || is.na(args$kith) || args$kith < 1L) {
      stop("`kith` must be NULL or one positive integer.", call. = FALSE)
    }
  }
  for (field in c("kernel", "metacellNormalization", "avg.in.data",
                  "verbose", "overwrite")) {
    if (!is.logical(args[[field]]) || length(args[[field]]) != 1L ||
        is.na(args[[field]])) {
      stop("`", field, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  cells <- as.character(colnames(object))
  meta <- object@meta.data[cells, , drop = FALSE]
  stratum <- interaction(
    meta[, c(condition_col, celltype_col), drop = FALSE],
    drop = TRUE, lex.order = TRUE
  )
  stratum_size <- table(stratum)
  if (any(stratum_size < args$min_cells_per_stratum)) {
    stop(
      "Condition/cell-type strata below `min_cells_per_stratum`: ",
      paste(names(stratum_size)[stratum_size < args$min_cells_per_stratum],
            collapse = ", "),
      call. = FALSE
    )
  }
  contract <- .rc_condition_metacell_cache_contract(
    object, condition_col, celltype_col, rna_assay, atac_assay, args
  )
  .rc_validate_condition_metacell_cache(
    outdir, contract, overwrite = isTRUE(args$overwrite)
  )
  if (.rc_condition_metacell_has_checkpoints(outdir) &&
      !isTRUE(args$overwrite)) {
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
    grouped <- list(parent_hierarchies = list())
  } else {
    grouped <- .rc_build_grouped_wnn_membership(
      object = object,
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      rna_reduction = args$rna_reduction,
      atac_reduction = args$atac_reduction,
      rna_dims = args$rna_dims,
      atac_dims = args$atac_dims,
      gamma = args$gamma,
      seed = args$seed,
      k.knn = args$k.knn,
      kith = args$kith,
      kernel = args$kernel,
      graph.name = args$graph.name,
      verbose = args$verbose
    )
    membership <- grouped$membership
    source_index <- match(membership$cell_id, rownames(object@meta.data))
    if (anyNA(source_index)) {
      stop("Grouped membership cannot be aligned to Seurat metadata.",
           call. = FALSE)
    }
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
      stop("Grouped WNN produced impure metacells: ",
           paste(utils::head(impure, 10L), collapse = ", "), call. = FALSE)
    }
    stratum_mc <- table(interaction(
      mc_meta[, c(condition_col, celltype_col), drop = FALSE],
      drop = TRUE, lex.order = TRUE
    ))
    if (any(stratum_mc < args$min_metacells_per_stratum)) {
      stop("SuperCell produced too few metacells in strata: ",
           paste(names(stratum_mc)[stratum_mc < args$min_metacells_per_stratum],
                 collapse = ", "), call. = FALSE)
    }
    mc_meta$low_power_metacell <- mc_meta$n_cells < args$min_metacell_size
    mc_meta$requested_gamma <- args$gamma
    mc_meta$pooling_scope <- "celltype_grouped_joint_condition_WNN"
    mc_meta$celltype_role <- "one_independent_WNN_graph_per_cell_type"
    mc_meta$condition_role <- "post_WNN_clustering_membership_split"
    aggregated <- .rc_aggregate_metacell_counts(
      object = object,
      membership = membership,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      rna_reduction = args$rna_reduction,
      atac_reduction = args$atac_reduction,
      rna_dims = args$rna_dims,
      atac_dims = args$atac_dims,
      seed = args$seed,
      metacellNormalization = args$metacellNormalization,
      avg.in.data = args$avg.in.data,
      verbose = args$verbose,
      parent_hierarchies = grouped$parent_hierarchies
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
    metacell_objects = list(grouped_wnn = aggregated$object),
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
    pooling_scope = "celltype_grouped_joint_condition_WNN",
    cache_contract = contract,
    input_design = list(
      metacell_purity_grouping = c(condition_col, celltype_col),
      graph_grouping = celltype_col,
      condition_pooling = condition_col,
      native_supercell_api = "SCimplify_by_graph_group",
      graph_group_argument = "cell.graph.group",
      condition_argument = "cell.split.condition",
      graph_method = "multimodal_WNN",
      clustering_method = "walktrap_cut_at",
      aggregation_method = "SCimplify_for_Seurat_with_membership",
      graph_scope = "one_independent_WNN_graph_per_cell_type",
      condition_scope = "all_conditions_joint_within_cell_type_graph",
      membership_split_timing = "after_joint_WNN_graph_clustering",
      modality_weighting = "adaptive_WNN_within_cell_type",
      temporary_combined_stratum = FALSE,
      gamma = args$gamma,
      k.knn = args$k.knn,
      kernel = args$kernel,
      inference_policy = paste(
        "Each broad cell type receives one independent multimodal WNN graph;",
        "all conditions jointly determine WNN neighbours and Walktrap clusters;",
        "condition splits membership only after clustering"
      ),
      sample_metadata = "not_used_or_retained"
    )
  )
}
