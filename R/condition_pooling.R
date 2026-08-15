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
    stop("Metacell fingerprinting requires a two-dimensional matrix.",
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
    verbose = FALSE
  )
}

.rc_condition_metacell_cache_contract <- function(
    object, condition_col, celltype_col, rna_assay, atac_assay,
    metacell_args) {
  cells <- as.character(colnames(object))
  meta_index <- match(cells, rownames(object@meta.data))
  if (anyNA(meta_index)) {
    stop("Seurat cell order cannot be aligned to metadata for fingerprinting.",
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
    schema_version = "regcompass_shared_walktrap_condition_cut_cache_v1",
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    native_supercell_api = "SCimplify_by_graph_group",
    graph_group_argument = "cell.graph.group",
    condition_argument = "cell.split.condition",
    condition_partition = "hierarchy_constrained",
    partition_schema_version = "shared_walktrap_condition_cut_v1",
    graph_scope = "one_independent_WNN_graph_per_cell_type",
    condition_scope =
      "shared_WNN_and_Walktrap_with_condition_specific_hierarchy_cut",
    membership_split_timing = "condition_specific_cut_of_shared_walktrap_hierarchy",
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

.rc_require_supercell_api <- function() {
  if (!requireNamespace("SuperCell", quietly = TRUE)) {
    stop(
      "Package 'SuperCell' is required. Install 1667857557/SuperCell_Seurat_V4 >= 2.2.0.",
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
    "kernel", "gamma", "condition.partition", "min.metacell.size",
    "min.metacells.per.condition", "graph.name", "assay", "reduction",
    "dims", "seed"
  )
  aggregate_required <- c(
    "seurat", "assay", "reduction", "dims", "membership", "return.seurat"
  )
  missing_grouped <- setdiff(grouped_required, names(formals(grouped)))
  missing_aggregate <- setdiff(aggregate_required, names(formals(aggregate)))
  if (length(missing_grouped) || length(missing_aggregate)) {
    stop(
      "Installed SuperCell API is incompatible with RegCompass hierarchy-constrained Stage 2. Missing grouped formals: ",
      paste(missing_grouped, collapse = ", "),
      "; missing aggregation formals: ",
      paste(missing_aggregate, collapse = ", "), ".",
      call. = FALSE
    )
  }
  list(grouped = grouped, aggregate = aggregate)
}

.rc_call_supercell <- function(fun, args) {
  formal_names <- names(formals(fun))
  unsupported <- setdiff(names(args), formal_names)
  if (length(unsupported) && !"..." %in% formal_names) {
    stop(
      "Installed SuperCell cannot accept required argument(s): ",
      paste(unsupported, collapse = ", "), ".", call. = FALSE
    )
  }
  do.call(fun, args)
}

.rc_build_grouped_wnn_membership <- function(
    object, condition_col, celltype_col, rna_assay, atac_assay,
    rna_reduction, atac_reduction, rna_dims, atac_dims,
    gamma, seed, k.knn, kith, kernel, graph.name, verbose,
    min_metacell_size, min_metacells_per_stratum) {
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
    condition.partition = "hierarchy_constrained",
    min.metacell.size = as.integer(min_metacell_size),
    min.metacells.per.condition = as.integer(min_metacells_per_stratum),
    graph.name = graph.name,
    assay = c(rna_assay, atac_assay),
    reduction = list(rna_reduction, atac_reduction),
    dims = list(as.integer(rna_dims), as.integer(atac_dims)),
    seed = as.integer(seed),
    return.group.results = FALSE,
    verbose = isTRUE(verbose)
  ))
  membership <- result$membership_table
  required_membership <- c(
    "cell_id", "metacell_id", "parent_metacell_id", "graph_group",
    "condition", "partition_policy", "shared_cut_k", "shared_community_id"
  )
  if (!is.data.frame(membership) ||
      !all(required_membership %in% colnames(membership))) {
    stop(
      "SCimplify_by_graph_group() returned an incompatible hierarchy membership table.",
      call. = FALSE
    )
  }
  if (!identical(as.character(result$partition_policy), "hierarchy_constrained") ||
      !identical(as.character(result$partition_schema_version),
                 "shared_walktrap_condition_cut_v1")) {
    stop("SuperCell did not execute the required hierarchy-constrained policy.",
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
  if (any(as.character(membership$graph_group) !=
          as.character(meta[[celltype_col]])) ||
      any(as.character(membership$condition) !=
          as.character(meta[[condition_col]])) ||
      any(as.character(membership$partition_policy) !=
          "hierarchy_constrained")) {
    stop(
      "SuperCell hierarchy membership provenance disagrees with RegCompass cell metadata.",
      call. = FALSE
    )
  }
  list(
    membership = membership,
    parent_hierarchies = result$h_membership %||% list(),
    partition_diagnostics = result$partition_diagnostics %||% data.frame(),
    partition_policy = result$partition_policy,
    partition_schema_version = result$partition_schema_version,
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

.rc_positive_integer_control <- function(value, name) {
  numeric <- suppressWarnings(as.numeric(value))
  if (length(numeric) != 1L || !is.finite(numeric) || numeric < 1 ||
      abs(numeric - round(numeric)) > sqrt(.Machine$double.eps)) {
    stop("`", name, "` must be one positive integer.", call. = FALSE)
  }
  as.integer(round(numeric))
}

.rc_make_condition_celltype_metacells <- function(
    object, outdir,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    metacell_args = list()) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
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
    "rna_assay", "atac_assay", "membership", "assay",
    "reduction", "dims", "label", "return.seurat",
    "cell.graph.group", "cell.split.condition", "condition.partition",
    "min.metacell.size", "min.metacells.per.condition"
  ))
  if (length(reserved)) {
    stop("`metacell_args` cannot override managed fields: ",
         paste(reserved, collapse = ", "), call. = FALSE)
  }
  retired <- intersect(names(metacell_args), c(
    "do.approx", "approx.N", "block.size", "igraph.clustering",
    "embedding_scaling", "overwrite", "return.graph", "repair_mode",
    "repair_max_iter", "affinity_repair"
  ))
  if (length(retired)) {
    stop(
      "Retired Stage 2 controls are no longer accepted: ",
      paste(retired, collapse = ", "),
      ". Stage 2 uses one shared Walktrap hierarchy per cell type and no post-hoc graph repair.",
      call. = FALSE
    )
  }
  args <- modifyList(.rc_condition_metacell_defaults(), metacell_args)
  integer_controls <- c(
    "gamma", "seed", "min_cells_per_stratum", "min_metacell_size",
    "min_metacells_per_stratum", "k.knn"
  )
  for (field in integer_controls) {
    args[[field]] <- .rc_positive_integer_control(args[[field]], field)
  }
  if (!is.null(args$kith)) {
    args$kith <- .rc_positive_integer_control(args$kith, "kith")
  }
  for (field in c("kernel", "metacellNormalization", "avg.in.data",
                  "verbose")) {
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
  feasibility_floor <- args$min_metacell_size * args$min_metacells_per_stratum
  if (any(stratum_size < feasibility_floor)) {
    stop(
      "Condition/cell-type strata are mathematically infeasible for the requested hard metacell constraints (N < min_metacell_size * min_metacells_per_stratum): ",
      paste(names(stratum_size)[stratum_size < feasibility_floor],
            collapse = ", "), call. = FALSE
    )
  }
  contract <- .rc_condition_metacell_cache_contract(
    object, condition_col, celltype_col, rna_assay, atac_assay, args
  )

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
    verbose = args$verbose,
    min_metacell_size = args$min_metacell_size,
    min_metacells_per_stratum = args$min_metacells_per_stratum
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
    stop("Hierarchy-constrained SuperCell produced impure metacells: ",
         paste(utils::head(impure, 10L), collapse = ", "), call. = FALSE)
  }
  if (any(mc_meta$n_cells < args$min_metacell_size)) {
    bad <- mc_meta$metacell_id[mc_meta$n_cells < args$min_metacell_size]
    stop(
      "Hierarchy-constrained SuperCell violated the hard minimum metacell size: ",
      paste(utils::head(bad, 10L), collapse = ", "), call. = FALSE
    )
  }
  stratum_mc <- table(interaction(
    mc_meta[, c(condition_col, celltype_col), drop = FALSE],
    drop = TRUE, lex.order = TRUE
  ))
  if (any(stratum_mc < args$min_metacells_per_stratum)) {
    stop("Hierarchy-constrained SuperCell produced too few metacells in strata: ",
         paste(names(stratum_mc)[stratum_mc < args$min_metacells_per_stratum],
               collapse = ", "), call. = FALSE)
  }
  mc_meta$requested_gamma <- args$gamma
  mc_meta$min_metacell_size <- args$min_metacell_size
  mc_meta$min_metacells_per_stratum <- args$min_metacells_per_stratum
  mc_meta$pooling_scope <- "celltype_shared_WNN_condition_hierarchy_cut"
  mc_meta$celltype_role <- "one_independent_WNN_and_Walktrap_per_cell_type"
  mc_meta$condition_role <- "condition_specific_cut_of_shared_walktrap_hierarchy"
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
  aggregated$object@misc$regcompass_supercell_partition <- list(
    policy = grouped$partition_policy,
    schema_version = grouped$partition_schema_version,
    diagnostics = grouped$partition_diagnostics
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(membership, file.path(outdir, "membership.tsv.gz"))
  .rc_write_tsv_gz(mc_meta, file.path(outdir, "metacell_metadata.tsv.gz"))
  if (is.data.frame(grouped$partition_diagnostics) &&
      nrow(grouped$partition_diagnostics)) {
    .rc_write_tsv_gz(
      grouped$partition_diagnostics,
      file.path(outdir, "metacell_partition_diagnostics.tsv.gz")
    )
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
    metacell_object = aggregated$object,
    metacell_meta = mc_meta,
    membership = membership,
    partition_diagnostics = grouped$partition_diagnostics,
    partition_policy = grouped$partition_policy,
    partition_schema_version = grouped$partition_schema_version,
    celltype_composition = celltype_composition,
    celltype_composition_summary = celltype_summary,
    condition_col = condition_col,
    celltype_col = celltype_col,
    selected_cell_types = unique(as.character(mc_meta[[celltype_col]])),
    pooling_scope = "celltype_shared_WNN_condition_hierarchy_cut",
    cache_contract = contract,
    input_design = list(
      metacell_purity_grouping = c(condition_col, celltype_col),
      graph_grouping = celltype_col,
      condition_pooling = condition_col,
      native_supercell_api = "SCimplify_by_graph_group",
      graph_group_argument = "cell.graph.group",
      condition_argument = "cell.split.condition",
      condition_partition = "hierarchy_constrained",
      partition_schema_version = grouped$partition_schema_version,
      graph_method = "multimodal_WNN",
      clustering_method = "one_shared_walktrap_hierarchy_per_cell_type",
      final_partition_method =
        "condition_specific_finest_feasible_cut_of_shared_hierarchy",
      aggregation_method = "SCimplify_for_Seurat_with_membership",
      graph_scope = "one_independent_WNN_graph_per_cell_type",
      condition_scope =
        "all_conditions_joint_for_WNN_and_Walktrap_then_condition_specific_hierarchy_cut",
      membership_split_timing = "during_final_shared_hierarchy_cut_selection",
      modality_weighting = "adaptive_WNN_within_cell_type",
      hard_min_metacell_size = args$min_metacell_size,
      hard_min_metacells_per_stratum = args$min_metacells_per_stratum,
      temporary_combined_stratum = FALSE,
      gamma = args$gamma,
      k.knn = args$k.knn,
      kernel = args$kernel,
      inference_policy = paste(
        "Each broad cell type receives one independent multimodal WNN graph and Walktrap hierarchy;",
        "conditions share that graph/hierarchy; final condition-pure metacells are selected as",
        "condition-specific feasible cuts subject to hard size/count constraints"
      ),
      sample_metadata = "not_used_or_retained"
    )
  )
}
