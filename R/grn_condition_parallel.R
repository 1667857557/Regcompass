# Task helpers used directly by the canonical condition-GRN implementation.

.rc_condition_network_label <- function(x) {
  value <- gsub("[^[:alnum:]_.-]+", "_", as.character(x))
  value[!nzchar(value)] <- "unnamed"
  value
}

.rc_condition_parallel_plan <- function(
    metadata, condition_types, condition_col, celltype_col, min_cells) {
  condition_types <- unique(as.character(condition_types))
  if (!length(condition_types)) {
    stop("No condition-GRN cell type was supplied.", call. = FALSE)
  }
  plans <- vector("list", length(condition_types))
  names(plans) <- condition_types
  for (type in condition_types) {
    type_cells <- rownames(metadata)[
      as.character(metadata[[celltype_col]]) == type
    ]
    if (!length(type_cells)) {
      stop("Condition-GRN cell type `", type, "` contains no cells.",
           call. = FALSE)
    }
    levels <- unique(as.character(metadata[type_cells, condition_col]))
    counts <- vapply(levels, function(level) {
      sum(as.character(metadata[type_cells, condition_col]) == level)
    }, integer(1))
    undersized <- levels[counts < as.integer(min_cells)]
    if (length(undersized)) {
      detail <- paste0(undersized, "=", counts[undersized], collapse = ", ")
      stop(
        "Cell type `", type, "` has condition(s) below min_cells: ",
        detail, call. = FALSE
      )
    }
    if (length(levels) < 2L) {
      stop(
        "Condition-GRN cell type `", type,
        "` must retain at least two conditions.", call. = FALSE
      )
    }
    cells_by_condition <- stats::setNames(lapply(levels, function(level) {
      type_cells[as.character(metadata[type_cells, condition_col]) == level]
    }), levels)
    plans[[type]] <- list(
      cell_type = type,
      conditions = levels,
      cells_by_condition = cells_by_condition,
      global_cells = unlist(cells_by_condition, use.names = FALSE)
    )
  }
  plans
}

.rc_condition_prepare_celltype_task <- function(
    task, atac_assay, rna_assay, pando_initiate_args,
    pando_motif_args, pfm, genome) {
  if (!is.list(task) || !inherits(task$object, "Seurat")) {
    stop("Invalid condition-GRN cell-type preparation task.", call. = FALSE)
  }
  filtered <- .rc_drop_zero_count_atac_features(
    task$object, atac_assay,
    paste0("Pando condition GRN for ", task$cell_type)
  )
  one <- filtered$object
  init <- list(object = one, peak_assay = atac_assay, rna_assay = rna_assay)
  init[names(pando_initiate_args)] <- NULL
  grn <- do.call(Pando::initiate_grn, c(init, pando_initiate_args))

  motif_args <- .rc_regcompass_motif_args(pando_motif_args)
  if (is.list(motif_args) && !is.null(motif_args$cache_dir)) {
    motif_args$cache_dir <- file.path(
      motif_args$cache_dir, .rc_safe_path_component(task$cell_type)
    )
  }
  motif <- list(object = grn, pfm = pfm, genome = genome)
  motif[names(motif_args)] <- NULL
  grn <- do.call(Pando::find_motifs, c(motif, motif_args))
  list(
    cell_type = task$cell_type,
    grn = grn,
    n_cells = ncol(one),
    n_removed_atac_features = filtered$n_removed %||% 0L
  )
}

.rc_condition_discovery_task <- function(
    task, target_genes, pando_infer_args) {
  if (!is.list(task) || !inherits(task$grn, "GRNData")) {
    stop("Invalid condition-GRN candidate-discovery task.", call. = FALSE)
  }
  edge <- Pando::discover_grn_edges(
    object = task$grn,
    genes = target_genes,
    cells = task$cells,
    source_label = task$source_label,
    source_type = task$source_type,
    tf_cor = pando_infer_args$tf_cor,
    peak_cor = pando_infer_args$peak_cor,
    rna_layer = pando_infer_args$rna_layer %||% "data",
    peak_layer = pando_infer_args$peak_layer %||% "data",
    peak_value_type = pando_infer_args$peak_value_type %||% "normalized",
    parallel = FALSE,
    verbose = FALSE
  )
  list(
    cell_type = task$cell_type,
    condition = task$condition,
    source_type = task$source_type,
    edge = edge
  )
}

.rc_condition_fit_task <- function(task, pando_infer_args) {
  if (!is.list(task) || !inherits(task$grn, "GRNData")) {
    stop("Invalid condition-GRN fixed-dictionary task.", call. = FALSE)
  }
  fitted <- Pando::fit_grn_from_edges(
    object = task$grn,
    edge_dictionary = task$dictionary,
    cells = task$cells,
    condition_label = task$condition,
    network_name = task$network_name,
    adjust_method = "BH",
    padj_threshold = 0.05,
    rank_action = pando_infer_args$rank_action,
    min_residual_df = pando_infer_args$min_residual_df,
    rna_layer = pando_infer_args$rna_layer %||% "data",
    peak_layer = pando_infer_args$peak_layer %||% "data",
    peak_value_type = pando_infer_args$peak_value_type %||% "normalized",
    parallel = FALSE,
    overwrite = TRUE,
    verbose = FALSE
  )
  network <- Pando::GetNetwork(fitted, network = task$network_name)
  list(
    cell_type = task$cell_type,
    condition = task$condition,
    network_name = task$network_name,
    network = network,
    coefficients = as.data.frame(stats::coef(network), stringsAsFactors = FALSE),
    fit = as.data.frame(Pando::gof(network), stringsAsFactors = FALSE)
  )
}
