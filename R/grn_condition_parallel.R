# Task helpers used directly by the canonical condition-GRN implementation.

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
  n_cells <- ncol(one)
  n_removed <- filtered$n_removed %||% 0L
  task$object <- NULL
  filtered <- NULL
  init <- list(object = one, peak_assay = atac_assay, rna_assay = rna_assay)
  init[names(pando_initiate_args)] <- NULL
  grn <- do.call(Pando::initiate_grn, c(init, pando_initiate_args))
  init <- NULL
  one <- NULL
  invisible(gc(verbose = FALSE, full = TRUE))

  motif_args <- .rc_regcompass_motif_args(pando_motif_args)
  if (is.list(motif_args) && !is.null(motif_args$cache_dir)) {
    motif_args$cache_dir <- file.path(
      motif_args$cache_dir, .rc_safe_path_component(task$cell_type)
    )
  }
  motif <- list(object = grn, pfm = pfm, genome = genome)
  motif[names(motif_args)] <- NULL
  grn <- do.call(Pando::find_motifs, c(motif, motif_args))
  motif <- NULL
  motif_args <- NULL
  invisible(gc(verbose = FALSE, full = TRUE))
  list(
    cell_type = task$cell_type,
    grn = grn,
    n_cells = n_cells,
    n_removed_atac_features = n_removed
  )
}

.rc_condition_multitask_fit_task <- function(
    task, target_genes, condition_col, celltype_col, min_cells,
    pando_infer_args, inner_parallel = FALSE, PANDO_BPPARAM = NULL) {
  if (!is.list(task) || !inherits(task$grn, "GRNData") ||
      !is.character(task$cell_type) || length(task$cell_type) != 1L) {
    stop("Invalid condition-GRN multi-task fit task.", call. = FALSE)
  }
  cell_type <- task$cell_type
  ridge_control <- pando_infer_args$condition_ridge_control %||% list()
  threshold <- suppressWarnings(as.numeric(pando_infer_args$padj_threshold))
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop("Condition Pando padj_threshold must be in (0, 1).",
         call. = FALSE)
  }
  args <- list(
    object = task$grn,
    cell_type_col = celltype_col,
    condition_col = condition_col,
    cell_type = cell_type,
    genes = target_genes,
    network_name = "regcompass_condition_grn",
    rna_layer = pando_infer_args$rna_layer %||% "data",
    peak_layer = pando_infer_args$peak_layer %||% "data",
    peak_value_type = pando_infer_args$peak_value_type %||% "normalized",
    tf_cor = pando_infer_args$tf_cor,
    peak_cor = pando_infer_args$peak_cor,
    min_cells_per_condition = as.integer(min_cells),
    small_condition_action = "error",
    adjust_method = "BH",
    padj_threshold = threshold,
    rank_action = pando_infer_args$rank_action,
    min_residual_df = pando_infer_args$min_residual_df,
    parallel = isTRUE(inner_parallel),
    parallel_scope = "target",
    overwrite = TRUE,
    fallback_args = list(condition_ridge_control = ridge_control),
    verbose = FALSE
  )
  if (isTRUE(inner_parallel) && !is.null(PANDO_BPPARAM) &&
      !identical(PANDO_BPPARAM, FALSE)) {
    args$BPPARAM <- PANDO_BPPARAM
  }
  fitted <- do.call(Pando::infer_condition_grn, args)
  args <- NULL
  task$grn <- NULL
  invisible(gc(verbose = FALSE, full = TRUE))
  fits <- Pando::condition_grn_fit(fitted)
  if (inherits(fits, "ConditionGRNFit")) {
    fit <- fits
  } else if (is.list(fits) && cell_type %in% names(fits)) {
    fit <- fits[[cell_type]]
  } else if (is.list(fits) && length(fits) == 1L &&
             inherits(fits[[1L]], "ConditionGRNFit")) {
    fit <- fits[[1L]]
  } else {
    stop("Pando multi-task condition fit was not returned for cell type `",
         cell_type, "`.", call. = FALSE)
  }
  if (!identical(as.character(fit$cell_type), cell_type)) {
    stop("Pando multi-task condition fit returned the wrong cell type.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(as.numeric(fit$padj_threshold), threshold))) {
    stop("Pando returned a condition fit with the wrong BH threshold.",
         call. = FALSE)
  }
  list(
    cell_type = cell_type,
    grn = fitted,
    fit = fit
  )
}
