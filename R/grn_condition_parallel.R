# Task helpers used directly by the canonical condition-GRN implementation.

.rc_condition_parallel_plan <- function(
    metadata, condition_types, condition_col, celltype_col, min_cells) {
  condition_types <- unique(as.character(condition_types))
  if (!length(condition_types)) {
    stop("No condition-GRN cell type was supplied.", call. = FALSE)
  }
  condition_value <- metadata[[condition_col]]
  declared_levels <- if (is.factor(condition_value)) {
    levels(condition_value)
  } else {
    unique(as.character(condition_value))
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
    observed <- unique(as.character(metadata[type_cells, condition_col]))
    condition_levels <- declared_levels[declared_levels %in% observed]
    counts <- vapply(condition_levels, function(level) {
      sum(as.character(metadata[type_cells, condition_col]) == level)
    }, integer(1))
    undersized <- condition_levels[counts < as.integer(min_cells)]
    if (length(undersized)) {
      detail <- paste0(undersized, "=", counts[undersized], collapse = ", ")
      stop(
        "Cell type `", type, "` has condition(s) below min_cells: ",
        detail, call. = FALSE
      )
    }
    if (length(condition_levels) < 2L) {
      stop(
        "Condition-GRN cell type `", type,
        "` must retain at least two conditions.", call. = FALSE
      )
    }
    cells_by_condition <- stats::setNames(
      lapply(condition_levels, function(level) {
        type_cells[as.character(metadata[type_cells, condition_col]) == level]
      }),
      condition_levels
    )
    plans[[type]] <- list(
      cell_type = type,
      conditions = condition_levels,
      reference_condition = condition_levels[[1L]],
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

.rc_pando_quiet_target_bpparam <- function(BPPARAM) {
  if (is.null(BPPARAM) || identical(BPPARAM, FALSE) ||
      !requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) {
    return(BPPARAM)
  }
  if (isTRUE(BiocParallel::bpisup(BPPARAM))) return(BPPARAM)
  workers <- max(1L, as.integer(BiocParallel::bpnworkers(BPPARAM)))
  quiet <- if (methods::is(BPPARAM, "SnowParam")) {
    BiocParallel::SnowParam(
      workers = workers,
      type = "SOCK",
      progressbar = FALSE,
      exportglobals = TRUE,
      exportvariables = TRUE
    )
  } else if (methods::is(BPPARAM, "MulticoreParam")) {
    BiocParallel::MulticoreParam(
      workers = workers,
      progressbar = FALSE
    )
  } else {
    BPPARAM
  }
  if (!identical(quiet, BPPARAM)) {
    attr(quiet, "regcompass_worker_limit") <-
      .rc_bpparam_worker_limit(BPPARAM, default = workers)
    config <- attr(BPPARAM, "regcompass_parallel_config", exact = TRUE)
    if (!is.null(config)) {
      attr(quiet, "regcompass_parallel_config") <- config
    }
  }
  quiet
}

.rc_condition_ridge_fit_task <- function(
    task, target_genes, condition_col, celltype_col, min_cells,
    pando_infer_args, inner_parallel = FALSE, PANDO_BPPARAM = NULL) {
  if (!is.list(task) || !inherits(task$grn, "GRNData") ||
      !is.character(task$cell_type) || length(task$cell_type) != 1L) {
    stop("Invalid condition-GRN E-star fit task.", call. = FALSE)
  }
  cell_type <- task$cell_type
  threshold <- suppressWarnings(as.numeric(pando_infer_args$padj_threshold))
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop("Condition Pando padj_threshold must be in (0, 1).",
         call. = FALSE)
  }
  reference_condition <- pando_infer_args$reference_condition %||% NULL
  show_progress <- .rc_progress_enabled(
    getOption("RegCompassR.progress", TRUE)
  )
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
    reference_condition = reference_condition,
    parallel = isTRUE(inner_parallel),
    parallel_scope = "target",
    overwrite = TRUE,
    verbose = show_progress
  )
  if (isTRUE(inner_parallel) && !is.null(PANDO_BPPARAM) &&
      !identical(PANDO_BPPARAM, FALSE)) {
    args$BPPARAM <- .rc_pando_quiet_target_bpparam(PANDO_BPPARAM)
  }
  if (show_progress) {
    message(
      "RegCompass grn condition detail | cell_type=", cell_type,
      ";phase=pando_Estar_z025_separate_inference_pipeline",
      ";targets_requested=", length(target_genes),
      ";reference_condition=",
      if (is.null(reference_condition)) "<first-retained>" else reference_condition,
      ";target_parallel=", isTRUE(inner_parallel),
      ";workers=", if (!is.null(args$BPPARAM) &&
          !identical(args$BPPARAM, FALSE)) {
        .rc_bpparam_worker_limit(args$BPPARAM, default = 1L)
      } else 1L
    )
  }
  fitted <- tryCatch(
    do.call(Pando::infer_condition_grn, args),
    error = function(error) {
      stop(
        "Condition-GRN E-star/separate-inference fit failed for cell type `",
        cell_type, "` during Pando target-level execution: ",
        conditionMessage(error), call. = FALSE
      )
    }
  )
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
    stop("Pando condition E-star fit was not returned for cell type `",
         cell_type, "`.", call. = FALSE)
  }
  if (!identical(as.character(fit$cell_type), cell_type)) {
    stop("Pando condition fit returned the wrong cell type.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(as.numeric(fit$padj_threshold), threshold))) {
    stop("Pando returned a condition fit with the wrong BH threshold.",
         call. = FALSE)
  }
  if (!is.null(reference_condition) &&
      !identical(as.character(fit$reference_condition),
                 as.character(reference_condition))) {
    stop("Pando returned a condition fit with the wrong predefined reference.",
         call. = FALSE)
  }
  if (!identical(as.character(fit$model_schema),
                 .RC_PANDO_CONDITION_GRN_MODEL_SCHEMA) ||
      !identical(as.character(fit$fit_engine),
                 .RC_PANDO_CONDITION_GRN_ENGINE) ||
      !identical(as.character(fit$inference_schema),
                 .RC_PANDO_CONDITION_INFERENCE_SCHEMA) ||
      !identical(as.character(fit$projection_policy),
                 .RC_PANDO_CONDITION_PROJECTION_POLICY) ||
      !is.list(fit$deviation_penalty) ||
      !identical(as.character(fit$deviation_penalty$family),
                 .RC_PANDO_CONDITION_PENALTY_FAMILY) ||
      !isTRUE(all.equal(
        as.numeric(fit$deviation_penalty$z),
        .RC_PANDO_CONDITION_SCHEME_E_Z,
        tolerance = 1e-15
      ))) {
    stop(
      "Pando returned a condition fit that is not the canonical E-star ",
      "z=0.25 separated-inference contract.", call. = FALSE
    )
  }
  invisible(.rc_require_pando_condition_grn_fit(fit))
  if (show_progress) {
    message(
      "RegCompass grn condition detail | cell_type=", cell_type,
      ";phase=pando_Estar_z025_separate_inference_complete",
      ";reference_condition=", as.character(fit$reference_condition),
      ";candidate_edges=", as.integer(fit$candidate_edge_count %||% NA_integer_),
      ";fit_edges=", as.integer(fit$fit_dictionary_edge_count %||% NA_integer_),
      ";condition_inference_estimable_rows=",
      sum(fit$coefficients$condition_inference_estimable %in% TRUE),
      ";edge_inference_estimable=",
      sum(fit$edge_inference$edge_inference_estimable %in% TRUE),
      ";regcompass_edges=",
      sum(fit$edge_inference$edge_supported %in% TRUE),
      ";targets_fitted=", length(unique(as.character(fit$target_genes)))
    )
  }
  list(cell_type = cell_type, grn = fitted, fit = fit)
}
