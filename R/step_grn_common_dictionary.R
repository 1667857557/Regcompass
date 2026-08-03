# Stage 1 dispatcher for standard and common-dictionary Pando GRNs.

#' Infer regulatory evidence with automatic mode selection
#'
#' Stage 1 filters the analysis cell set before normalization. With at least two
#' retained conditions in a broad cell type, Pando performs pooled and
#' condition-specific candidate discovery, freezes the exact TF-peak-target
#' union, and fits the same unscaled Gaussian identity model in every condition.
#' With no condition or one effective condition, the original per-cell-type
#' Pando workflow is used without constructing condition coefficients.
#'
#' @param fragment_files Preserve existing Signac fragment references when TRUE.
#' The default FALSE clears them before Stage 1 because the workflow uses the
#' in-memory peak matrix and genome sequence.
#' @export
rc_regcompass_step_grn <- function(
    object, gem, outdir, genome,
    pfm = NULL,
    species = c("auto", "human", "mouse"),
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    pando_args = list(),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start(
    "grn", outdir, progress, total_parts = 12L
  )
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)

  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!is.list(pando_args)) {
    stop("`pando_args` must be a list.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  if (identical(BPPARAM, TRUE)) {
    stop("`BPPARAM = TRUE` is invalid.", call. = FALSE)
  }

  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  preserve_fragments <- .rc_validate_stage1_fragment_policy(fragment_files)
  n_input <- ncol(object)

  cell_set <- .rc_build_stage_analysis_cell_set(
    object = object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    pando_args = pando_args
  )
  object <- cell_set$object
  pando_args <- cell_set$pando_args
  cell_type <- cell_set$retained_cell_types

  cell_set$diagnostics$threshold_source <- "fixed_pando_args_min_cells"
  object@misc$regcompass_stage1_group_filter <- cell_set$diagnostics
  object@misc$regcompass_stage1_min_cells_contract <- list(
    min_cells = cell_set$min_cells,
    source = "pando_args$min_cells",
    fixed = TRUE,
    analysis_mode = cell_set$analysis_mode,
    threshold_scope = if (identical(cell_set$analysis_mode, "condition_grn")) {
      "condition_x_cell_type_independent"
    } else {
      "cell_type"
    },
    condition_levels = cell_set$condition_levels,
    retained_cell_types = cell_set$retained_cell_types,
    skipped_condition_cell_types = cell_set$skipped_condition_cell_types,
    applied_before_normalization = TRUE,
    passed_to_pando = TRUE
  )

  if (!preserve_fragments) {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  }
  zero_filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Stage 1 min_cells prefilter"
  )
  object <- zero_filtered$object
  object@misc$regcompass_stage1_zero_peak_filter <- zero_filtered$diagnostics
  object@misc$regcompass_stage1_fragment_policy <- list(
    fragment_files = preserve_fragments,
    policy = if (preserve_fragments) "preserve" else "clear_before_stage1"
  )

  motif_args <- pando_args$pando_motif_args %||% list()
  if (!is.list(motif_args)) {
    stop("`pando_motif_args` must be a list.", call. = FALSE)
  }
  if (.rc_pando_supports_motif_cache()) {
    motif_args$cache_dir <- motif_args$cache_dir %||%
      file.path(outdir, "motif_cache")
    motif_args$reuse_cache <- motif_args$reuse_cache %||% TRUE
  }
  pando_args$pando_motif_args <- motif_args

  design <- .rc_resolve_condition_design(object, condition_col)
  object <- design$object
  effective_condition_col <- design$condition_col
  if (!identical(design$analysis_mode, cell_set$analysis_mode)) {
    stop(
      "Stage 1 cell filtering and condition design resolved different analysis modes.",
      call. = FALSE
    )
  }

  object <- .rc_normalize_single_cell_grn_object(
    object,
    condition_col = effective_condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )

  reserved <- intersect(names(pando_args), c(
    "object", "gem", "outdir", "genome", "pfm", "species",
    "condition_col", "celltype_col", "cell_type", "rna_assay", "atac_assay",
    "BPPARAM", "parallel", "fragment_files"
  ))
  if (length(reserved)) {
    stop(
      "`pando_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "), call. = FALSE
    )
  }

  infer_args <- pando_args$pando_infer_args %||% list()
  if (!is.list(infer_args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  extra_args <- pando_args
  extra_args$pando_infer_args <- NULL
  defaults <- list(
    object = object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    condition_col = effective_condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  call_args <- c(
    defaults[setdiff(names(defaults), names(extra_args))],
    extra_args
  )

  if (identical(design$analysis_mode, "condition_grn")) {
    standard_only <- intersect(
      names(call_args), c("min_abs_estimate", "min_model_rsq")
    )
    if (length(standard_only)) {
      message(
        "Ignoring standard-Pando-only edge filters in multi-condition mode: ",
        paste(standard_only, collapse = ", "),
        ". The condition penalty is fixed to estimable BH padj < 0.05."
      )
      call_args[standard_only] <- NULL
    }
    retired <- intersect(names(infer_args), c(
      "candidate_screen", "condition_mix", "condition_weight", "alpha",
      "nlambda", "lambda", "lambda_min_ratio", "outer_nfolds",
      "inner_nfolds", "lambda_selection", "scale", "engine_control",
      "comparison_conditions", "active_tol", "max_iter", "tol_objective",
      "tol_coef", "seed", "method", "sample_col", "cv_block_col"
    ))
    if (length(retired)) {
      stop(
        "Retired condition-GRN parameter(s): ",
        paste(retired, collapse = ", "),
        ". Use tf_cor, peak_cor, adjust_method='BH', padj_threshold=0.05, rank_action and min_residual_df.",
        call. = FALSE
      )
    }
    infer_args <- utils::modifyList(list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L,
      verbose = .rc_progress_enabled(progress)
    ), infer_args)
    if (!identical(toupper(as.character(infer_args$adjust_method)), "BH") ||
        !isTRUE(all.equal(as.numeric(infer_args$padj_threshold), 0.05))) {
      stop("Canonical RegCompass condition effects require BH padj < 0.05.",
           call. = FALSE)
    }
    call_args$pando_infer_args <- infer_args
    call_args$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE
    call_args$progress_monitor <- monitor
    .rc_step_monitor_event(
      monitor, "pando_configuration",
      "configured two-stage exact-edge union and fixed-dictionary GLMs",
      current = 4L,
      context = list(
        tf_cor = infer_args$tf_cor,
        peak_cor = infer_args$peak_cor,
        adjust_method = "BH",
        padj_threshold = 0.05,
        rank_action = infer_args$rank_action,
        scale = FALSE,
        interaction = ":"
      )
    )
    grn_result <- .rc_with_step_diagnostics(
      do.call(.rc_fit_condition_grns_by_cell_type, call_args), monitor
    )
  } else {
    infer_args$verbose <- infer_args$verbose %||%
      .rc_progress_enabled(progress)
    call_args$pando_infer_args <- infer_args
    call_args$parallel <- isTRUE(parallel)
    call_args$progress_monitor <- monitor
    .rc_step_monitor_event(
      monitor, "standard_pando",
      "dispatching original per-cell-type Pando GRN workflow", current = 5L
    )
    grn_result <- .rc_with_step_diagnostics(
      do.call(.rc_fit_standard_pando_by_cell_type, call_args), monitor
    )
  }

  grn_result$analysis_mode <- design$analysis_mode
  grn_result$requested_condition_col <- design$requested_condition_col
  grn_result$effective_condition_col <- effective_condition_col
  grn_result$condition_levels <- design$condition_levels
  grn_result$fallback_reason <- design$fallback_reason
  grn_result$rna_assay <- rna_assay
  grn_result$atac_assay <- atac_assay

  answer <- list(
    grn_result = grn_result,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    cell_filter = list(
      source = "stage1_min_cells_before_normalization",
      min_cells = cell_set$min_cells,
      n_input_cells = as.integer(n_input),
      n_retained_cells = as.integer(length(cell_set$retained_cells)),
      n_removed_cells = as.integer(n_input - length(cell_set$retained_cells)),
      retained_cells = cell_set$retained_cells,
      retained_cell_types = cell_set$retained_cell_types,
      skipped_condition_cell_types = cell_set$skipped_condition_cell_types,
      diagnostics = cell_set$diagnostics,
      analysis_mode = cell_set$analysis_mode,
      condition_levels = cell_set$condition_levels
    ),
    params = list(
      requested_condition_col = design$requested_condition_col,
      condition_col = effective_condition_col,
      condition_levels = design$condition_levels,
      analysis_mode = design$analysis_mode,
      fallback_reason = design$fallback_reason,
      celltype_col = celltype_col,
      cell_type = cell_set$retained_cell_types,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      fragment_files = preserve_fragments,
      pando_args = c(extra_args, list(pando_infer_args = infer_args)),
      parallel = parallel,
      species = species,
      n_input_cells = as.integer(n_input),
      n_stage_cells = as.integer(length(cell_set$retained_cells)),
      cell_set_contract = "stage1_exact_cell_ids_v1"
    )
  )
  class(answer) <- c("regcompass_grn_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(answer, file.path(outdir, "step_grn.rds"))
  answer
}
