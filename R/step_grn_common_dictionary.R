# Stage 1 dispatcher for common-dictionary condition GRNs.

#' Infer Pando regulatory evidence with automatic mode selection
#'
#' Two or more condition levels use global plus condition candidate discovery,
#' exact TF-peak-target union and fixed-dictionary condition GLMs. A missing or
#' single-level condition uses original per-cell-type Pando and creates no
#' condition coefficient.
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
    pando_args = list(),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start(
    "grn", outdir, progress, total_parts = 12L
  )
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  if (!is.list(pando_args)) stop("`pando_args` must be a list.", call. = FALSE)
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  if (identical(BPPARAM, TRUE)) {
    stop("`BPPARAM = TRUE` is invalid.", call. = FALSE)
  }
  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  design <- .rc_resolve_condition_design(object, condition_col)
  object <- design$object
  effective_condition_col <- design$condition_col
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
    "BPPARAM", "parallel"
  ))
  if (length(reserved)) {
    stop("`pando_args` cannot override workflow fields: ",
         paste(reserved, collapse = ", "), call. = FALSE)
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
    retired <- intersect(names(infer_args), c(
      "candidate_screen", "condition_mix", "condition_weight", "alpha",
      "nlambda", "lambda", "lambda_min_ratio", "outer_nfolds",
      "inner_nfolds", "lambda_selection", "scale", "engine_control",
      "comparison_conditions", "active_tol", "max_iter", "tol_objective",
      "tol_coef", "seed", "method", "sample_col", "cv_block_col"
    ))
    if (length(retired)) {
      stop(
        "Retired condition-GRN parameter(s): ", paste(retired, collapse = ", "),
        ". Use tf_cor, peak_cor, adjust_method='BH', padj_threshold=0.05, rank_action and min_residual_df.",
        call. = FALSE
      )
    }
    infer_args$adjust_method <- infer_args$adjust_method %||% "BH"
    infer_args$padj_threshold <- infer_args$padj_threshold %||% 0.05
    infer_args$rank_action <- infer_args$rank_action %||% "mark"
    infer_args$min_residual_df <- infer_args$min_residual_df %||% 1L
    infer_args$verbose <- infer_args$verbose %||%
      .rc_progress_enabled(progress)
    call_args$pando_infer_args <- infer_args
    call_args$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE
    call_args$progress_monitor <- monitor
    .rc_step_monitor_event(
      monitor, "pando_configuration",
      "configured two-stage exact-edge union and fixed-dictionary GLMs",
      current = 4L,
      context = list(
        tf_cor = infer_args$tf_cor %||% 0.1,
        peak_cor = infer_args$peak_cor %||% 0,
        adjust_method = infer_args$adjust_method,
        padj_threshold = infer_args$padj_threshold,
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
    params = list(
      requested_condition_col = design$requested_condition_col,
      condition_col = effective_condition_col,
      condition_levels = design$condition_levels,
      analysis_mode = design$analysis_mode,
      fallback_reason = design$fallback_reason,
      celltype_col = celltype_col,
      cell_type = cell_type,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      pando_args = c(extra_args, list(pando_infer_args = infer_args)),
      parallel = parallel,
      species = species
    )
  )
  class(answer) <- c("regcompass_grn_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(answer, file.path(outdir, "step_grn.rds"))
  answer
}
