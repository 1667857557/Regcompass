# Stage 1 dispatcher for standard and common-dictionary Pando GRNs.

#' Infer regulatory evidence with automatic mode selection
#'
#' Stage 1 filters the analysis cell set before normalization. Broad cell types
#' with at least two retained conditions use Pando's pooled/global plus
#' condition-local exact-edge union dictionary, fixed E-star z=0.25 production
#' fit, fusion-component joint-SE inference, condition-by-target BH, and
#' any-condition exact-edge RegCompass handoff. Cell types with one effective
#' condition use the independent K=1 standard Pando route.
#'
#' @param pando_args Pando configuration list. `min_cells` defaults to `500L`.
#'   Conditional inference controls supplied through `pando_infer_args` are
#'   `tf_cor`, `peak_cor`, `padj_threshold`, `rank_action`,
#'   `min_residual_df`, and optional `reference_condition`, plus the canonical
#'   layer fields. `reference_condition` is a predefined experimental-design
#'   coordinate for the K-condition contrast tree; it must be retained in every
#'   conditional cell type and must not be selected after inspecting GRN
#'   results. When omitted, Pando uses and records the first retained condition.
#'   Conditional ridge-CV, alternative-z and fusion-ratio controls are removed.
#'   Standard Pando retains its separate standard-route controls, including
#'   `ridge_control` when used.
#' @param target_rsq_threshold Full-data target R-squared threshold used by the
#'   standard Pando route. Defaults to `0.05`. Conditional E-star/JSE keeps the
#'   same target R-squared quantity as a diagnostic only; it does not gate the
#'   conditional exact-edge handoff.
#' @param workers Total RegCompass worker cap, default 10. Stage 1 uses one
#'   parallel level at a time and reuses this budget for target-level Pando work.
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
    target_rsq_threshold = 0.05,
    workers = 10L,
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
  target_rsq_threshold <- .rc_target_rsq_threshold(target_rsq_threshold)
  old_target_rsq_option <- getOption("RegCompassR.target_rsq_threshold")
  options(RegCompassR.target_rsq_threshold = target_rsq_threshold)
  on.exit({
    if (is.null(old_target_rsq_option)) {
      options(RegCompassR.target_rsq_threshold = NULL)
    } else {
      options(RegCompassR.target_rsq_threshold = old_target_rsq_option)
    }
  }, add = TRUE)

  parallel_plan <- .rc_stage_parallel_plan(workers, argument = "workers")
  on.exit(.rc_release_bpparam(parallel_plan$BPPARAM), add = TRUE)
  parallel <- parallel_plan$parallel
  BPPARAM <- parallel_plan$BPPARAM

  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
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
  condition_types <- cell_set$condition_pando_cell_types
  standard_types <- cell_set$standard_pando_cell_types

  cell_set$diagnostics$threshold_source <- "configurable_pando_args_min_cells"
  object@misc$regcompass_stage1_group_filter <- cell_set$diagnostics
  object@misc$regcompass_stage1_min_cells_contract <- list(
    min_cells = cell_set$min_cells,
    default_min_cells = .rc_stage1_min_cells_default,
    source = "pando_args$min_cells",
    fixed = FALSE,
    configurable = TRUE,
    analysis_mode = cell_set$analysis_mode,
    threshold_scope = if (length(condition_types)) {
      "condition_x_cell_type"
    } else {
      "cell_type"
    },
    condition_levels = cell_set$condition_levels,
    retained_cell_types = cell_set$retained_cell_types,
    condition_pando_cell_types = condition_types,
    standard_pando_cell_types = standard_types,
    applied_before_normalization = TRUE,
    passed_to_pando = TRUE
  )

  object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  zero_filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Stage 1 min_cells prefilter"
  )
  object <- zero_filtered$object
  object@misc$regcompass_stage1_zero_peak_filter <- zero_filtered$diagnostics
  object@misc$regcompass_stage1_fragment_policy <- list(
    policy = "clear_before_stage1",
    reason = "workflow uses in-memory ATAC counts and genome sequence"
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
    "target_rsq_threshold", "BPPARAM", "parallel", "workers"
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
  routed_infer_args <- .rc_route_pando_infer_args(
    infer_args,
    condition_types = condition_types,
    standard_types = standard_types
  )
  condition_infer_args <- routed_infer_args$condition
  standard_infer_args <- routed_infer_args$standard
  standard_infer_args$verbose <- standard_infer_args$verbose %||%
    .rc_progress_enabled(progress)

  dispatch_extra_args <- extra_args
  resources <- .rc_materialize_pando_resources(
    pfm = pfm,
    species = species,
    pando_initiate_args = dispatch_extra_args$pando_initiate_args %||%
      list(exclude_exons = TRUE),
    pando_motif_args = dispatch_extra_args$pando_motif_args %||% list()
  )
  pfm <- resources$pfm
  dispatch_extra_args$pando_initiate_args <- resources$pando_initiate_args
  dispatch_extra_args$pando_motif_args <- resources$pando_motif_args
  .rc_step_monitor_event(
    monitor,
    "pando_resource_materialization",
    "materialized Pando reference resources before worker dispatch",
    current = 5L,
    context = list(
      pfm = !is.null(pfm),
      regions = !is.null(dispatch_extra_args$pando_initiate_args$regions),
      motif_tfs = !is.null(dispatch_extra_args$pando_motif_args$motif_tfs),
      worker_limit = parallel_plan$workers,
      backend = parallel_plan$config$actual_backend,
      target_rsq_threshold = target_rsq_threshold
    )
  )

  grn_result <- .rc_with_step_diagnostics(
    .rc_fit_pando_by_celltype_route(
      object = object, gem = gem, outdir = outdir, genome = genome,
      pfm = pfm, species = species,
      condition_col = effective_condition_col,
      celltype_col = celltype_col,
      condition_types = condition_types,
      standard_types = standard_types,
      rna_assay = rna_assay, atac_assay = atac_assay,
      extra_args = dispatch_extra_args,
      condition_infer_args = condition_infer_args,
      standard_infer_args = standard_infer_args,
      parallel = parallel, BPPARAM = BPPARAM,
      progress_monitor = monitor
    ),
    monitor
  )

  grn_result$analysis_mode <- cell_set$analysis_mode
  grn_result$requested_condition_col <- design$requested_condition_col
  grn_result$effective_condition_col <- effective_condition_col
  grn_result$condition_levels <- design$condition_levels
  grn_result$fallback_reason <- if (identical(
    cell_set$analysis_mode, "mixed_pando"
  )) "cell_type_specific_condition_count" else design$fallback_reason
  grn_result$rna_assay <- rna_assay
  grn_result$atac_assay <- atac_assay
  grn_result$target_rsq_threshold <- target_rsq_threshold
  grn_result$pando_infer_argument_routing <- routed_infer_args$diagnostics

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
      condition_pando_cell_types = condition_types,
      standard_pando_cell_types = standard_types,
      diagnostics = cell_set$diagnostics,
      analysis_mode = cell_set$analysis_mode,
      condition_levels = cell_set$condition_levels
    ),
    params = list(
      requested_condition_col = design$requested_condition_col,
      condition_col = effective_condition_col,
      condition_levels = design$condition_levels,
      analysis_mode = cell_set$analysis_mode,
      fallback_reason = grn_result$fallback_reason,
      cell_type_analysis_mode = grn_result$cell_type_analysis_mode,
      celltype_col = celltype_col,
      cell_type = cell_set$retained_cell_types,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      target_rsq_threshold = target_rsq_threshold,
      pando_args = c(extra_args, list(
        pando_infer_args = infer_args,
        condition_pando_infer_args = condition_infer_args,
        standard_pando_infer_args = standard_infer_args,
        infer_argument_routing = routed_infer_args$diagnostics
      )),
      pando_resource_materialization = list(
        scope = "controller_before_parallel_dispatch",
        pfm = TRUE,
        regions = TRUE,
        motif_tfs = TRUE,
        stored_in_step_params = FALSE
      ),
      workers = parallel_plan$workers,
      parallel = parallel,
      parallel_backend = parallel_plan$config$actual_backend,
      pando_execution_plan = grn_result$pando_execution_plan,
      species = species,
      n_input_cells = as.integer(n_input),
      n_stage_cells = as.integer(length(cell_set$retained_cells)),
      cell_set_contract = "stage1_exact_cell_ids"
    )
  )
  class(answer) <- c("regcompass_grn_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(answer, file.path(outdir, "step_grn.rds"))
  answer
}
