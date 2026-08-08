#' Build RNA plus Pando reaction support
#'
#' Each cell type uses its Stage 1 Pando route. Standard mode uses
#' the original Pando full-fit TF-by-ATAC projection and calculates no condition
#' coefficients. Both modes use exact native SuperCell membership. Parallel work
#' is bounded by one global `workers` budget and automatically shrinks to the
#' number of independent tasks.
#'
#' @param workers Global RegCompass worker upper bound. `NULL` uses
#' `options(RegCompassR.workers)`, scheduler allocation, or detected cores.
#' @export
rc_regcompass_step_layer1 <- function(
    grn, metacells, meta_modules, gem, outdir,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    workers = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer1", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  gpr_and_method <- match.arg(gpr_and_method)
  worker_config <- rc_parallel_config(workers = workers, backend = "auto")
  BPPARAM <- if (worker_config$workers > 1L &&
                 !identical(worker_config$actual_backend, "serial")) {
    .rc_task_bpparam(workers = worker_config$worker_budget)
  } else {
    FALSE
  }
  on.exit(.rc_release_bpparam(BPPARAM), add = TRUE)

  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!identical(metacells$params, meta_modules$workflow_params)) {
    stop("Metacell and meta-module stages use different workflow settings.",
         call. = FALSE)
  }
  if (!identical(grn$params$analysis_mode, metacells$params$analysis_mode)) {
    stop("GRN and metacell stages resolved different analysis modes.",
         call. = FALSE)
  }
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_require_stage_gem(grn, gem, "grn")
  params <- metacells$params
  layer1 <- .rc_cell_first_projection_layer1(
    grn_result = grn$grn_result,
    metacell_object = metacells$metacell_object,
    membership = metacells$pooled$membership,
    gem = gem,
    metacell_meta = metacells$pooled$metacell_meta,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    gpr_and_method = gpr_and_method,
    gene_half_saturation = gene_half_saturation,
    parallel = worker_config$workers > 1L &&
      !identical(worker_config$actual_backend, "serial"),
    BPPARAM = BPPARAM
  )
  layer1$workflow_params <- params
  layer1$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  layer1$parallel_execution <- list(
    workers = worker_config$worker_budget,
    backend = worker_config$actual_backend,
    policy = "task_count_capped_by_global_worker_budget"
  )
  class(layer1) <- c("regcompass_layer1_step", "list")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer1 <- .rc_step_monitor_finish(layer1, monitor)
  saveRDS(layer1, file.path(outdir, "step_layer1.rds"))
  layer1
}
