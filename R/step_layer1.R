#' Build RNA plus Pando reaction support
#'
#' Each cell type uses its Stage 1 Pando route. Standard mode uses
#' the original Pando full-fit TF-by-ATAC projection and calculates no condition
#' coefficients. Both modes use exact native SuperCell membership. Parallel
#' projection jobs use at most `workers` processes and automatically shrink to
#' the number of independent tasks available.
#'
#' @param workers Total RegCompass worker cap, default 10. Windows uses SOCK/Snow
#'   workers and Linux/macOS uses multicore workers. Stage 4 never starts more
#'   workers than there are independent jobs. Users may increase or decrease the
#'   cap explicitly.
#' @export
rc_regcompass_step_layer1 <- function(
    grn, metacells, meta_modules, gem, outdir,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    workers = 10L,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer1", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  gpr_and_method <- match.arg(gpr_and_method)
  parallel_plan <- .rc_stage_parallel_plan(workers, argument = "workers")
  on.exit(.rc_release_bpparam(parallel_plan$BPPARAM), add = TRUE)
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
  layer1 <- .rc_cell_first_projection_layer1_v6(
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
    parallel = parallel_plan$parallel,
    BPPARAM = parallel_plan$BPPARAM
  )
  layer1$workflow_params <- params
  layer1$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  layer1$parallel_contract <- list(
    worker_limit = parallel_plan$workers,
    backend = parallel_plan$config$actual_backend,
    dynamic_task_sizing = TRUE
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
