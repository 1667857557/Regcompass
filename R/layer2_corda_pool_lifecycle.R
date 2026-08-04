# Keep one package-managed worker pool alive across all CORDA stages.

.rc_regcompass_step_layer2_corda_pool_base <- rc_regcompass_step_layer2

rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  model_params <- if (is.list(layer2_args)) {
    layer2_args$model_params %||% list()
  } else {
    list()
  }
  requested <- as.character(model_params$model_completion %||% "fastcore")
  is_corda <- length(requested) == 1L && !is.na(requested) &&
    requested %in% c("corda", "corda_like")
  managed <- FALSE
  if (isTRUE(is_corda) && isTRUE(parallel) && is.null(BPPARAM)) {
    BPPARAM <- rc_default_bpparam()
    if (!is.null(BPPARAM) &&
        requireNamespace("BiocParallel", quietly = TRUE)) {
      BiocParallel::bpstart(BPPARAM)
      managed <- TRUE
      on.exit({
        .rc_release_bpparam(BPPARAM)
        invisible(gc(verbose = FALSE, full = TRUE))
      }, add = TRUE)
    }
  }
  answer <- .rc_regcompass_step_layer2_corda_pool_base(
    layer1 = layer1,
    meta_modules = meta_modules,
    gem = gem,
    medium_scenarios = medium_scenarios,
    outdir = outdir,
    model_mode = model_mode,
    layer2_args = layer2_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
  answer$params$corda_package_managed_worker_pool <- managed
  if (isTRUE(is_corda)) {
    rc_export_microcompass(answer, outdir)
    saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  }
  answer
}
