# Keep one worker pool alive across all CORDA stages.

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
  pool_origin <- "not_used"
  pool_started_here <- FALSE
  if (isTRUE(is_corda) && isTRUE(parallel)) {
    if (is.null(BPPARAM)) {
      BPPARAM <- rc_default_bpparam()
      pool_origin <- if (is.null(BPPARAM)) "serial_fallback" else "package_default"
    } else {
      pool_origin <- "caller_supplied"
    }
    if (!is.null(BPPARAM)) {
      if (!requireNamespace("BiocParallel", quietly = TRUE) ||
          !methods::is(BPPARAM, "BiocParallelParam")) {
        stop(
          "CORDA parallel execution requires a BiocParallelParam object.",
          call. = FALSE
        )
      }
      was_started <- isTRUE(BiocParallel::bpisup(BPPARAM))
      if (!was_started) {
        thread_state <- .rc_set_internal_single_thread()
        BiocParallel::bpstart(BPPARAM)
        pool_started_here <- TRUE
        on.exit({
          .rc_release_bpparam(BPPARAM)
          .rc_restore_internal_threads(thread_state)
          invisible(gc(verbose = FALSE, full = TRUE))
        }, add = TRUE)
      }
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
  answer$params$corda_worker_pool_origin <- pool_origin
  answer$params$corda_worker_pool_started_by_layer2 <- pool_started_here
  if (isTRUE(is_corda)) {
    answer$params$structural_completion <- "corda"
    answer$params$structural_completion_algorithm <-
      "Schultz_Qutub_CORDA_2016_three_stage_dependency_assessment"
    answer$method <- paste(
      "microCOMPASS directional LP on cell-type-specific medium models",
      "reconstructed by original three-stage CORDA"
    )
    rc_export_microcompass(answer, outdir)
    saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  }
  answer
}
