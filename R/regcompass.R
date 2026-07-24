#' Run the canonical GRN-first RegCompass workflow
#'
#' @export
rc_run_regcompass <- function(
    object, gem, outdir, pfm, genome,
    fragment_files = FALSE,
    sample_col = NULL,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    model_mode = c("meta_module_gem", "full_gem"),
    medium_scenarios = NULL,
    metacell_args = list(),
    layer1_args = list(),
    pando_args = list(),
    layer2_args = list(),
    upstream_workers = NULL,
    layer2_workers = NULL,
    parallel_backend = c("auto", "serial", "snow", "multicore"),
    progress = getOption("RegCompassR.progress", TRUE),
    species = c("auto", "human", "mouse")) {
  model_mode <- match.arg(model_mode)
  parallel_backend <- match.arg(parallel_backend)
  progress <- .rc_progress_enabled(progress)
  old_progress_option <- options(RegCompassR.progress = progress)
  on.exit(do.call(options, old_progress_option), add = TRUE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  total_timer <- .rc_timing_start("total_workflow")
  overall <- .rc_progress_new(6L, "RegCompass total", progress)
  workflow_status <- "error"
  on.exit({
    if (!identical(workflow_status, "success")) {
      total_error <- .rc_timing_finish(
        total_timer, status = "error", outdir = outdir
      )
      .rc_write_execution_timing(total_error, outdir)
      .rc_progress_done(overall, "error")
    }
  }, add = TRUE)

  bundles <- list(
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args,
    layer2_args = layer2_args
  )
  invalid_bundles <- names(bundles)[!vapply(bundles, is.list, logical(1))]
  if (length(invalid_bundles)) {
    stop(
      "Workflow argument bundles must be lists: ",
      paste(invalid_bundles, collapse = ", "),
      call. = FALSE
    )
  }
  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  if (is.null(medium_scenarios)) {
    medium_scenarios <- rc_make_medium_scenarios(
      gem, scenario = "physiologic", species = species
    )
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  if (is.null(metacell_args$gamma)) metacell_args$gamma <- 75L
  saveRDS(gem$model_info %||% list(), file.path(outdir, "00_model_info.rds"))
  saveRDS(medium_scenarios, file.path(outdir, "00_medium_scenarios.rds"))

  upstream_config <- rc_parallel_config(
    workers = upstream_workers, backend = parallel_backend
  )
  layer2_config <- rc_parallel_config(
    workers = layer2_workers, backend = parallel_backend
  )
  upstream <- .rc_phase_bpparam(upstream_workers, parallel_backend)
  on.exit(.rc_release_bpparam(upstream), add = TRUE)

  stage_rows <- list()
  run_stage <- function(index, label, expression) {
    timer <- .rc_timing_start(label)
    value <- force(expression)
    timing <- .rc_timing_finish(timer, status = "success")
    stage_rows[[length(stage_rows) + 1L]] <<- timing
    .rc_progress_update(overall, index, label)
    value
  }

  step1 <- run_stage(1L, "single_cell_grn", rc_regcompass_step_grn(
    object = object, gem = gem,
    outdir = file.path(outdir, "01_single_cell_grn"),
    pfm = pfm, genome = genome,
    condition_col = condition_col, celltype_col = celltype_col,
    rna_assay = rna_assay, atac_assay = atac_assay,
    pando_args = pando_args, parallel = TRUE, BPPARAM = upstream,
    progress = progress
  ))
  step2 <- run_stage(2L, "condition_metacells", rc_regcompass_step_metacells(
    object = object, outdir = file.path(outdir, "02_condition_metacells"),
    sample_col = sample_col, condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay, fragment_files = fragment_files,
    metacell_args = metacell_args, progress = progress
  ))

  step3_args <- layer1_args
  local_fastcore_args <- step3_args$local_fastcore_args %||% list()
  if (is.null(local_fastcore_args$parallel)) {
    local_fastcore_args$parallel <- !identical(
      upstream_config$actual_backend, "serial"
    )
  }
  if (is.null(local_fastcore_args$workers)) {
    local_fastcore_args$workers <- upstream_config$workers
  }
  if (is.null(local_fastcore_args$backend)) {
    local_fastcore_args$backend <- upstream_config$actual_backend
  }
  step3_args$local_fastcore_args <- local_fastcore_args
  step3 <- run_stage(3L, "meta_modules", rc_regcompass_step_meta_modules(
    grn = step1, metacells = step2, gem = gem,
    outdir = file.path(outdir, "03_meta_modules"),
    layer1_args = step3_args, progress = progress
  ))
  step4 <- run_stage(4L, "layer1", rc_regcompass_step_layer1(
    metacells = step2, meta_modules = step3, gem = gem,
    outdir = file.path(outdir, "04_layer1"),
    regulatory_alpha = layer1_args$regulatory_alpha %||% 1,
    tau = layer1_args$tau %||% 0.20,
    gene_half_saturation = layer1_args$gene_half_saturation %||%
      getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE, BPPARAM = upstream, progress = progress
  ))
  .rc_release_bpparam(upstream)
  upstream <- FALSE

  lp_param <- .rc_phase_bpparam(layer2_workers, parallel_backend)
  step5 <- run_stage(5L, "layer2", tryCatch(
    rc_regcompass_step_layer2(
      layer1 = step4, meta_modules = step3, gem = gem,
      medium_scenarios = medium_scenarios,
      outdir = file.path(outdir, "05_layer2"),
      model_mode = model_mode, layer2_args = layer2_args,
      parallel = TRUE, BPPARAM = lp_param, progress = progress
    ),
    finally = .rc_release_bpparam(lp_param)
  ))
  step6 <- run_stage(6L, "results", rc_regcompass_step_results(
    grn = step1, metacells = step2, meta_modules = step3,
    layer1 = step4, layer2 = step5, gem = gem,
    outdir = file.path(outdir, "06_results"), species = species,
    progress = progress
  ))

  total_row <- .rc_timing_finish(total_timer, status = "success")
  execution_timing <- do.call(rbind, c(stage_rows, list(total_row)))
  .rc_write_execution_timing(execution_timing, outdir)

  result <- step6
  result$timing <- list(
    stages = execution_timing[
      execution_timing$stage != "total_workflow", , drop = FALSE
    ],
    total = total_row
  )
  result$params$execution_mode <- "one_shot"
  result$params$parallel_backend_requested <- parallel_backend
  result$params$parallel_backend_resolved <- list(
    upstream = upstream_config$actual_backend,
    layer2 = layer2_config$actual_backend
  )
  result$params$upstream_workers <- upstream_config$workers
  result$params$layer2_workers <- layer2_config$workers
  result$params$operating_system <- .Platform$OS.type
  saveRDS(result, file.path(outdir, "06_results", "regcompass_result.rds"))
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  workflow_status <- "success"
  .rc_progress_done(overall, paste0("complete in ", total_row$elapsed_hms))
  result
}
