#' Run the canonical significant-Pando-target RegCompass workflow
#'
#' Stage 3 constructs biological meta-modules and a deduplicated merged
#' reaction catalogue. Within each condition by cell-type group, Human-GEM
#' metabolic genes with significant Pando TF-peak-gene evidence define the
#' supported gene set, and complete GPR branches define core reactions. Stage 3
#' does not run FASTCORE and does not create a GEM. With
#' `model_mode = "meta_module_gem"`, Stage 5 constructs one union GEM per medium
#' scenario and performs the only FASTCORE completion on that shared structure.
#'
#' @param meta_module_args Stage 3 annotation-expansion arguments:
#'   `subsystem_table`, `expansion_mode`, and `max_iterations`.
#' @param layer1_args Stage 4 integrated-evidence arguments:
#'   `regulatory_alpha`, `gpr_and_method`, and `gene_half_saturation`.
#'   `gpr_and_method` accepts COMPASS-compatible `"min"`, `"median"`, or
#'   `"mean"`; RegCompass defaults to `"min"`.
#' @param upstream_workers Worker count for GRN inference and Layer 1
#'   reaction-expression calculation. Defaults to 6. Set to 1 for serial
#'   upstream execution.
#' @param layer2_workers Worker count for Layer 2 LP scoring. Defaults to 30.
#'   Set to 1 for serial Layer 2 execution.
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
    meta_module_args = list(),
    layer2_args = list(),
    upstream_workers = 6L,
    layer2_workers = 30L,
    species = c("auto", "human", "mouse"),
    progress = getOption("RegCompassR.progress", TRUE)) {
  model_mode <- match.arg(model_mode)
  progress <- .rc_progress_enabled(progress)
  old_progress_option <- options(RegCompassR.progress = progress)
  on.exit(do.call(options, old_progress_option), add = TRUE)

  thread_state <- .rc_set_internal_single_thread()
  on.exit(.rc_restore_internal_threads(thread_state), add = TRUE)

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
    meta_module_args = meta_module_args,
    layer2_args = layer2_args
  )
  invalid_bundles <- names(bundles)[
    !vapply(bundles, is.list, logical(1))
  ]
  if (length(invalid_bundles)) {
    stop(
      "Workflow argument bundles must be lists: ",
      paste(invalid_bundles, collapse = ", "),
      call. = FALSE
    )
  }
  unknown_meta_module <- setdiff(
    names(meta_module_args),
    c("subsystem_table", "expansion_mode", "max_iterations")
  )
  if (length(unknown_meta_module)) {
    stop(
      "Unknown `meta_module_args` fields: ",
      paste(unknown_meta_module, collapse = ", "),
      call. = FALSE
    )
  }
  unknown_layer1 <- setdiff(
    names(layer1_args),
    c("regulatory_alpha", "gpr_and_method", "gene_half_saturation")
  )
  if (length(unknown_layer1)) {
    stop(
      "Unknown `layer1_args` fields: ",
      paste(unknown_layer1, collapse = ", "),
      ". The retired `tau`/Boltzmann GPR-AND API has been removed.",
      call. = FALSE
    )
  }
  gpr_and_method <- match.arg(
    as.character(layer1_args$gpr_and_method %||% "min"),
    c("min", "median", "mean")
  )

  pando_infer_args <- pando_args$pando_infer_args %||% list()
  if (!is.list(pando_infer_args)) {
    stop("`pando_args$pando_infer_args` must be a list.", call. = FALSE)
  }
  if (!is.null(pando_infer_args$parallel) &&
      !identical(pando_infer_args$parallel, FALSE)) {
    warning(
      "Ignoring `pando_args$pando_infer_args$parallel`; canonical outer ",
      "parallelism requires each Pando fit to remain single-threaded.",
      call. = FALSE
    )
  }
  pando_infer_args$parallel <- FALSE
  pando_args$pando_infer_args <- pando_infer_args

  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  if (is.null(medium_scenarios)) {
    medium_scenarios <- rc_make_medium_scenarios(
      gem, scenario = "physiologic", species = species
    )
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  if (is.null(metacell_args$gamma)) metacell_args$gamma <- 30L
  saveRDS(gem$model_info %||% list(), file.path(outdir, "00_model_info.rds"))
  saveRDS(medium_scenarios, file.path(outdir, "00_medium_scenarios.rds"))

  upstream_config <- .rc_stage_worker_config(
    upstream_workers, argument = "upstream_workers"
  )
  layer2_config <- .rc_stage_worker_config(
    layer2_workers, argument = "layer2_workers"
  )

  stage_rows <- list()
  run_stage <- function(index, label, expression) {
    timer <- .rc_timing_start(label)
    value <- force(expression)
    timing <- .rc_timing_finish(timer, status = "success")
    stage_rows[[length(stage_rows) + 1L]] <<- timing
    .rc_progress_update(overall, index, label)
    invisible(gc(verbose = FALSE, full = TRUE))
    value
  }

  step1 <- run_stage(
    1L,
    "single_cell_grn",
    .rc_with_stage_workers(
      upstream_config$workers,
      argument = "upstream_workers",
      FUN = function(upstream, config) {
        rc_regcompass_step_grn(
          object = object,
          gem = gem,
          outdir = file.path(outdir, "01_single_cell_grn"),
          pfm = pfm,
          genome = genome,
          condition_col = condition_col,
          celltype_col = celltype_col,
          rna_assay = rna_assay,
          atac_assay = atac_assay,
          pando_args = pando_args,
          parallel = !identical(config$actual_backend, "serial"),
          BPPARAM = upstream,
          progress = progress
        )
      }
    )
  )

  step2 <- run_stage(
    2L,
    "condition_metacells",
    rc_regcompass_step_metacells(
      object = object,
      outdir = file.path(outdir, "02_condition_metacells"),
      sample_col = sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      fragment_files = fragment_files,
      metacell_args = metacell_args,
      progress = progress
    )
  )

  step3 <- run_stage(
    3L,
    "meta_modules",
    rc_regcompass_step_meta_modules(
      grn = step1,
      metacells = step2,
      gem = gem,
      outdir = file.path(outdir, "03_meta_modules"),
      meta_module_args = meta_module_args,
      progress = progress
    )
  )

  step4 <- run_stage(
    4L,
    "layer1",
    .rc_with_stage_workers(
      upstream_config$workers,
      argument = "upstream_workers",
      FUN = function(upstream, config) {
        rc_regcompass_step_layer1(
          metacells = step2,
          meta_modules = step3,
          gem = gem,
          outdir = file.path(outdir, "04_layer1"),
          regulatory_alpha = layer1_args$regulatory_alpha %||% 1,
          gpr_and_method = gpr_and_method,
          gene_half_saturation = layer1_args$gene_half_saturation %||%
            getOption("RegCompassR.cpm_half_saturation", 1),
          parallel = !identical(config$actual_backend, "serial"),
          BPPARAM = upstream,
          progress = progress
        )
      }
    )
  )

  step5 <- run_stage(
    5L,
    "layer2",
    .rc_with_stage_workers(
      layer2_config$workers,
      argument = "layer2_workers",
      FUN = function(layer2_param, config) {
        rc_regcompass_step_layer2(
          layer1 = step4,
          meta_modules = step3,
          gem = gem,
          medium_scenarios = medium_scenarios,
          outdir = file.path(outdir, "05_layer2"),
          model_mode = model_mode,
          layer2_args = layer2_args,
          parallel = !identical(config$actual_backend, "serial"),
          BPPARAM = layer2_param,
          progress = progress
        )
      }
    )
  )

  step6 <- run_stage(
    6L,
    "results",
    rc_regcompass_step_results(
      grn = step1,
      metacells = step2,
      meta_modules = step3,
      layer1 = step4,
      layer2 = step5,
      gem = gem,
      outdir = file.path(outdir, "06_results"),
      species = species,
      progress = progress
    )
  )

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
  result$params$parallel_backend_requested <- "auto"
  result$params$parallel_backend_resolved <- list(
    upstream = upstream_config$actual_backend,
    layer2 = layer2_config$actual_backend
  )
  result$params$upstream_workers <- upstream_config$workers
  result$params$layer2_workers <- layer2_config$workers
  result$params$internal_threads_per_task <- 1L
  result$params$pando_internal_parallel <- FALSE
  result$params$parallel_worker_lifecycle <-
    "stage_scoped_create_start_stop_release_full_gc"
  result$params$parallel_stage_groups <- list(
    upstream = c("single_cell_grn", "layer1"),
    layer2 = "layer2"
  )
  result$params$operating_system <- .Platform$OS.type
  saveRDS(result, file.path(outdir, "06_results", "regcompass_result.rds"))
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  workflow_status <- "success"
  .rc_progress_done(overall, paste0("complete in ", total_row$elapsed_hms))
  result
}
