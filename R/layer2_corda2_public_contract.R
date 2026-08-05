# Public Layer 2 contract for the corrected Python CORDA2 implementation.

.rc_build_celltype_medium_union_gem_cache_before_corda2_public <-
  .rc_build_celltype_medium_union_gem_cache

.rc_corda2_move_cache_files <- function(cache) {
  if (!length(cache)) return(cache)
  old_files <- unique(vapply(cache, function(entry) {
    as.character(entry$file %||% "")
  }, character(1)))
  old_files <- old_files[nzchar(old_files)]
  file_map <- stats::setNames(old_files, old_files)
  for (old_file in old_files) {
    old_dir <- dirname(old_file)
    new_dir <- file.path(dirname(old_dir), "corda2")
    dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)
    new_file <- file.path(new_dir, basename(old_file))
    if (!identical(normalizePath(old_file, mustWork = FALSE),
                   normalizePath(new_file, mustWork = FALSE))) {
      moved <- isTRUE(file.rename(old_file, new_file))
      if (!moved) {
        copied <- isTRUE(file.copy(old_file, new_file, overwrite = TRUE))
        if (!copied) {
          stop("Could not move CORDA2 cache model to its dedicated directory.",
               call. = FALSE)
        }
        unlink(old_file)
      }
    }
    file_map[[old_file]] <- new_file
  }
  for (name in names(cache)) {
    old_file <- as.character(cache[[name]]$file)
    new_file <- unname(file_map[[old_file]])
    cache[[name]]$file <- new_file
    cache[[name]]$file_checksum <- unname(tools::md5sum(new_file))
    cache[[name]]$build_strategy <-
      "celltype_medium_corrected_python_corda2"
  }
  summary <- attr(cache, "summary")
  if (is.data.frame(summary) && nrow(summary)) {
    summary$file <- unname(file_map[as.character(summary$file)])
    summary$file_checksum <- unname(tools::md5sum(summary$file))
    summary$build_strategy <-
      "celltype_medium_corrected_python_corda2"
    summary$completion_stage <-
      "corrected_python_CORDA2_after_confidence_mapping"
    summary$completion_method <- "corda2"
    attr(cache, "summary") <- summary
  }
  attr(cache, "completion_method") <- "corda2"
  attr(cache, "fastcore_parallel_task") <- "not_applicable_to_corda2"
  cache
}

.rc_build_celltype_medium_union_gem_cache <- function(...) {
  context <- .rc_layer2_completion_context
  is_corda2 <- isTRUE(context$active) &&
    identical(context$model_completion, "corda") &&
    identical(
      as.character(context$corda_options$algorithm %||% ""),
      "resendislab_python_CORDA2_corrected_redundant_path_assessment"
    )
  cache <- .rc_build_celltype_medium_union_gem_cache_before_corda2_public(...)
  if (!is_corda2) return(cache)
  .rc_corda2_move_cache_files(cache)
}

.rc_regcompass_step_layer2_before_corda2_public <- rc_regcompass_step_layer2

rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  model_mode <- match.arg(model_mode)
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  original_params <- layer2_args$model_params %||% list()
  requested <- as.character(
    original_params$model_completion %||% "fastcore"
  )
  is_corda2 <- length(requested) == 1L && !is.na(requested) &&
    requested %in% c("corda2", "corda", "corda_like")
  translated <- original_params
  if (is_corda2) {
    translated$model_completion <- "corda"
    translated$corda_gamma <-
      original_params$corda2_penalty_factor %||%
      original_params$corda_penalty_factor %||%
      original_params$corda_gamma %||% 100
    translated$corda_kappa <-
      original_params$corda2_cost_increase %||%
      original_params$corda_cost_increase %||%
      original_params$corda_kappa %||% 1.01
    translated$corda_epsilon <-
      original_params$corda2_target_flux %||%
      original_params$corda_tflux %||%
      original_params$corda_epsilon %||% 1
    translated$corda_n <-
      original_params$corda2_redundancies %||%
      original_params$corda_n %||% 3L
    translated$corda_p <-
      original_params$corda2_support %||%
      original_params$corda_support %||%
      original_params$corda_p %||% 5L
    translated$corda_flux_tolerance <-
      original_params$corda2_flux_tolerance %||%
      original_params$corda_flux_tolerance %||% 1e-8
    translated[c(
      "corda2_penalty_factor", "corda_penalty_factor",
      "corda2_cost_increase", "corda_cost_increase",
      "corda2_target_flux", "corda_tflux",
      "corda2_redundancies", "corda2_support", "corda_support",
      "corda2_flux_tolerance"
    )] <- NULL
  }
  layer2_args$model_params <- translated
  answer <- .rc_regcompass_step_layer2_before_corda2_public(
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
  if (is_corda2) {
    options <- .rc_layer2_corda_options(original_params)
    answer$params$model_completion <- "corda2"
    answer$params$structural_completion <- "corda2"
    answer$params$structural_completion_algorithm <-
      "resendislab_python_CORDA2_corrected_redundant_path_assessment"
    answer$params$corda2_reference_commit <-
      options$python_reference_commit
    answer$params$corda2_redundancies <- options$redundancies
    answer$params$corda2_support <- options$support
    answer$params$corda2_penalty_factor <- options$penalty_factor
    answer$params$corda2_cost_increase <- options$cost_increase
    answer$params$corda2_target_flux <- options$target_flux
    answer$completion_contract$model_completion <- "corda2"
    answer$completion_contract$algorithm <-
      "resendislab_python_CORDA2_corrected_redundant_path_assessment"
    answer$completion_contract$python_reference <- list(
      repository = "resendislab/corda",
      commit = options$python_reference_commit,
      class = "CORDA2"
    )
    answer$completion_contract$confidence_levels <- c(
      absent = -1L, unknown = 0L, low = 1L, medium = 2L, high = 3L
    )
    answer$completion_contract$redundant_path_search <- list(
      maximum_paths = options$redundancies,
      cost_increase = options$cost_increase,
      absent_penalty_factor = options$penalty_factor,
      support_threshold = options$support,
      target_flux = options$target_flux
    )
    answer$completion_contract$intentional_corrections <-
      options$intentional_corrections
    answer$union_gem_policy <- paste(
      "one corrected Python-CORDA2 reconstruction per cell type and medium;",
      "five directional confidence levels, redundant-path cost increases,",
      "absent support counting, independent medium testing and final free",
      "reaction completion; only included HC core directions are scored"
    )
    answer$method <- paste(
      "microCOMPASS directional LP on cell-type-specific medium models",
      "reconstructed by corrected Python CORDA2"
    )
    rc_export_microcompass(answer, outdir)
    model_dir <- file.path(outdir, "03_models")
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      answer$completion_contract,
      file.path(model_dir, "model_completion_contract.rds")
    )
    saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  }
  answer
}
