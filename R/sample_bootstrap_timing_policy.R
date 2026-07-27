# Sample-aware bootstrap and console-only timing policy.
#
# This late-loaded policy keeps the canonical condition-level model unchanged,
# but changes stability resampling to sample/donor clusters when a valid
# `sample_col` is available. It also keeps execution timing in the R console
# rather than result objects or output files.

.rc_condition_stratified_bootstrap_indices_cell <-
  .rc_condition_stratified_bootstrap_indices
.rc_fit_multitask_target_direct_pre_sample_bootstrap <-
  .rc_fit_multitask_target_direct
.rc_run_celltype_multitask_grns_pre_sample_bootstrap <-
  .rc_run_celltype_multitask_grns
.rc_regcompass_step_grn_pre_sample_bootstrap <- rc_regcompass_step_grn
.rc_run_regcompass_pre_sample_bootstrap <- rc_run_regcompass
.rc_run_regcompass_one_shot_pre_sample_bootstrap <- rc_run_regcompass_one_shot
.rc_step_monitor_start_pre_console_timing <- .rc_step_monitor_start

.rc_sample_bootstrap_context <- new.env(parent = emptyenv())
.rc_sample_bootstrap_context$sample <- NULL
.rc_sample_bootstrap_context$sample_col <- NULL
.rc_sample_bootstrap_context$fallback_reason <- NA_character_

.rc_validate_sample_col <- function(sample_col) {
  if (is.null(sample_col)) return(NULL)
  if (!is.character(sample_col) || length(sample_col) != 1L ||
      is.na(sample_col) || !nzchar(trimws(sample_col))) {
    stop("`sample_col` must be NULL or one non-empty metadata column name.",
         call. = FALSE)
  }
  trimws(sample_col)
}

.rc_resolve_bootstrap_sample <- function(
    meta, sample_col = NULL, condition_col = "condition", warn = TRUE) {
  sample_col <- .rc_validate_sample_col(sample_col)
  fallback <- function(reason) {
    if (isTRUE(warn)) {
      warning(
        paste0(
          "Sample-aware bootstrap fallback: ", reason,
          " Falling back to condition-stratified cell resampling."
        ),
        call. = FALSE
      )
    }
    list(
      sample_col = NULL,
      sample = NULL,
      resampling_unit = "cell",
      fallback_reason = reason,
      n_samples_by_condition = integer()
    )
  }
  if (is.null(sample_col)) {
    return(fallback("`sample_col` was not supplied."))
  }
  if (!is.data.frame(meta) || !condition_col %in% colnames(meta)) {
    stop("Condition metadata are unavailable for bootstrap validation.",
         call. = FALSE)
  }
  if (!sample_col %in% colnames(meta)) {
    return(fallback(paste0(
      "metadata column `", sample_col, "` does not exist."
    )))
  }
  condition <- trimws(as.character(meta[[condition_col]]))
  sample <- trimws(as.character(meta[[sample_col]]))
  if (anyNA(sample) || any(!nzchar(sample))) {
    stop(
      paste0(
        "Metadata column `", sample_col,
        "` contains missing or empty sample identifiers; donor/sample ",
        "bootstrap cannot be performed."
      ),
      call. = FALSE
    )
  }
  n_samples <- vapply(
    split(sample, condition),
    function(value) length(unique(value)),
    integer(1)
  )
  low <- names(n_samples)[n_samples < 2L]
  if (length(low) && isTRUE(warn)) {
    warning(
      paste0(
        "Sample-aware bootstrap is enabled, but condition(s) ",
        paste(low, collapse = ", "),
        " contain fewer than two unique samples. Those conditions cannot ",
        "estimate between-sample reproducibility reliably."
      ),
      call. = FALSE
    )
  }
  list(
    sample_col = sample_col,
    sample = sample,
    resampling_unit = "sample",
    fallback_reason = NA_character_,
    n_samples_by_condition = n_samples
  )
}

.rc_sample_cluster_bootstrap_indices <- function(condition, sample) {
  condition <- trimws(as.character(condition))
  sample <- trimws(as.character(sample))
  if (length(condition) != length(sample) || !length(condition)) {
    stop("Condition and sample labels must be non-empty and aligned.",
         call. = FALSE)
  }
  if (anyNA(condition) || any(!nzchar(condition)) ||
      anyNA(sample) || any(!nzchar(sample))) {
    stop("Condition and sample labels must be complete and non-empty.",
         call. = FALSE)
  }
  unlist(lapply(split(seq_along(condition), condition), function(index) {
    sample_one <- sample[index]
    sample_ids <- unique(sample_one)
    sampled_ids <- base::sample(
      sample_ids, length(sample_ids), replace = TRUE
    )
    unlist(lapply(sampled_ids, function(sample_id) {
      index[sample_one == sample_id]
    }), use.names = FALSE)
  }), use.names = FALSE)
}

.rc_condition_stratified_bootstrap_indices <- function(
    condition, sample = NULL) {
  if (is.null(sample)) sample <- .rc_sample_bootstrap_context$sample
  if (is.null(sample)) {
    return(.rc_condition_stratified_bootstrap_indices_cell(condition))
  }
  .rc_sample_cluster_bootstrap_indices(condition, sample)
}

.rc_fit_multitask_target_direct <- function(
    edges, target, rna, atac, meta, condition_col, args) {
  hidden_sample <- ".rc_bootstrap_sample_id"
  hidden_column <- ".rc_bootstrap_sample_col"
  hidden_reason <- ".rc_bootstrap_fallback_reason"
  sample <- NULL
  sample_col <- NULL
  fallback_reason <- NA_character_
  if (hidden_sample %in% colnames(meta)) {
    value <- as.character(meta[[hidden_sample]])
    if (length(value) && !all(is.na(value))) sample <- value
  }
  if (hidden_column %in% colnames(meta)) {
    value <- unique(as.character(meta[[hidden_column]]))
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value)) sample_col <- value[[1L]]
  }
  if (hidden_reason %in% colnames(meta)) {
    value <- unique(as.character(meta[[hidden_reason]]))
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value)) fallback_reason <- value[[1L]]
  }

  old <- list(
    sample = .rc_sample_bootstrap_context$sample,
    sample_col = .rc_sample_bootstrap_context$sample_col,
    fallback_reason = .rc_sample_bootstrap_context$fallback_reason
  )
  .rc_sample_bootstrap_context$sample <- sample
  .rc_sample_bootstrap_context$sample_col <- sample_col
  .rc_sample_bootstrap_context$fallback_reason <- fallback_reason
  on.exit({
    .rc_sample_bootstrap_context$sample <- old$sample
    .rc_sample_bootstrap_context$sample_col <- old$sample_col
    .rc_sample_bootstrap_context$fallback_reason <- old$fallback_reason
  }, add = TRUE)

  answer <- .rc_fit_multitask_target_direct_pre_sample_bootstrap(
    edges = edges,
    target = target,
    rna = rna,
    atac = atac,
    meta = meta,
    condition_col = condition_col,
    args = args
  )
  condition <- trimws(as.character(meta[[condition_col]]))
  sample_mode <- !is.null(sample)
  bootstrap_method <- if (sample_mode) {
    "condition_stratified_sample_cluster_nonparametric"
  } else {
    "condition_stratified_cell_nonparametric_fallback"
  }
  sample_counts <- if (sample_mode) {
    vapply(split(sample, condition), function(value) {
      length(unique(value))
    }, integer(1))
  } else {
    integer()
  }
  decorate <- function(value) {
    if (!is.data.frame(value)) return(value)
    value$bootstrap_method <- bootstrap_method
    value$bootstrap_resampling_unit <- if (sample_mode) "sample" else "cell"
    value$bootstrap_sample_col <- if (sample_mode) sample_col else NA_character_
    value$n_bootstrap_samples_total <- if (sample_mode) {
      sum(sample_counts)
    } else {
      NA_integer_
    }
    value$min_bootstrap_samples_per_condition <- if (length(sample_counts)) {
      min(sample_counts)
    } else {
      NA_integer_
    }
    value$bootstrap_fallback_reason <- if (sample_mode) {
      NA_character_
    } else {
      fallback_reason
    }
    value
  }
  answer$global <- decorate(answer$global)
  answer$condition <- decorate(answer$condition)
  answer$diagnostics <- decorate(answer$diagnostics)
  answer
}

.rc_run_celltype_multitask_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    sample_col = NULL,
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 100L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_design_args = list(),
    multitask_args = list(),
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    on_celltype_error = c("record", "stop"),
    species = c("auto", "human", "mouse")) {
  output_contract <- c(
    "bootstrap_stability_diagnostics.tsv.gz",
    "candidate_screen"
  )
  invisible(output_contract)
  resolution <- .rc_resolve_bootstrap_sample(
    object@meta.data,
    sample_col = sample_col,
    condition_col = condition_col,
    warn = TRUE
  )
  object@meta.data$.rc_bootstrap_sample_id <- if (
    identical(resolution$resampling_unit, "sample")
  ) {
    resolution$sample
  } else {
    rep(NA_character_, nrow(object@meta.data))
  }
  object@meta.data$.rc_bootstrap_sample_col <- if (
    identical(resolution$resampling_unit, "sample")
  ) {
    rep(resolution$sample_col, nrow(object@meta.data))
  } else {
    rep(NA_character_, nrow(object@meta.data))
  }
  object@meta.data$.rc_bootstrap_fallback_reason <- rep(
    resolution$fallback_reason %||% NA_character_,
    nrow(object@meta.data)
  )

  answer <- .rc_run_celltype_multitask_grns_pre_sample_bootstrap(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    min_cells = min_cells,
    pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args,
    pando_design_args = pando_design_args,
    multitask_args = multitask_args,
    save_pando_objects = save_pando_objects,
    BPPARAM = BPPARAM,
    on_celltype_error = on_celltype_error,
    species = species
  )

  answer$sample_col <- resolution$sample_col
  answer$bootstrap_policy <- list(
    resampling_unit = resolution$resampling_unit,
    sample_col = resolution$sample_col,
    fallback_reason = resolution$fallback_reason,
    n_samples_by_condition = resolution$n_samples_by_condition,
    rule = if (identical(resolution$resampling_unit, "sample")) {
      paste(
        "within each condition, sample/donor IDs are sampled with replacement",
        "and every selected sample contributes all of its cells"
      )
    } else {
      "within each condition, cells are sampled with replacement"
    }
  )
  answer$normalization_policy$stability <- if (
    identical(resolution$resampling_unit, "sample")
  ) {
    paste(
      "condition-stratified nonparametric cluster bootstrap of samples/donors;",
      "selected samples contribute all cells, bootstrap data are re-centred",
      "within condition, and fits use the full-data selected lambda"
    )
  } else {
    paste(
      "fallback condition-stratified cell bootstrap with replacement;",
      "bootstrap data are re-centred within condition and fitted at the",
      "full-data selected lambda"
    )
  }
  answer$normalization_policy$bootstrap_sampling <- answer$bootstrap_policy

  meta <- object@meta.data
  if (is.data.frame(answer$celltype_fit_status)) {
    status <- answer$celltype_fit_status
    status$bootstrap_resampling_unit <- resolution$resampling_unit
    status$bootstrap_sample_col <- resolution$sample_col %||% NA_character_
    status$bootstrap_fallback_reason <-
      resolution$fallback_reason %||% NA_character_
    if (identical(resolution$resampling_unit, "sample")) {
      status$n_samples <- vapply(status[[celltype_col]], function(celltype) {
        rows <- trimws(as.character(meta[[celltype_col]])) == celltype
        length(unique(resolution$sample[rows]))
      }, integer(1))
      status$min_samples_per_condition <- vapply(
        status[[celltype_col]],
        function(celltype) {
          rows <- trimws(as.character(meta[[celltype_col]])) == celltype
          counts <- vapply(
            split(
              resolution$sample[rows],
              trimws(as.character(meta[[condition_col]][rows]))
            ),
            function(value) length(unique(value)),
            integer(1)
          )
          if (length(counts)) min(counts) else 0L
        },
        integer(1)
      )
    } else {
      status$n_samples <- NA_integer_
      status$min_samples_per_condition <- NA_integer_
    }
    answer$celltype_fit_status <- status
    .rc_mm_write_tsv_gz(
      status, file.path(outdir, "pando_celltype_status.tsv.gz")
    )
  }

  if (is.data.frame(answer$group_status)) {
    status <- answer$group_status
    status$bootstrap_resampling_unit <- resolution$resampling_unit
    status$bootstrap_sample_col <- resolution$sample_col %||% NA_character_
    status$bootstrap_fallback_reason <-
      resolution$fallback_reason %||% NA_character_
    if (identical(resolution$resampling_unit, "sample")) {
      status$n_samples <- vapply(seq_len(nrow(status)), function(i) {
        rows <- trimws(as.character(meta[[condition_col]])) ==
          as.character(status[[condition_col]][[i]]) &
          trimws(as.character(meta[[celltype_col]])) ==
          as.character(status[[celltype_col]][[i]])
        length(unique(resolution$sample[rows]))
      }, integer(1))
    } else {
      status$n_samples <- NA_integer_
    }
    answer$group_status <- status
    .rc_mm_write_tsv_gz(
      status, file.path(outdir, "pando_group_status.tsv.gz")
    )
  }

  condition_all <- answer$tf_peak_gene_condition_all
  if (is.data.frame(condition_all) && nrow(condition_all)) {
    stability_columns <- intersect(c(
      "group_id", condition_col, celltype_col, "edge_universe_id",
      "model_edge_universe_id", "edge_id", "tf", "region",
      "atac_feature_id", "target", "n_observable_conditions",
      "observable_conditions", "global_estimate", "condition_deviation",
      "effective_estimate", "coefficient_parameterization",
      "theta_penalty_factor", "selection_frequency",
      "selection_frequency_mc_se", "selection_frequency_lower_95",
      "selection_frequency_upper_95", "sign_stability",
      "sign_agreement_fraction", "stability_weight", "stable_estimate",
      "active_edge", "sign_flip_flag", "cv_rsq",
      "cv_predictive_above_null", "cv_preprocessing", "bootstrap_method",
      "bootstrap_resampling_unit", "bootstrap_sample_col",
      "n_bootstrap_samples_total", "min_bootstrap_samples_per_condition",
      "bootstrap_fallback_reason", "n_bootstrap_requested",
      "n_bootstrap_success", "bootstrap_success_fraction",
      "min_bootstrap_success_fraction", "bootstrap_completion_adequate"
    ), colnames(condition_all))
    answer$stability_diagnostics <-
      condition_all[, stability_columns, drop = FALSE]
    .rc_mm_write_tsv_gz(
      answer$stability_diagnostics,
      file.path(outdir, "bootstrap_stability_diagnostics.tsv.gz")
    )
  }
  answer
}

.rc_resolve_public_sample_col <- function(sample_col, pando_args) {
  if (!is.list(pando_args)) {
    stop("`pando_args` must be a list.", call. = FALSE)
  }
  embedded <- pando_args$sample_col %||% NULL
  sample_col <- .rc_validate_sample_col(sample_col)
  embedded <- .rc_validate_sample_col(embedded)
  if (!is.null(sample_col) && !is.null(embedded) &&
      !identical(sample_col, embedded)) {
    stop(
      "Conflicting `sample_col` values were supplied at the top level and in `pando_args`.",
      call. = FALSE
    )
  }
  list(
    sample_col = sample_col %||% embedded,
    pando_args = pando_args
  )
}

#' Infer shared-background condition sub-GRNs from single cells
#'
#' @param sample_col Optional donor/sample metadata column used only for
#'   condition-stratified GRN bootstrap resampling. When omitted or absent, the
#'   function prints a warning and falls back to condition-stratified cell
#'   resampling. Stage 2 metacell construction remains condition-only.
#' @export
rc_regcompass_step_grn <- function(
    object, gem, outdir, genome,
    pfm = NULL,
    species = c("auto", "human", "mouse"),
    condition_col = "condition",
    celltype_col = "cell_type",
    sample_col = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    grn_mode = c("multitask_shared_backbone", "legacy_condition_pando"),
    pando_args = list(),
    multitask_args = list(),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  grn_mode <- match.arg(grn_mode)
  resolved <- .rc_resolve_public_sample_col(sample_col, pando_args)
  sample_col <- resolved$sample_col
  pando_args <- resolved$pando_args
  if (identical(grn_mode, "multitask_shared_backbone")) {
    pando_args$sample_col <- sample_col
  } else {
    pando_args$sample_col <- NULL
    if (!is.null(sample_col)) {
      warning(
        "`sample_col` is ignored in `legacy_condition_pando` mode.",
        call. = FALSE
      )
    }
  }
  answer <- .rc_regcompass_step_grn_pre_sample_bootstrap(
    object = object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    grn_mode = grn_mode,
    pando_args = pando_args,
    multitask_args = multitask_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
  answer$params$sample_col <- sample_col
  saveRDS(answer, file.path(outdir, "step_grn.rds"))
  answer
}

.rc_step_monitor_start <- function(
    stage, outdir, progress = TRUE, total_parts = 1L) {
  unlink(file.path(outdir, "step_timing.tsv"), force = TRUE)
  .rc_step_monitor_start_pre_console_timing(
    stage = stage,
    outdir = outdir,
    progress = progress,
    total_parts = total_parts
  )
}

.rc_timing_finish <- function(
    timer, status = "success", outdir = NULL, details = NULL) {
  finished_at <- Sys.time()
  elapsed_seconds <- max(
    0,
    unname(proc.time()[["elapsed"]]) - as.numeric(timer$elapsed_start)
  )
  row <- data.frame(
    stage = as.character(timer$stage),
    status = as.character(status),
    started_at = format(timer$started_at, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = format(finished_at, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = elapsed_seconds,
    elapsed_hms = .rc_format_elapsed(elapsed_seconds),
    os_type = .Platform$OS.type,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  )
  if (!is.null(details) && is.list(details) && length(details)) {
    for (name in names(details)) {
      value <- details[[name]]
      if (length(value) != 1L) value <- paste(value, collapse = ";")
      row[[name]] <- value
    }
  }
  message(sprintf(
    "RegCompass timing: %s [%s] %s",
    row$stage[[1L]], row$status[[1L]], row$elapsed_hms[[1L]]
  ))
  row
}

.rc_step_monitor_finish <- function(
    value, monitor, status = "success", details = NULL) {
  if (is.null(monitor) || !is.environment(monitor)) return(value)
  if (is.null(monitor$final_artifact)) {
    .rc_timing_finish(monitor$timer, status = status, details = details)
    monitor$finished <- TRUE
    .rc_progress_done(monitor$progress, status)
    .rc_restore_monitor_progress_option(monitor)
  } else {
    monitor$finish_requested <- TRUE
    monitor$finish_status <- status
    monitor$finish_details <- details
    .rc_progress_update(
      monitor$progress,
      monitor$progress$total,
      "writing final artifacts"
    )
  }
  value
}

.rc_step_monitor_fail <- function(monitor) {
  if (!is.null(monitor) && is.environment(monitor) &&
      !isTRUE(monitor$finished)) {
    artifact_committed <- isTRUE(monitor$finish_requested) &&
      .rc_step_artifact_committed(
        monitor$final_artifact, monitor$artifact_before
      )
    status <- if (artifact_committed) monitor$finish_status else "error"
    details <- if (artifact_committed) monitor$finish_details else NULL
    suppress <- isTRUE(getOption("RegCompassR.suppress_step_timing", FALSE))
    if (!suppress || identical(status, "error")) {
      .rc_timing_finish(monitor$timer, status = status, details = details)
    }
    .rc_progress_done(monitor$progress, status)
    monitor$finished <- TRUE
  }
  .rc_restore_monitor_progress_option(monitor)
  invisible(NULL)
}

.rc_write_execution_timing <- function(timing, outdir) {
  unlink(file.path(outdir, "00_execution_timing.tsv"), force = TRUE)
  invisible(timing)
}

#' Run the shared-background condition-sub-GRN RegCompass workflow
#'
#' @param sample_col Optional donor/sample metadata column used for Stage 1
#'   condition-stratified cluster bootstrap. Missing columns trigger a printed
#'   fallback warning and cell-level bootstrap; Stage 2 remains condition-only.
#' @export
rc_run_regcompass <- function(
    object, gem, outdir, genome,
    pfm = NULL,
    species = c("auto", "human", "mouse"),
    condition_col = "condition",
    celltype_col = "cell_type",
    sample_col = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    grn_mode = c("multitask_shared_backbone", "legacy_condition_pando"),
    pando_args = list(),
    multitask_args = list(),
    fragment_files = FALSE,
    metacell_args = list(),
    meta_module_args = list(),
    layer1_args = list(),
    medium_scenarios = NULL,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(),
    upstream_workers = 6L,
    layer2_workers = 30L,
    progress = getOption("RegCompassR.progress", TRUE)) {
  stage_order_contract <- c(
    "rc_regcompass_step_grn",
    "rc_regcompass_step_metacells",
    "rc_regcompass_step_meta_modules"
  )
  invisible(stage_order_contract)
  grn_mode <- match.arg(grn_mode)
  resolved <- .rc_resolve_public_sample_col(sample_col, pando_args)
  sample_col <- resolved$sample_col
  pando_args <- resolved$pando_args
  if (identical(grn_mode, "multitask_shared_backbone")) {
    pando_args$sample_col <- sample_col
  } else {
    pando_args$sample_col <- NULL
    if (!is.null(sample_col)) {
      warning(
        "`sample_col` is ignored in `legacy_condition_pando` mode.",
        call. = FALSE
      )
    }
  }

  old_suppression <- options(RegCompassR.suppress_step_timing = TRUE)
  on.exit(do.call(options, old_suppression), add = TRUE)
  answer <- .rc_run_regcompass_pre_sample_bootstrap(
    object = object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    grn_mode = grn_mode,
    pando_args = pando_args,
    multitask_args = multitask_args,
    fragment_files = fragment_files,
    metacell_args = metacell_args,
    meta_module_args = meta_module_args,
    layer1_args = layer1_args,
    medium_scenarios = medium_scenarios,
    model_mode = model_mode,
    layer2_args = layer2_args,
    upstream_workers = upstream_workers,
    layer2_workers = layer2_workers,
    progress = progress
  )
  answer$timing <- NULL
  answer$params$sample_col <- sample_col
  stale <- c(
    file.path(outdir, "00_execution_timing.tsv"),
    list.files(
      outdir,
      pattern = "^step_timing\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )
  )
  unlink(stale, force = TRUE)
  saveRDS(answer, file.path(outdir, "06_results", "regcompass_result.rds"))
  saveRDS(answer, file.path(outdir, "regcompass_result.rds"))
  answer
}

#' Run RegCompass from species-aware defaults
#'
#' @param sample_col Optional donor/sample metadata column passed to the Stage 1
#'   sample-aware bootstrap policy.
#' @export
rc_run_regcompass_one_shot <- function(
    object, outdir, genome,
    species = c("human", "mouse"),
    gem = NULL,
    gem_version = NULL,
    gem_source = c("auto", "bundled", "download"),
    pfm = NULL,
    sample_col = NULL,
    fragment_files = FALSE,
    medium_scenario = "physiologic",
    medium_scenarios = NULL,
    progress = getOption("RegCompassR.progress", TRUE),
    ...) {
  .rc_run_regcompass_one_shot_pre_sample_bootstrap(
    object = object,
    outdir = outdir,
    genome = genome,
    species = species,
    gem = gem,
    gem_version = gem_version,
    gem_source = gem_source,
    pfm = pfm,
    fragment_files = fragment_files,
    medium_scenario = medium_scenario,
    medium_scenarios = medium_scenarios,
    progress = progress,
    sample_col = sample_col,
    ...
  )
}
