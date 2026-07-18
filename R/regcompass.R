.rc_resolve_workflow_or_method <- function(
    layer1_args = list(), strict_biological_defaults = TRUE) {
  if (!is.list(layer1_args)) stop("`layer1_args` must be a list.", call. = FALSE)
  if (!is.null(layer1_args$or_method)) {
    return(match.arg(
      as.character(layer1_args$or_method),
      c("max", "sum_sqrtK", "prob_or", "sum")
    ))
  }
  if (isTRUE(strict_biological_defaults)) "max" else "sum_sqrtK"
}

.rc_workflow_or_method_source <- function(
    layer1_args = list(), strict_biological_defaults = TRUE) {
  if (!is.null(layer1_args$or_method)) return("explicit_layer1_args")
  if (isTRUE(strict_biological_defaults)) {
    "strict_biological_default"
  } else {
    "legacy_sensitivity_default"
  }
}

#' Run the canonical RegCompass workflow
#'
#' Builds strict condition-by-sample-by-cell-type metacells, estimates signed
#' Pando regulatory evidence, constructs either a shared full GEM or a
#' FASTCORE-completed shared meta-module GEM, and runs directional
#' minimum-evidence-discordance LPs. Human-GEM and Mouse-GEM are supported
#' without cross-species gene conversion.
#'
#' @param object A Seurat object containing RNA and ATAC assays.
#' @param gem A species GEM prepared by `rc_prepare_gem()`.
#' @param outdir Persistent output directory.
#' @param pfm Motif position-frequency matrices passed to Pando.
#' @param genome Genome object matching `species` and the ATAC coordinates.
#' @param fragment_files Fragment-file manifest/path(s), or `FALSE` to skip
#'   fragment aggregation and use ATAC peak raw counts from `object` directly.
#' @param species `"auto"`, `"human"`, or `"mouse"`.
#' @param sample_col,condition_col,celltype_col Metadata columns.
#' @param rna_assay,atac_assay Assay names.
#' @param model_mode `"meta_module_gem"` or `"full_gem"`.
#' @param medium_scenarios Shared medium table. The default is the
#'   species-matched literature-backed physiological environment.
#' @param metacell_args,layer1_args,pando_args,layer2_args Named argument lists.
#' @param upstream_workers,layer2_workers Worker counts.
#' @param parallel_backend Parallel backend.
#' @param strict_biological_defaults Use `promiscuity=none`, `AND=min`, and
#'   `OR=max` unless explicitly overridden.
#' @param inference_unit Primary scoring unit. Sample-by-cell-type is recommended
#'   for biological inference; metacells are exploratory observations.
#' @return A RegCompass result list; canonical RDS files and model cache files are
#'   written below `outdir`.
#' @export
rc_run_regcompass <- function(
    object, gem, outdir, pfm, genome,
    fragment_files = NULL,
    sample_col = "sample_id",
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
    strict_biological_defaults = TRUE,
    inference_unit = c("sample_celltype", "metacell"),
    species = c("auto", "human", "mouse")) {
  species <- .rc_infer_gem_species(gem, species)
  model_mode <- match.arg(model_mode)
  parallel_backend <- match.arg(parallel_backend)
  inference_unit <- match.arg(inference_unit)
  if (!is.logical(strict_biological_defaults) ||
      length(strict_biological_defaults) != 1L ||
      is.na(strict_biological_defaults)) {
    stop("`strict_biological_defaults` must be TRUE or FALSE.", call. = FALSE)
  }
  bundles <- list(
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args,
    layer2_args = layer2_args
  )
  invalid <- names(bundles)[!vapply(bundles, is.list, logical(1))]
  if (length(invalid)) {
    stop(
      "Workflow argument bundles must be lists: ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }
  or_method <- .rc_resolve_workflow_or_method(
    layer1_args, strict_biological_defaults
  )
  or_method_source <- .rc_workflow_or_method_source(
    layer1_args, strict_biological_defaults
  )
  if (isTRUE(strict_biological_defaults)) {
    layer1_args$promiscuity_mode <- layer1_args$promiscuity_mode %||% "none"
    layer1_args$and_method <- layer1_args$and_method %||% "min"
  } else {
    layer1_args$promiscuity_mode <- layer1_args$promiscuity_mode %||% "sqrt"
    layer1_args$and_method <- layer1_args$and_method %||% "boltzmann"
  }
  layer1_args$or_method <- or_method
  layer1_args$promiscuity_mode <- match.arg(
    layer1_args$promiscuity_mode,
    c("none", "sqrt", "linear")
  )
  layer1_args$and_method <- match.arg(
    layer1_args$and_method,
    c("min", "boltzmann", "mean")
  )
  layer1_args$or_method <- match.arg(
    layer1_args$or_method,
    c("max", "sum_sqrtK", "prob_or", "sum")
  )
  if (is.null(medium_scenarios)) {
    medium_scenarios <- rc_make_medium_scenarios(
      gem,
      scenario = "physiologic",
      species = species
    )
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(outdir, "04_model_cache", model_mode)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  layer2_args$model_params <- layer2_args$model_params %||% list()
  if (!is.list(layer2_args$model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  supplied_cache <- layer2_args$model_params$cache_dir %||% cache_dir
  if (!identical(normalizePath(supplied_cache, mustWork = FALSE),
                 normalizePath(cache_dir, mustWork = FALSE))) {
    stop(
      "The integrated workflow owns `layer2_args$model_params$cache_dir`; use the persistent path below `outdir`.",
      call. = FALSE
    )
  }
  layer2_args$model_params$cache_dir <- cache_dir
  previous <- options(RegCompassR.inference_unit = inference_unit)
  on.exit(options(previous), add = TRUE)
  answer <- .rc_run_regcompass_engine(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    model_mode = model_mode,
    medium_scenarios = medium_scenarios,
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args,
    layer2_args = layer2_args,
    upstream_workers = upstream_workers,
    layer2_workers = layer2_workers,
    parallel_backend = parallel_backend
  )
  answer$schema_version <- "regcompass_v5_species_canonical"
  answer$species <- species
  answer$model_info <- gem$model_info
  answer$layer1$C_abs <- answer$layer1$C_abs %||% answer$layer1$C_rel
  answer$layer1$gene_activity_absolute <- answer$layer1$global_gene_score
  answer$layer1$gene_state_relative <- .rc_weighted_gene_score(
    answer$layer1$rna_metacell_logcpm,
    answer$layer1$sample_balance_weights,
    mode = "relative"
  )
  answer$layer1$capacity_calibration_scope <-
    "zero_preserving_absolute_activity_with_weighted_stratified_q95_diagnostics"
  answer$layer1$reaction_confidence_source <-
    "pando_signed_tf_peak_gene_regulatory_support"
  answer$params$species <- species
  answer$params$model_source <- gem$model_info$source %||% NA_character_
  answer$params$model_version <- gem$model_info$version %||% NA_character_
  answer$params$strict_biological_defaults <- strict_biological_defaults
  answer$params$gpr_promiscuity_mode <- layer1_args$promiscuity_mode
  answer$params$gpr_and_method <- layer1_args$and_method
  answer$params$gpr_or_method <- or_method
  answer$params$gpr_or_method_source <- or_method_source
  answer$params$q95_role <- "diagnostic_only"
  answer$params$sample_balance_role <-
    answer$layer1$sample_balance_role %||%
    "q95_and_relative_state_diagnostics_only"
  answer$params$sample_balance_estimand <-
    answer$layer1$sample_balance_estimand %||%
    "equal biological-sample mass globally and within each Q95 stratum"
  answer$params$q95_n0 <- answer$layer1$calibration_params$q95_n0
  answer$params$q95_stratum_col <-
    answer$layer1$calibration_params$q95_stratum_col
  answer$params$model_cache_dir <- cache_dir
  answer$params$medium_policy <-
    "species_matched_literature_catalog_with_original_gem_bound_intersection"
  answer$params$inference_unit <- inference_unit
  answer$params$penalty_unit <- inference_unit
  saveRDS(gem$model_info, file.path(outdir, "00_model_info.rds"))
  saveRDS(medium_scenarios, file.path(outdir, "01_medium_scenarios.rds"))
  saveRDS(answer, file.path(outdir, "regcompass_global_metacell_result.rds"))
  saveRDS(answer, file.path(outdir, "regcompass_result.rds"))
  answer
}

#' Run the RegCompass workflow
#'
#' Builds strict-stratum RNA+ATAC metacells, runs Pando and local meta-module
#' completion, computes sample-balanced global and within-stratum Q95
#' diagnostics without rescaling absolute metacell activity, constructs one
#' shared GEM, and performs directional sample-by-cell-type inference.
.rc_run_regcompass_engine <- function(object, gem, outdir, pfm, genome,
                               fragment_files = NULL,
                               sample_col = "sample_id",
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
                               parallel_backend = c("auto", "serial", "snow", "multicore")) {
  model_mode <- match.arg(model_mode)
  parallel_backend <- match.arg(parallel_backend)
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  if (is.null(gem$gpr_table)) stop("`gem` must contain `gpr_table`.", call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  group_cols <- c(condition_col, sample_col, celltype_col)
  meta <- object@meta.data
  missing_cols <- setdiff(group_cols, colnames(meta))
  if (length(missing_cols)) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  all_samples <- unique(as.character(meta[[sample_col]]))
  if (anyNA(all_samples) || any(!nzchar(all_samples))) stop("Sample IDs must be non-missing and non-empty.", call. = FALSE)

  minimum_cells <- as.integer(metacell_args$min_cells_per_stratum %||% 100L)
  group_ids <- rc_make_stratum_id(meta, group_cols)
  counts <- table(group_ids)
  retained_ids <- names(counts[counts >= minimum_cells])
  if (!length(retained_ids)) stop("No strict strata satisfy `min_cells_per_stratum`.", call. = FALSE)
  retained_samples <- unique(as.character(meta[[sample_col]][group_ids %in% retained_ids]))
  samples_without_retained_strata <- setdiff(all_samples, retained_samples)
  if (length(samples_without_retained_strata)) {
    stop("Every biological sample must contribute at least one retained strict stratum. Missing samples: ",
         paste(samples_without_retained_strata, collapse = ", "), call. = FALSE)
  }
  excluded_ids <- setdiff(names(counts), retained_ids)
  status_dir <- file.path(outdir, "00_strata")
  dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(
    data.frame(group_id = names(counts), n_cells = as.integer(counts),
               retained = names(counts) %in% retained_ids, stringsAsFactors = FALSE),
    file.path(status_dir, "stratum_filter.tsv.gz")
  )

  upstream_dir <- file.path(outdir, "01_stratum_workflows")
  upstream_param <- .rc_phase_bpparam(upstream_workers, parallel_backend)
  run_one <- function(group_id) {
    tryCatch(
      .rc_run_regcompass_stratum(
        object = object,
        group_id = group_id,
        group_cols = group_cols,
        gem = gem,
        outdir = upstream_dir,
        pfm = pfm,
        genome = genome,
        fragment_files = fragment_files,
        sample_col = sample_col,
        condition_col = condition_col,
        celltype_col = celltype_col,
        rna_assay = rna_assay,
        atac_assay = atac_assay,
        metacell_args = metacell_args,
        layer1_args = layer1_args,
        pando_args = pando_args
      ),
      error = function(error) {
        list(group_id = group_id, status = "failed", artifact_file = NA_character_,
             n_cells = as.integer(counts[[group_id]]), n_metacells = NA_integer_,
             error_class = class(error)[[1L]], error_message = conditionMessage(error))
      }
    )
  }
  upstream_results <- tryCatch(
    rc_parallel_lapply(retained_ids, run_one, BPPARAM = upstream_param),
    finally = .rc_release_bpparam(upstream_param)
  )
  rm(upstream_param)
  invisible(gc(verbose = FALSE))

  retained_status <- .rc_bind_frames_fill(lapply(upstream_results, as.data.frame))
  upstream_status <- retained_status
  if (length(excluded_ids)) {
    excluded_status <- data.frame(
      group_id = excluded_ids, status = "skipped_too_few_cells", artifact_file = NA_character_,
      n_cells = as.integer(counts[excluded_ids]), n_metacells = NA_integer_,
      error_class = NA_character_, error_message = NA_character_, stringsAsFactors = FALSE
    )
    upstream_status <- .rc_bind_frames_fill(list(retained_status, excluded_status))
  }
  .rc_write_tsv_gz(upstream_status, file.path(status_dir, "stratum_workflow_status.tsv.gz"))

  completed_ids <- as.character(retained_status$group_id[retained_status$status == "ok"])
  skipped_metacell_ids <- as.character(retained_status$group_id[retained_status$status == "skipped_too_few_metacells"])
  failed_ids <- as.character(retained_status$group_id[retained_status$status == "failed"])
  artifact_files <- retained_status$artifact_file[match(completed_ids, retained_status$group_id)]
  missing_artifacts <- completed_ids[is.na(artifact_files) | !file.exists(artifact_files)]
  completed_samples <- unique(as.character(meta[[sample_col]][group_ids %in% completed_ids]))
  samples_without_completed_strata <- setdiff(all_samples, completed_samples)
  if (length(failed_ids) || length(missing_artifacts) || !length(completed_ids) || length(samples_without_completed_strata)) {
    barrier <- data.frame(
      stage = "upstream_complete_barrier",
      passed = FALSE,
      n_required_samples = length(all_samples),
      n_retained_strata = length(retained_ids),
      n_completed_strata = length(completed_ids),
      n_skipped_too_few_metacells = length(skipped_metacell_ids),
      failed_strata = paste(failed_ids, collapse = ";"),
      skipped_too_few_metacells = paste(skipped_metacell_ids, collapse = ";"),
      missing_artifacts = paste(missing_artifacts, collapse = ";"),
      missing_samples = paste(samples_without_completed_strata, collapse = ";"),
      stringsAsFactors = FALSE
    )
    .rc_write_tsv_gz(barrier, file.path(status_dir, "upstream_barrier.tsv.gz"))
    stop(
      "Global recalibration and union-GEM construction were blocked because one or more retained strata failed, no analyzable strata remained, or a biological sample had no analyzable stratum. Inspect stratum_workflow_status.tsv.gz and upstream_barrier.tsv.gz.",
      call. = FALSE
    )
  }
  barrier <- data.frame(
    stage = "upstream_complete_barrier",
    passed = TRUE,
    n_required_samples = length(all_samples),
    n_retained_strata = length(retained_ids),
    n_completed_strata = length(completed_ids),
    n_skipped_too_few_metacells = length(skipped_metacell_ids),
    failed_strata = "",
    skipped_too_few_metacells = paste(skipped_metacell_ids, collapse = ";"),
    missing_artifacts = "",
    missing_samples = "",
    stringsAsFactors = FALSE
  )
  .rc_write_tsv_gz(barrier, file.path(status_dir, "upstream_barrier.tsv.gz"))

  artifacts <- lapply(artifact_files, readRDS)
  artifact_group_ids <- vapply(artifacts, function(x) as.character(x$group_id), character(1))
  if (!identical(artifact_group_ids, completed_ids)) {
    stop("Upstream artifacts are incomplete or out of analyzable strict-stratum order.", call. = FALSE)
  }
  artifact_contract <- vapply(
    artifacts,
    .rc_validate_stratum_artifact_contract,
    character(1)
  )
  invalid_artifact <- artifact_contract != "ok"
  if (any(invalid_artifact)) {
    stop(
      "One or more upstream artifacts failed input/output contract validation: ",
      paste(unique(artifact_contract[invalid_artifact]), collapse = ", "),
      call. = FALSE
    )
  }

  single_cell_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  grn_meta_modules <- .rc_merge_stratum_meta_modules(artifacts)
  layer1 <- .rc_merge_stratum_layer1(
    artifacts = artifacts,
    gem = gem,
    single_cell_genes = single_cell_genes,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  missing_global_samples <- setdiff(all_samples, unique(as.character(layer1$unit_meta[[sample_col]])))
  if (length(missing_global_samples)) {
    stop("Global Layer 1 is missing biological samples after the upstream barrier: ",
         paste(missing_global_samples, collapse = ", "), call. = FALSE)
  }
  saveRDS(layer1, file.path(outdir, "02_global_layer1.rds"))
  saveRDS(grn_meta_modules, file.path(outdir, "03_global_meta_modules.rds"))
  rm(artifacts)
  invisible(gc(verbose = FALSE))

  reserved_layer2 <- intersect(
    names(layer2_args),
    c("layer1", "gem", "mode", "unit", "reaction_membership", "core_reactions",
      "target_reactions", "medium_scenarios", "sample_col", "condition_col",
      "celltype_col", "BPPARAM", "parallel")
  )
  if (length(reserved_layer2)) {
    stop("`layer2_args` cannot override integrated workflow fields: ", paste(reserved_layer2, collapse = ", "), call. = FALSE)
  }
  layer2_param <- .rc_phase_bpparam(layer2_workers, parallel_backend)
  layer2_defaults <- list(
    layer1 = layer1,
    gem = gem,
    target_reactions = grn_meta_modules$global_core_reactions$reaction_id,
    medium_scenarios = medium_scenarios,
    mode = model_mode,
    reaction_membership = if (identical(model_mode, "meta_module_gem")) grn_meta_modules$global_reaction_membership else NULL,
    core_reactions = if (identical(model_mode, "meta_module_gem")) grn_meta_modules$global_core_reactions else NULL,
    unit = "metacell",
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    parallel = TRUE,
    BPPARAM = layer2_param
  )
  layer2_defaults[names(layer2_args)] <- NULL
  microcompass <- tryCatch(
    do.call(rc_run_microcompass, c(layer2_defaults, layer2_args)),
    finally = .rc_release_bpparam(layer2_param)
  )
  rm(layer2_param)
  invisible(gc(verbose = FALSE))

  result <- list(
    schema_version = "regcompass_v3_global_metacell",
    model_mode = model_mode,
    strict_stratum_cols = group_cols,
    upstream_status = upstream_status,
    upstream_barrier = barrier,
    layer1 = layer1,
    grn_meta_modules = grn_meta_modules,
    microcompass = microcompass,
    params = list(
      shared_gem = TRUE,
      penalty_unit = "metacell",
      capacity_calibration_scope = layer1$capacity_calibration_scope,
      sample_balanced_gene_score = isTRUE(layer1$calibration_params$sample_balance),
      sample_balanced_q95 = isTRUE(layer1$calibration_params$sample_balance),
      expression_batch_correction = layer1$calibration_params$expression_batch_correction,
      local_fastcore_before_global_union = identical(
        grn_meta_modules$global_union_source,
        "deduplicated_local_fastcore_completed_meta_modules"
      ),
      global_fastcore_repair = "conditional_on_global_core_directional_incompleteness",
      upstream_parallel_unit = "condition_sample_celltype",
      global_stage_requires_all_retained_strata = TRUE,
      global_stage_requires_all_samples = TRUE,
      layer2_parallel_unit = "shared_model_by_metacell",
      parallel_backend = parallel_backend,
      upstream_workers = upstream_workers,
      layer2_workers = layer2_workers
    )
  )
  saveRDS(result, file.path(outdir, "regcompass_global_metacell_result.rds"))
  result
}
