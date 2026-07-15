#' Run the integrated RegCompass workflow
#'
#' Each retained condition-by-sample-by-cell-type stratum is processed by one
#' upstream worker from metacell construction through Pando and meta-module
#' inference. The resulting reaction envelopes are unioned into one shared GEM,
#' and a fresh worker pool evaluates every metacell with its own penalty vector.
#' @export
rc_run_regcompass <- function(object, gem, outdir, pfm, genome,
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
                               on_stratum_error = c("record", "stop")) {
  model_mode <- match.arg(model_mode)
  parallel_backend <- match.arg(parallel_backend)
  on_stratum_error <- match.arg(on_stratum_error)
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  if (is.null(gem$gpr_table)) stop("`gem` must contain `gpr_table`.", call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  group_cols <- c(condition_col, sample_col, celltype_col)
  meta <- object@meta.data
  missing_cols <- setdiff(group_cols, colnames(meta))
  if (length(missing_cols)) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  minimum_cells <- as.integer(metacell_args$min_cells_per_stratum %||% 100L)
  group_ids <- rc_make_stratum_id(meta, group_cols)
  counts <- table(group_ids)
  retained_ids <- names(counts[counts >= minimum_cells])
  if (!length(retained_ids)) stop("No strict strata satisfy `min_cells_per_stratum`.", call. = FALSE)
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
        if (identical(on_stratum_error, "stop")) stop(error)
        list(group_id = group_id, status = "failed", artifact_file = NA_character_,
             n_cells = as.integer(counts[[group_id]]), n_metacells = NA_integer_,
             error_class = class(error)[[1L]], error_message = conditionMessage(error))
      }
    )
  }
  upstream_results <- rc_parallel_lapply(retained_ids, run_one, BPPARAM = upstream_param)
  .rc_release_bpparam(upstream_param)
  rm(upstream_param)
  invisible(gc(verbose = FALSE))

  upstream_status <- .rc_bind_frames_fill(lapply(upstream_results, as.data.frame))
  if (length(excluded_ids)) {
    excluded_status <- data.frame(
      group_id = excluded_ids, status = "skipped_too_few_cells", artifact_file = NA_character_,
      n_cells = as.integer(counts[excluded_ids]), n_metacells = NA_integer_,
      error_class = NA_character_, error_message = NA_character_, stringsAsFactors = FALSE
    )
    upstream_status <- .rc_bind_frames_fill(list(upstream_status, excluded_status))
  }
  .rc_write_tsv_gz(upstream_status, file.path(status_dir, "stratum_workflow_status.tsv.gz"))
  artifact_files <- upstream_status$artifact_file[upstream_status$status == "ok"]
  artifact_files <- artifact_files[!is.na(artifact_files) & file.exists(artifact_files)]
  if (!length(artifact_files)) stop("All strict-stratum workflows failed; inspect stratum_workflow_status.tsv.gz.", call. = FALSE)
  artifacts <- lapply(artifact_files, readRDS)

  single_cell_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  layer1 <- .rc_merge_stratum_layer1(
    artifacts = artifacts,
    gem = gem,
    single_cell_genes = single_cell_genes,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  grn_meta_modules <- .rc_merge_stratum_meta_modules(artifacts)
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
  microcompass <- do.call(rc_run_microcompass, c(layer2_defaults, layer2_args))
  .rc_release_bpparam(layer2_param)
  rm(layer2_param)
  invisible(gc(verbose = FALSE))

  result <- list(
    schema_version = "regcompass_v2_global_metacell",
    model_mode = model_mode,
    strict_stratum_cols = group_cols,
    upstream_status = upstream_status,
    layer1 = layer1,
    grn_meta_modules = grn_meta_modules,
    microcompass = microcompass,
    params = list(
      shared_gem = TRUE,
      penalty_unit = "metacell",
      upstream_parallel_unit = "condition_sample_celltype",
      layer2_parallel_unit = "target_direction_by_metacell",
      parallel_backend = parallel_backend,
      upstream_workers = upstream_workers,
      layer2_workers = layer2_workers
    )
  )
  saveRDS(result, file.path(outdir, "regcompass_global_metacell_result.rds"))
  result
}
