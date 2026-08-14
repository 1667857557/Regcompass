#' Run the RegCompass workflow with cell-type-specific Pando routing
#'
#' Stage 1 resolves the Pando route independently for each retained broad cell
#' type. Cell types with at least two retained conditions use the common-
#' dictionary condition GRN; cell types with one retained condition use standard
#' Pando. RegCompass applies a separate target-model quality gate using the
#' selected-lambda final full-data Pando R-squared; OOF R-squared is diagnostic.
#'
#' Stage 2 builds one independent multimodal WNN graph per broad cell type with
#' all conditions jointly present during WNN construction and Walktrap
#' clustering. Condition splits parent memberships only after clustering. If a
#' minimum final metacell size is requested, repair uses that exact original
#' shared WNN within the same condition and broad cell type. The repaired
#' membership is the only membership used by downstream aggregation/projection.
#'
#' Layer 2 selects exactly one structural route. `model_mode =
#' "meta_module_gem"` uses CORDA2 by default. Set
#' `layer2_args$model_params$model_completion = "fastcore"` for the
#' supplementary FASTCORE route. `model_mode = "full_gem"` retains the complete
#' GEM and applies medium exchange bounds without structural reconstruction.
#'
#' @param object Paired-cell Seurat RNA+ATAC object.
#' @param gem Prepared genome-scale metabolic model.
#' @param outdir Persistent output directory.
#' @param genome Genome object matching ATAC coordinates and regulatory regions.
#' @param pfm Optional motif collection.
#' @param species `"auto"`, `"human"`, or `"mouse"`.
#' @param condition_col Condition metadata column, or `NULL` for standard Pando.
#' @param celltype_col Broad-cell-type metadata column.
#' @param cell_type Optional broad-cell-type subset.
#' @param rna_assay RNA assay name.
#' @param atac_assay ATAC assay name.
#' @param fragment_files Optional Stage 2 raw ATAC fragment input.
#' @param pando_args Stage 1 Pando configuration.
#' @param target_rsq_threshold RegCompass target-model quality threshold on the
#' selected-lambda final full-data Pando R-squared. Default `0.05`. This is
#' applied after Pando edge significance and does not redefine Pando
#' `active`/`significant`; OOF R-squared remains diagnostic only.
#' @param metacell_args Stage 2 WNN/metacell controls. When
#' `min_metacell_size > 1`, `min_merge_affinity` must be explicit.
#' @param meta_module_args Stage 3 controls.
#' @param layer1_args Layer 1 controls.
#' @param medium_scenarios Shared medium table. If omitted, species defaults are
#' used. Built-in concentration challenges do not automatically infer uptake
#' flux bounds from mM values.
#' @param model_mode Structural route: `"meta_module_gem"` or `"full_gem"`.
#' @param layer2_args Layer 2 controls.
#' @param workers Total RegCompass worker cap, default 10.
#' @param progress Show and persist stage progress.
#' @return A complete RegCompass result.
#' @export
rc_run_regcompass <- function(
    object, gem, outdir, genome,
    pfm = NULL,
    species = c("auto", "human", "mouse"),
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = NULL,
    pando_args = list(),
    target_rsq_threshold = 0.05,
    metacell_args = list(),
    meta_module_args = list(),
    layer1_args = list(),
    medium_scenarios = NULL,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(),
    workers = 10L,
    progress = getOption("RegCompassR.progress", TRUE)) {
  model_mode <- match.arg(model_mode)
  target_rsq_threshold <- .rc_target_rsq_threshold(target_rsq_threshold)
  bundles <- list(
    pando_args = pando_args,
    metacell_args = metacell_args,
    meta_module_args = meta_module_args,
    layer1_args = layer1_args,
    layer2_args = layer2_args
  )
  invalid <- names(bundles)[!vapply(bundles, is.list, logical(1))]
  if (length(invalid)) {
    stop("Workflow argument bundles must be lists: ",
         paste(invalid, collapse = ", "), call. = FALSE)
  }
  allowed_layer1 <- c("gpr_and_method", "gene_half_saturation")
  unknown_layer1 <- setdiff(names(layer1_args), allowed_layer1)
  if (length(unknown_layer1)) {
    stop("Unknown `layer1_args`: ", paste(unknown_layer1, collapse = ", "),
         call. = FALSE)
  }
  worker_config <- .rc_stage_worker_config(workers, argument = "workers")
  worker_limit <- worker_config$worker_limit
  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  if (is.null(medium_scenarios)) {
    default_scenario <- if (identical(species, "human")) {
      "normal_human_plasma"
    } else {
      "mouse_plasma"
    }
    medium_scenarios <- rc_make_medium_scenarios(
      gem,
      scenario = default_scenario,
      species = species
    )
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  step1 <- rc_regcompass_step_grn(
    object = object, gem = gem,
    outdir = file.path(outdir, "01_single_cell_grn"),
    genome = genome, pfm = pfm, species = species,
    condition_col = condition_col, celltype_col = celltype_col,
    cell_type = cell_type, rna_assay = rna_assay,
    atac_assay = atac_assay,
    pando_args = pando_args,
    target_rsq_threshold = target_rsq_threshold,
    workers = worker_limit,
    progress = progress
  )
  step2 <- rc_regcompass_step_metacells(
    object = object,
    outdir = file.path(outdir, "02_metacells"),
    condition_col = condition_col, celltype_col = celltype_col,
    cell_type = NULL, rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args,
    workers = worker_limit,
    progress = progress,
    grn = step1
  )
  step3 <- rc_regcompass_step_meta_modules(
    grn = step1, metacells = step2, gem = gem,
    outdir = file.path(outdir, "03_meta_modules"),
    meta_module_args = meta_module_args, progress = progress
  )
  step4 <- rc_regcompass_step_layer1(
    grn = step1, metacells = step2, meta_modules = step3, gem = gem,
    outdir = file.path(outdir, "04_layer1"),
    gpr_and_method = layer1_args$gpr_and_method %||% "min",
    gene_half_saturation = layer1_args$gene_half_saturation %||%
      getOption("RegCompassR.cpm_half_saturation", 1),
    workers = worker_limit,
    progress = progress
  )
  step5 <- rc_regcompass_step_layer2(
    layer1 = step4, meta_modules = step3, gem = gem,
    medium_scenarios = medium_scenarios,
    outdir = file.path(outdir, "05_layer2"),
    model_mode = model_mode, layer2_args = layer2_args,
    workers = worker_limit,
    progress = progress
  )
  result <- rc_regcompass_step_results(
    grn = step1, metacells = step2, meta_modules = step3,
    layer1 = step4, layer2 = step5, gem = gem,
    outdir = file.path(outdir, "06_results"),
    species = species, progress = progress
  )
  result$params$execution_mode <- "one_shot"
  result$params$analysis_mode <- step1$params$analysis_mode
  result$params$condition_coefficients_calculated <-
    isTRUE(step1$grn_result$condition_coefficients_calculated)
  result$params$cell_type_analysis_mode <-
    step1$grn_result$cell_type_analysis_mode
  result$params$requested_condition_col <- condition_col
  result$params$effective_condition_col <- step1$params$condition_col
  result$params$target_rsq_threshold <- target_rsq_threshold
  result$params$target_rsq_metric <- "selected_lambda_final_full_data_rsq"
  result$params$oof_rsq_role <- "diagnostic_only"
  design <- step2$pooled$input_design
  result$params$native_supercell_api <- design$native_supercell_api
  result$params$native_supercell_inputs <- c(
    graph_group = design$graph_group_argument,
    condition = design$condition_argument
  )
  result$params$metacell_graph_method <- design$graph_method
  result$params$metacell_graph_scope <- design$graph_scope
  result$params$metacell_condition_scope <- design$condition_scope
  result$params$metacell_membership_split_timing <-
    design$membership_split_timing
  result$params$metacell_repair_timing <- design$repair_timing
  result$params$metacell_repair_geometry <- design$repair_geometry
  result$params$metacell_min_merge_affinity <- design$min_merge_affinity
  result$params$metacell_unresolved_small_policy <-
    design$unresolved_small_policy
  result$params$metacell_modality_weighting <- design$modality_weighting
  result$params$metacell_atac_aggregation <- design$atac_aggregation_method
  result$params$fragment_files_supplied <- isTRUE(design$fragment_files_supplied)
  result$params$temporary_combined_stratum <- FALSE
  result$params$workers <- worker_limit
  result$params$requested_workers <- worker_config$requested_workers
  result$params$detected_cpu_capacity <- worker_config$detected_cpu_capacity
  result$params$reserved_cpus <- worker_config$reserved_cpus
  result$params$parallel_backend <- worker_config$actual_backend
  result$params$parallel_worker_policy <-
    "min(independent_tasks,requested_workers,max(1,detected_cpus-2))"
  result
}
