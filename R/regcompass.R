#' Run the RegCompass workflow with cell-type-specific Pando routing
#'
#' Stage 1 resolves the Pando route independently for each retained broad cell
#' type. Cell types with at least two retained conditions use the common-
#' dictionary condition GRN; cell types with one retained condition use standard
#' Pando and do not receive condition coefficients. Stage 2 builds one
#' independent multimodal WNN graph per broad cell type and keeps final
#' metacells condition-pure. Optional raw ATAC fragment files are aggregated to
#' the final metacell membership and recounted before metacell TF-IDF.
#'
#' Layer 2 selects exactly one structural route. `model_mode =
#' "meta_module_gem"` uses original MATLAB CORDA2 by default. Set
#' `layer2_args$model_params$model_completion = "fastcore"` for the
#' supplementary add-only FASTCORE route. `model_mode = "full_gem"` retains the
#' complete GEM, applies the shared medium as exchange-reaction bounds only, and
#' uses directional COMPASS vmax to identify infeasible targets; it does not
#' execute FASTCORE or CORDA2. All routes use the COMPASS reaction-cost scale,
#' with missing expression and structural reaction roles assigned cost 1.
#'
#' When `medium_scenarios` is omitted, Human-GEM uses
#' `"normal_human_plasma"` and Mouse-GEM uses `"mouse_plasma"`.
#'
#' One `workers` argument controls the maximum process count for all parallel
#' stages. The default is 10. RegCompass automatically selects SOCK/Snow workers
#' on Windows and Multicore workers on Linux/macOS. Two detected logical CPUs are
#' reserved globally, and each stage/sub-step uses only the number of independent
#' tasks it can actually execute concurrently.
#'
#' @param fragment_files Optional Stage 2 fragment input. Use `NULL`/`FALSE` to
#' aggregate the existing ATAC matrix, a fragment path (or named path vector for
#' multiple samples), or a data frame with `fragment_file`, `object_cell`, and
#' `fragment_barcode`. Fragment routing does not make sample a metacell stratum.
#' @param workers Total RegCompass worker cap, default 10. Users may increase or
#' decrease it. The effective limit is
#' `min(workers, max(1, detected logical CPUs - 2))`; the backend is selected
#' automatically for the operating system.
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
    metacell_args = list(),
    meta_module_args = list(),
    layer1_args = list(),
    medium_scenarios = NULL,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(),
    workers = 10L,
    progress = getOption("RegCompassR.progress", TRUE)) {
  model_mode <- match.arg(model_mode)
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
