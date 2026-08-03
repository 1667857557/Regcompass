#' Run the RegCompass workflow with automatic Pando mode selection
#'
#' Stage 1 selects condition-aware Pando only when at least two condition levels
#' are present. Otherwise it uses original Pando `infer_grn()` and calculates no
#' condition coefficients. Stage 2 builds one independent multimodal WNN graph
#' per broad cell type while pooling all conditions within that graph; adaptive
#' RNA/ATAC modality weights and neighbours are learned jointly, and condition
#' splits parent membership only after graph clustering.
#'
#' When `medium_scenarios` is omitted, Human-GEM uses
#' `"normal_human_plasma"` and Mouse-GEM uses `"mouse_plasma"`.
#'
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
    pando_args = list(),
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
  if ("tau" %in% names(layer1_args)) {
    stop(
      "The retired `tau`/Boltzmann Layer 1 transform is not supported.",
      call. = FALSE
    )
  }
  allowed_layer1 <- c(
    "projection_component", "comparison_support", "regulatory_alpha",
    "gpr_and_method", "gene_half_saturation"
  )
  unknown_layer1 <- setdiff(names(layer1_args), allowed_layer1)
  if (length(unknown_layer1)) {
    stop("Unknown `layer1_args`: ", paste(unknown_layer1, collapse = ", "),
         call. = FALSE)
  }
  if (!is.null(layer1_args$regulatory_alpha) &&
      !isTRUE(all.equal(as.numeric(layer1_args$regulatory_alpha), 1))) {
    stop("Canonical RegCompass requires `regulatory_alpha = 1`.",
         call. = FALSE)
  }
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
  saveRDS(gem$model_info %||% list(), file.path(outdir, "00_model_info.rds"))
  saveRDS(medium_scenarios, file.path(outdir, "00_medium_scenarios.rds"))

  upstream_config <- .rc_stage_worker_config(
    upstream_workers, argument = "upstream_workers"
  )
  layer2_config <- .rc_stage_worker_config(
    layer2_workers, argument = "layer2_workers"
  )
  step1 <- .rc_with_stage_workers(
    upstream_config$workers,
    argument = "upstream_workers",
    FUN = function(param, config) {
      rc_regcompass_step_grn(
        object = object, gem = gem,
        outdir = file.path(outdir, "01_single_cell_grn"),
        genome = genome, pfm = pfm, species = species,
        condition_col = condition_col, celltype_col = celltype_col,
        cell_type = cell_type, rna_assay = rna_assay,
        atac_assay = atac_assay, fragment_files = fragment_files,
        pando_args = pando_args,
        parallel = !identical(config$actual_backend, "serial"),
        BPPARAM = param, progress = progress
      )
    }
  )
  step2 <- rc_regcompass_step_metacells(
    object = object,
    outdir = file.path(outdir, "02_metacells"),
    condition_col = condition_col, celltype_col = celltype_col,
    cell_type = NULL, rna_assay = rna_assay,
    atac_assay = atac_assay, fragment_files = fragment_files,
    metacell_args = metacell_args, progress = progress,
    grn = step1
  )
  step3 <- rc_regcompass_step_meta_modules(
    grn = step1, metacells = step2, gem = gem,
    outdir = file.path(outdir, "03_meta_modules"),
    meta_module_args = meta_module_args, progress = progress
  )
  step4 <- .rc_with_stage_workers(
    upstream_config$workers,
    argument = "upstream_workers",
    FUN = function(param, config) {
      rc_regcompass_step_layer1(
        grn = step1, metacells = step2, meta_modules = step3, gem = gem,
        outdir = file.path(outdir, "04_layer1"),
        projection_component = layer1_args$projection_component %||% "condition",
        comparison_support = layer1_args$comparison_support %||% "auto",
        regulatory_alpha = 1,
        gpr_and_method = layer1_args$gpr_and_method %||% "min",
        gene_half_saturation = layer1_args$gene_half_saturation %||%
          getOption("RegCompassR.cpm_half_saturation", 1),
        parallel = !identical(config$actual_backend, "serial"),
        BPPARAM = param, progress = progress
      )
    }
  )
  step5 <- .rc_with_stage_workers(
    layer2_config$workers,
    argument = "layer2_workers",
    FUN = function(param, config) {
      rc_regcompass_step_layer2(
        layer1 = step4, meta_modules = step3, gem = gem,
        medium_scenarios = medium_scenarios,
        outdir = file.path(outdir, "05_layer2"),
        model_mode = model_mode, layer2_args = layer2_args,
        parallel = !identical(config$actual_backend, "serial"),
        BPPARAM = param, progress = progress
      )
    }
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
    identical(step1$params$analysis_mode, "condition_grn")
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
  result$params$temporary_combined_stratum <- FALSE
  result$params$upstream_workers <- upstream_config$workers
  result$params$layer2_workers <- layer2_config$workers
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  result
}
