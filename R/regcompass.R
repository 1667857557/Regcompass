#' Run the canonical GRN-first RegCompass workflow
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
    species = c("auto", "human", "mouse")) {
  model_mode <- match.arg(model_mode)
  parallel_backend <- match.arg(parallel_backend)
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
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(gem$model_info %||% list(), file.path(outdir, "00_model_info.rds"))
  saveRDS(medium_scenarios, file.path(outdir, "00_medium_scenarios.rds"))
  upstream <- .rc_phase_bpparam(upstream_workers, parallel_backend)
  on.exit(.rc_release_bpparam(upstream), add = TRUE)
  step1 <- rc_regcompass_step_grn(
    object = object, gem = gem,
    outdir = file.path(outdir, "01_single_cell_grn"),
    pfm = pfm, genome = genome,
    condition_col = condition_col, celltype_col = celltype_col,
    rna_assay = rna_assay, atac_assay = atac_assay,
    pando_args = pando_args, parallel = TRUE, BPPARAM = upstream
  )
  step2 <- rc_regcompass_step_metacells(
    object = object, outdir = file.path(outdir, "02_condition_metacells"),
    sample_col = sample_col, condition_col = condition_col,
    celltype_col = celltype_col, rna_assay = rna_assay,
    atac_assay = atac_assay, fragment_files = fragment_files,
    metacell_args = metacell_args
  )
  step3 <- rc_regcompass_step_meta_modules(
    grn = step1, metacells = step2, gem = gem,
    outdir = file.path(outdir, "03_meta_modules"),
    layer1_args = layer1_args
  )
  step4 <- rc_regcompass_step_layer1(
    metacells = step2, meta_modules = step3, gem = gem,
    outdir = file.path(outdir, "04_layer1"),
    regulatory_alpha = layer1_args$regulatory_alpha %||% 1,
    tau = layer1_args$tau %||% 0.20,
    gene_half_saturation = layer1_args$gene_half_saturation %||%
      getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE, BPPARAM = upstream
  )
  .rc_release_bpparam(upstream)
  upstream <- FALSE
  lp_param <- .rc_phase_bpparam(layer2_workers, parallel_backend)
  step5 <- tryCatch(
    rc_regcompass_step_layer2(
      layer1 = step4, meta_modules = step3, gem = gem,
      medium_scenarios = medium_scenarios,
      outdir = file.path(outdir, "05_layer2"),
      model_mode = model_mode, layer2_args = layer2_args,
      parallel = TRUE, BPPARAM = lp_param
    ),
    finally = .rc_release_bpparam(lp_param)
  )
  result <- rc_regcompass_step_results(
    grn = step1, metacells = step2, meta_modules = step3,
    layer1 = step4, layer2 = step5, gem = gem,
    outdir = file.path(outdir, "06_results"), species = species
  )
  result$params$execution_mode <- "one_shot"
  result$params$parallel_backend <- parallel_backend
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  result
}
