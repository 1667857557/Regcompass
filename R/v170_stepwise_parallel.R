# Stepwise v1.7 corrections are implemented under unique internal names and
# assigned explicitly to the public API, avoiding Collate-order redefinitions.

.rc_regcompass_step_metacells_base <- rc_regcompass_step_metacells
.rc_regcompass_step_metacells_v170 <- function(
    object, outdir,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list()) {
  # Emit the biological-replication warning before the legacy implementation
  # performs the same design check with its historical warning handler.
  .rc_condition_pool_design_summary(
    object@meta.data,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    strict_biological_defaults = FALSE
  )
  .rc_regcompass_step_metacells_base(
    object = object,
    outdir = outdir,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args
  )
}
rc_regcompass_step_metacells <- .rc_regcompass_step_metacells_v170

.rc_regcompass_step_meta_modules_base <- rc_regcompass_step_meta_modules
.rc_regcompass_step_meta_modules_v170 <- function(
    metacells, gem, outdir, pfm, genome,
    pando_args = list(),
    layer1_args = list(),
    parallel = TRUE,
    BPPARAM = NULL) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  answer <- .rc_regcompass_step_meta_modules_base(
    metacells = metacells,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    pando_args = pando_args,
    layer1_args = layer1_args,
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  answer$params$parallel <- parallel
  answer$params$bpparam_class <- if (is.null(BPPARAM)) {
    "auto_or_sequential"
  } else if (identical(BPPARAM, FALSE)) {
    "sequential"
  } else {
    class(BPPARAM)[[1L]]
  }
  saveRDS(answer, file.path(outdir, "step_meta_modules.rds"))
  answer
}
rc_regcompass_step_meta_modules <- .rc_regcompass_step_meta_modules_v170

.rc_regcompass_step_layer1_v170 <- function(
    metacells, meta_modules, gem, outdir,
    regulatory_alpha = 1,
    tau = 0.20,
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE,
    BPPARAM = NULL) {
  if (!inherits(metacells, "regcompass_metacell_step")) {
    stop(
      "`metacells` must be the output of `rc_regcompass_step_metacells()`.",
      call. = FALSE
    )
  }
  if (!inherits(meta_modules, "regcompass_meta_module_step")) {
    stop(
      "`meta_modules` must be the output of `rc_regcompass_step_meta_modules()`.",
      call. = FALSE
    )
  }
  if (!identical(metacells$params, meta_modules$workflow_params)) {
    stop(
      "Metacell and meta-module stages use different workflow metadata settings.",
      call. = FALSE
    )
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  params <- metacells$params
  layer1 <- .rc_build_condition_pooled_layer1(
    metacell_object = metacells$metacell_object,
    meta_modules = meta_modules$condition_modules,
    gem = gem,
    metacell_meta = metacells$pooled$metacell_meta,
    sample_col = params$sample_col,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    atac_assay = params$atac_assay,
    regulatory_alpha = regulatory_alpha,
    gpr_tau = tau,
    gene_half_saturation = gene_half_saturation,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  if (!identical(
    colnames(layer1$reaction_expression),
    as.character(layer1$unit_meta$pool_id)
  )) {
    stop(
      "Layer 1 reaction expression and unit metadata are not ordered identically.",
      call. = FALSE
    )
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(layer1, file.path(outdir, "step_layer1.rds"))
  layer1
}
rc_regcompass_step_layer1 <- .rc_regcompass_step_layer1_v170

.rc_regcompass_step_layer2_base <- rc_regcompass_step_layer2
.rc_regcompass_step_layer2_v170 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(),
    parallel = TRUE,
    BPPARAM = NULL) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  warning(
    paste(
      "Metacell-level scores are descriptive pseudo-observations and are not",
      "independent biological replicates."
    ),
    call. = FALSE
  )
  .rc_regcompass_step_layer2_base(
    layer1 = layer1,
    meta_modules = meta_modules,
    gem = gem,
    medium_scenarios = medium_scenarios,
    outdir = outdir,
    model_mode = model_mode,
    layer2_args = layer2_args,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
}
rc_regcompass_step_layer2 <- .rc_regcompass_step_layer2_v170
