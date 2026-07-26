#' Build integrated RNA+ATAC reaction expression
#'
#' Uses condition and cell-type metadata from the validated metacell stage. The
#' canonical Layer 1 contract contains no biological-sample field.
#'
#' @export
rc_regcompass_step_layer1 <- function(
    metacells, meta_modules, gem, outdir,
    regulatory_alpha = 1,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer1", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  gpr_and_method <- match.arg(gpr_and_method)
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!identical(metacells$params, meta_modules$workflow_params)) {
    stop(
      "Metacell and meta-module stages use different workflow settings.",
      call. = FALSE
    )
  }
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  params <- metacells$params
  layer1 <- .rc_build_condition_pooled_layer1(
    metacell_object = metacells$metacell_object,
    meta_modules = meta_modules$condition_modules,
    gem = gem,
    metacell_meta = metacells$pooled$metacell_meta,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    atac_assay = params$atac_assay,
    regulatory_alpha = regulatory_alpha,
    gpr_and_method = gpr_and_method,
    gene_half_saturation = gene_half_saturation,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  layer1$workflow_params <- params
  layer1$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  class(layer1) <- c("regcompass_layer1_step", "list")
  .rc_validate_layer1_stage(
    layer1,
    workflow_params = params,
    gem = gem,
    argument = "layer1"
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer1 <- .rc_step_monitor_finish(layer1, monitor)
  saveRDS(layer1, file.path(outdir, "step_layer1.rds"))
  layer1
}
