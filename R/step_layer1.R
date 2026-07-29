#' Project condition GRN effects and build RNA+ATAC reaction support
#'
#' Uses Pando's outer-heldout single-cell TF-by-ATAC projections and exact
#' SuperCell cell-to-metacell membership. Interactions are never reconstructed
#' from metacell means.
#'
#' @export
rc_regcompass_step_layer1 <- function(
    grn, metacells, meta_modules, gem, outdir,
    projection_component = "condition",
    comparison_support = c(
      "auto", "pairwise_common", "global_common"
    ),
    regulatory_alpha = 0.5,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer1", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  if (!identical(projection_component, "condition")) {
    stop("Only the condition OOF projection can enter Layer 1.",
         call. = FALSE)
  }
  comparison_support <- match.arg(comparison_support)
  gpr_and_method <- match.arg(gpr_and_method)
  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
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
  .rc_require_stage_gem(grn, gem, "grn")
  layer1 <- .rc_cell_first_projection_layer1(
    grn_result = grn$grn_result,
    metacell_object = metacells$metacell_object,
    membership = metacells$pooled$membership,
    gem = gem,
    metacell_meta = metacells$pooled$metacell_meta,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    projection_component = projection_component,
    comparison_support = comparison_support,
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
