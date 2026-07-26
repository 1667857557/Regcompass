# Public Stage 3 contract adapter for the GRN evidence mode.

.rc_regcompass_step_meta_modules_core <- rc_regcompass_step_meta_modules

rc_regcompass_step_meta_modules <- function(
    grn, metacells, gem, outdir,
    meta_module_args = list(),
    progress = getOption("RegCompassR.progress", TRUE)) {
  answer <- .rc_regcompass_step_meta_modules_core(
    grn = grn,
    metacells = metacells,
    gem = gem,
    outdir = outdir,
    meta_module_args = meta_module_args,
    progress = progress
  )
  multitask <- identical(
    grn$grn_result$grn_mode %||% grn$params$grn_mode %||%
      "legacy_condition_pando",
    "multitask_shared_backbone"
  )
  if (multitask) {
    answer$params$core_definition <-
      "condition_celltype_stability_selected_targets_complete_gpr"
    answer$params$regulatory_evidence <- paste(
      "shared cell-type candidate universe with condition-effective",
      "global-plus-deviation coefficients and stability selection"
    )
    answer$params$condition_core_policy <- paste(
      "a reaction is core in one condition-celltype group only when at least",
      "one complete GPR AND branch is contained in that group's active target",
      "gene set"
    )
    answer$params$merge_creates_gem <- FALSE
    saveRDS(answer, file.path(outdir, "step_meta_modules.rds"))
  }
  answer
}
