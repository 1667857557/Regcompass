#' Construct complete-GPR cores and biological meta-modules
#'
#' Stage 3 defines biological reaction membership only. It does not run
#' FASTCORE and does not construct a GEM. Reactions from all modules are
#' deduplicated into a merged meta-module catalogue for downstream Layer 1 and
#' medium-specific union-GEM construction.
#'
#' @export
rc_regcompass_step_meta_modules <- function(
    grn, metacells, gem, outdir,
    layer1_args = list(),
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("meta_modules", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  if (!identical(.rc_workflow_signature(grn),
                 .rc_workflow_signature(metacells))) {
    stop(
      "GRN and metacell stages use different metadata or assay settings.",
      call. = FALSE
    )
  }
  if (!is.list(layer1_args)) {
    stop("`layer1_args` must be a list.", call. = FALSE)
  }
  obsolete <- intersect(
    names(layer1_args), c("local_fastcore", "local_fastcore_args")
  )
  if (length(obsolete)) {
    stop(
      paste0(
        "Local FASTCORE was removed. Delete obsolete `layer1_args` fields: ",
        paste(obsolete, collapse = ", "),
        ". Configure global FASTCORE through `layer2_args$model_params`."
      ),
      call. = FALSE
    )
  }
  .rc_require_stage_gem(grn, gem, "grn")
  validated_gem <- rc_validate_gem(gem)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  group_coverage <- .rc_validate_grn_metacell_group_coverage(
    grn_result = grn$grn_result,
    metacell_meta = metacells$pooled$metacell_meta,
    condition_col = metacells$params$condition_col,
    celltype_col = metacells$params$celltype_col
  )
  .rc_write_tsv_gz(
    group_coverage,
    file.path(outdir, "grn_metacell_group_coverage.tsv.gz")
  )
  condition_modules <- .rc_build_condition_meta_modules(
    grn$grn_result, gem, outdir, layer1_args
  )
  condition_modules$grn_metacell_group_coverage <- group_coverage
  if (!is.data.frame(condition_modules$reaction_membership) ||
      !nrow(condition_modules$reaction_membership)) {
    stop(
      "Meta-module construction produced no reaction membership.",
      call. = FALSE
    )
  }
  missing <- setdiff(
    unique(as.character(condition_modules$reaction_membership$reaction_id)),
    colnames(validated_gem$S)
  )
  if (length(missing)) {
    stop(
      "Meta-module reactions absent from the GEM: ",
      paste(utils::head(missing, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  merged_modules <- .rc_merge_stratum_meta_modules(list(list(
    group_id = "condition_pooled",
    grn_meta_modules = condition_modules
  )))
  if (!is.data.frame(merged_modules$merged_core_reactions) ||
      !nrow(merged_modules$merged_core_reactions)) {
    stop(
      "No complete-GPR merged core reactions remain after module merging.",
      call. = FALSE
    )
  }
  answer <- list(
    condition_modules = condition_modules,
    merged_modules = merged_modules,
    group_coverage = group_coverage,
    workflow_params = metacells$params,
    grn_params = grn$params,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      layer1_args = layer1_args,
      feasibility_completion = "layer2_medium_specific_only",
      merge_creates_gem = FALSE
    )
  )
  class(answer) <- c("regcompass_meta_module_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  saveRDS(
    condition_modules,
    file.path(outdir, "condition_meta_modules.rds")
  )
  saveRDS(
    merged_modules,
    file.path(outdir, "merged_meta_modules.rds")
  )
  saveRDS(answer, file.path(outdir, "step_meta_modules.rds"))
  answer
}
