# Canonical runner contract for the medium-specific union-GEM architecture.

.rc_reject_local_fastcore_api <- function(layer1_args) {
  if (!is.list(layer1_args)) {
    stop("`layer1_args` must be a list.", call. = FALSE)
  }
  obsolete <- intersect(
    names(layer1_args),
    c("local_fastcore", "local_fastcore_args")
  )
  if (length(obsolete)) {
    stop(
      "Local FASTCORE was removed. Delete `",
      paste(obsolete, collapse = "` and `"),
      "` from `layer1_args`. Configure the single medium-specific global ",
      "FASTCORE through `layer2_args$model_params` ",
      "(`completion_time_limit`, `fastcore_epsilon`, ",
      "`max_support_reactions`, and `strict`).",
      call. = FALSE
    )
  }
  layer1_args
}

.rc_run_regcompass_before_union_contract <- rc_run_regcompass

#' Run the canonical GRN-first RegCompass workflow
#'
#' Stage 3 constructs biological meta-modules and a merged reaction catalogue.
#' No local FASTCORE is run. In `model_mode = "meta_module_gem"`, Stage 5
#' constructs one union GEM per medium scenario and performs the only FASTCORE
#' completion on the merged biological reaction set.
#'
#' @param upstream_workers Worker count for GRN inference and Layer 1
#'   reaction-expression calculation. Defaults to 6. Set to 1 for serial
#'   upstream execution.
#' @param layer2_workers Worker count for Layer 2 LP scoring. Defaults to 30.
#'   Set to 1 for serial Layer 2 execution.
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
    upstream_workers = 6L,
    layer2_workers = 30L,
    species = c("auto", "human", "mouse"),
    progress = getOption("RegCompassR.progress", TRUE)) {
  layer1_args <- .rc_reject_local_fastcore_api(layer1_args)
  .rc_run_regcompass_before_union_contract(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    model_mode = model_mode,
    medium_scenarios = medium_scenarios,
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args,
    layer2_args = layer2_args,
    upstream_workers = upstream_workers,
    layer2_workers = layer2_workers,
    species = species,
    progress = progress
  )
}
