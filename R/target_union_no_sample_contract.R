# Remove historical sample placeholders from the second-pass target-union path.

.rc_build_target_union_model_cache_sample_core <-
  .rc_build_target_union_model_cache

.rc_build_target_union_model_cache <- function(
    layer2, target_reactions,
    target_direction = c("both", "forward", "reverse")) {
  cache <- .rc_build_target_union_model_cache_sample_core(
    layer2 = layer2,
    target_reactions = target_reactions,
    target_direction = target_direction
  )
  for (key in names(cache)) cache[[key]]$sample_id <- NULL
  cache
}

.rc_score_existing_union_cache_sample_core <- .rc_score_existing_union_cache

.rc_score_existing_union_cache <- function(
    layer1, gem, model_cache, condition_col, celltype_col,
    omega = 0.95,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8,
    parallel = TRUE, BPPARAM = NULL, ...) {
  dots <- list(...)
  unknown <- setdiff(names(dots), "sample_col")
  if (length(unknown)) {
    stop(
      "Unsupported target-union scoring arguments: ",
      paste(unknown, collapse = ", "), call. = FALSE
    )
  }
  answer <- .rc_score_existing_union_cache_sample_core(
    layer1 = layer1,
    gem = gem,
    model_cache = model_cache,
    condition_col = condition_col,
    sample_col = NULL,
    celltype_col = celltype_col,
    omega = omega,
    solver = solver,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  if (is.data.frame(answer$lp_diagnostics) &&
      "sample_id" %in% colnames(answer$lp_diagnostics)) {
    answer$lp_diagnostics$sample_id <- NULL
  }
  answer$params$unit <- "metacell"
  answer$params$aggregation <- "none"
  answer$params$sample_column <- NULL
  answer
}

.rc_regcompass_step_target_union_sample_core <-
  rc_regcompass_step_target_union

#' Score direct database-linked reactions without sample aggregation
#'
#' @inheritParams rc_regcompass_step_target_union
#' @export
rc_regcompass_step_target_union <- function(
    layer1, meta_modules, layer2, gem, outdir,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  meta_modules$workflow_params$sample_col <- NULL
  layer1$workflow_params$sample_col <- NULL
  layer2$workflow_params$sample_col <- NULL
  answer <- .rc_regcompass_step_target_union_sample_core(
    layer1 = layer1,
    meta_modules = meta_modules,
    layer2 = layer2,
    gem = gem,
    outdir = outdir,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match,
    layer2_args = layer2_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
  answer$workflow_params$sample_col <- NULL
  answer$microcompass$workflow_params$sample_col <- NULL
  answer
}
