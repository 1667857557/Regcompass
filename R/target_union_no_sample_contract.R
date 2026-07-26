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
