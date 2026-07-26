.rc_compact_microcompass <- function(layer2) {
  if (!is.list(layer2)) {
    stop("`layer2` must be a completed microCOMPASS result.", call. = FALSE)
  }
  required <- c(
    "score", "penalty", "vmax", "feasible", "target_direction",
    "unit_meta", "params", "model_mode"
  )
  missing <- setdiff(required, names(layer2))
  if (length(missing)) {
    stop(
      "Layer 2 lacks fields required for downstream condition analysis: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  fields <- c(
    "score", "relative_penalty_rank", "score_semantics",
    "noninformative_target", "primary_output", "primary_output_semantics",
    "penalty", "vmax", "feasible", "evaluated", "target_direction",
    "direction_diagnostics", "medium_scenarios", "model_mode",
    "model_cache_summary", "model_file_manifest", "unit_meta", "params",
    "method", "evidence_policy", "union_gem_policy"
  )
  out <- layer2[intersect(fields, names(layer2))]
  out$params$sample_column <- NULL
  out$params$unit <- "metacell"
  out$params$aggregation <- "none"
  attr(out, "omitted_stage5_fields") <- setdiff(names(layer2), names(out))
  class(out) <- unique(c(
    "regcompass_compact_microcompass", class(layer2), "list"
  ))
  out
}
