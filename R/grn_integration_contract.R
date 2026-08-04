# Final cross-package argument contract for Stage 1 Pando routing.

.rc_route_pando_infer_args_integration_base <-
  .rc_route_pando_infer_args

.rc_require_condition_pando_version <- function(condition_types) {
  if (!length(condition_types)) return(invisible(TRUE))
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Package 'Pando' is required for condition-GRN Stage 1.",
         call. = FALSE)
  }
  installed <- utils::packageVersion("Pando")
  required <- utils::package_version("2.0.3")
  if (installed < required) {
    stop(
      "Condition-GRN projection requires Pando >= 2.0.3 for ",
      "case-safe target matching and strict fitted-cell alignment; installed ",
      as.character(installed), ".", call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_normalize_condition_pando_layer_args <- function(
    args, condition_types = character()) {
  layer_fields <- c("rna_layer", "peak_layer", "peak_value_type")
  supplied <- intersect(names(args), layer_fields)
  if (!length(condition_types) || !length(supplied)) {
    return(list(args = args, supplied = character()))
  }
  expected <- list(
    rna_layer = "data",
    peak_layer = "data",
    peak_value_type = "normalized"
  )
  invalid <- vapply(supplied, function(field) {
    value <- args[[field]]
    !is.character(value) || length(value) != 1L || is.na(value) ||
      !identical(as.character(value), expected[[field]])
  }, logical(1))
  if (any(invalid)) {
    stop(
      "RegCompass Stage 1 currently requires Pando condition projection ",
      "layers rna_layer='data', peak_layer='data', and ",
      "peak_value_type='normalized'. Unsupported override(s): ",
      paste(supplied[invalid], collapse = ", "), ".",
      call. = FALSE
    )
  }
  args[supplied] <- NULL
  list(args = args, supplied = supplied)
}

.rc_route_pando_infer_args <- function(
    args, condition_types = character(), standard_types = character()) {
  .rc_require_condition_pando_version(condition_types)
  normalized <- .rc_normalize_condition_pando_layer_args(
    args, condition_types
  )
  answer <- .rc_route_pando_infer_args_integration_base(
    args = normalized$args,
    condition_types = condition_types,
    standard_types = standard_types
  )
  if (length(normalized$supplied)) {
    confirmed <- data.frame(
      route = "condition_grn",
      argument = normalized$supplied,
      action = "confirmed_canonical_default",
      stringsAsFactors = FALSE
    )
    answer$diagnostics <- .rc_bind_frames_fill(list(
      answer$diagnostics, confirmed
    ))
  }
  answer
}
