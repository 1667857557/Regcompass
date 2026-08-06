# Public parameter contract for the original MATLAB CORDA2 algorithm.

.rc_corda2_original_args <- function(model_params) {
  args <- model_params$corda2_args %||% list()
  if (!is.list(args)) {
    stop("`corda2_args` must be a named list.", call. = FALSE)
  }
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))))) {
    stop("`corda2_args` must be a named list.", call. = FALSE)
  }
  allowed <- c("MCxNCthresh", "constraint", "constrainby", "om", "ci")
  unknown <- setdiff(names(args), allowed)
  if (length(unknown)) {
    stop(
      "Unknown CORDA2 parameter(s) in `corda2_args`: ",
      paste(unknown, collapse = ", "),
      ". Allowed names are MCxNCthresh, constraint, constrainby, om and ci.",
      call. = FALSE
    )
  }
  MCxNCthresh <- .rc_corda_scalar(
    args$MCxNCthresh %||% 2,
    "corda2_args$MCxNCthresh", 0, Inf
  )
  constraint <- .rc_corda_scalar(
    args$constraint %||% 1,
    "corda2_args$constraint", 0, Inf
  )
  constrainby <- as.character(args$constrainby %||% "val")
  if (length(constrainby) != 1L || is.na(constrainby) ||
      !constrainby %in% c("val", "perc")) {
    stop("`corda2_args$constrainby` must be `val` or `perc`.",
         call. = FALSE)
  }
  om <- .rc_corda_scalar(args$om %||% 1e4, "corda2_args$om", 0, Inf)
  if (om <= 0) stop("`corda2_args$om` must be positive.", call. = FALSE)
  ci <- .rc_corda_scalar(args$ci %||% 0.01, "corda2_args$ci", 0, Inf)
  list(
    MCxNCthresh = MCxNCthresh,
    constraint = constraint,
    constrainby = constrainby,
    om = om,
    ci = ci
  )
}

.rc_layer2_corda_options <- function(model_params = list()) {
  if (!is.list(model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  requested <- as.character(model_params$model_completion %||% "corda2")
  if (length(requested) != 1L || is.na(requested)) {
    stop("`model_completion` must be `corda2` or `fastcore`.",
         call. = FALSE)
  }
  if (identical(requested, "fastcore")) {
    return(list(
      model_completion = "fastcore",
      requested_model_completion = "fastcore",
      algorithm = "add_only_compact_FASTCORE"
    ))
  }
  if (!requested %in% c("corda2", "corda")) {
    stop("`model_completion` must be `corda2` or `fastcore`.",
         call. = FALSE)
  }

  obsolete <- intersect(names(model_params), c(
    "corda2_redundancies", "corda2_penalty_factor", "corda2_support",
    "corda_penalty_factor", "corda_support", "corda_n", "corda_p",
    "corda2_target_flux", "corda_tflux", "corda_epsilon",
    "corda2_cost_increase", "corda_cost_increase", "corda_kappa",
    "corda2_flux_tolerance", "corda_flux_tolerance", "corda_seed",
    "corda_gamma", "corda_other_penalty", "corda_negative_penalty"
  ))
  if (length(obsolete)) {
    stop(
      "Unsupported Python-CORDA parameter(s): ",
      paste(obsolete, collapse = ", "),
      ". Use `corda2_args` with the original MATLAB CORDA2 names.",
      call. = FALSE
    )
  }

  original <- .rc_corda2_original_args(model_params)
  max_mc <- suppressWarnings(as.numeric(
    model_params$corda_max_medium_confidence_reactions %||% Inf
  ))
  if (length(max_mc) != 1L || is.na(max_mc) || max_mc < 0) {
    stop(
      "`corda_max_medium_confidence_reactions` must be non-negative or Inf.",
      call. = FALSE
    )
  }
  if (is.finite(max_mc)) max_mc <- as.integer(floor(max_mc))

  answer <- c(list(
    model_completion = "corda2",
    requested_model_completion = requested
  ), original, list(
    flux_threshold = 1e-7,
    baseline_cost = 1e-3,
    output_bound = 1000,
    medium_confidence_threshold = .rc_corda_scalar(
      model_params$corda_medium_confidence_threshold %||% 0.75,
      "corda_medium_confidence_threshold", 0, 1
    ),
    negative_confidence_threshold = .rc_corda_scalar(
      model_params$corda_negative_confidence_threshold %||% 0.10,
      "corda_negative_confidence_threshold", 0, 1
    ),
    regulatory_weight = .rc_corda_scalar(
      model_params$corda_regulatory_weight %||% 0.20,
      "corda_regulatory_weight", 0, 1
    ),
    include_evidence_outside_modules = .rc_corda_flag(
      model_params$corda_include_evidence_outside_modules %||% TRUE,
      "corda_include_evidence_outside_modules"
    ),
    max_medium_confidence_reactions = max_mc,
    evidence_definition = paste(
      "(1 - regulatory_weight) * max(within-cell-type RNA percentile,",
      "multiome percentile) + regulatory_weight * regulatory support"
    ),
    algorithm = "schultzdre_MATLAB_CORDA2_original_semantics",
    reference_repository = "schultzdre/Constraint-Based-Modeling",
    reference_file = "CORDA2.m",
    random_noise = FALSE,
    source_semantics = "original_matlab_corda2"
  ))
  if (answer$negative_confidence_threshold >
      answer$medium_confidence_threshold) {
    stop(
      "`corda_negative_confidence_threshold` must not exceed ",
      "`corda_medium_confidence_threshold`.",
      call. = FALSE
    )
  }
  answer
}
