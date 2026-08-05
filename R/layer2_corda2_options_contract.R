# Clean parameter contract for corrected Python CORDA2.

.rc_layer2_corda_options_algorithm_base <- .rc_layer2_corda_options

.rc_layer2_corda_options <- function(model_params = list()) {
  if (!is.list(model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  requested <- as.character(model_params$model_completion %||% "fastcore")
  if (length(requested) != 1L || is.na(requested)) {
    stop("`model_completion` must be `fastcore` or `corda2`.",
         call. = FALSE)
  }
  if (identical(requested, "fastcore")) {
    return(list(
      model_completion = "fastcore",
      requested_model_completion = "fastcore",
      algorithm = "add_only_compact_FASTCORE"
    ))
  }
  if (!requested %in% c("corda2", "corda", "corda_like")) {
    stop("`model_completion` must be `fastcore` or `corda2`.",
         call. = FALSE)
  }
  if ("corda_seed" %in% names(model_params)) {
    stop(
      "`corda_seed` is not a CORDA2 parameter. CORDA2 uses deterministic ",
      "multiplicative cost updates and no random-noise process.",
      call. = FALSE
    )
  }
  obsolete <- intersect(
    names(model_params),
    c("corda_other_penalty", "corda_negative_penalty")
  )
  if (length(obsolete)) {
    stop(
      "Obsolete weighted-FASTCORE parameter(s): ",
      paste(obsolete, collapse = ", "), ".",
      call. = FALSE
    )
  }
  scalar <- .rc_corda_scalar
  integer <- .rc_corda_integer
  flag <- .rc_corda_flag
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
  answer <- list(
    model_completion = "corda2",
    requested_model_completion = "corda2",
    penalty_factor = scalar(
      model_params$corda2_penalty_factor %||%
        model_params$corda_penalty_factor %||%
        model_params$corda_gamma %||% 100,
      "corda2_penalty_factor", 1, Inf
    ),
    cost_increase = scalar(
      model_params$corda2_cost_increase %||%
        model_params$corda_cost_increase %||%
        model_params$corda_kappa %||% 1.01,
      "corda2_cost_increase", 1, Inf
    ),
    target_flux = scalar(
      model_params$corda2_target_flux %||%
        model_params$corda_tflux %||%
        model_params$corda_epsilon %||% 1,
      "corda2_target_flux", .Machine$double.eps, Inf
    ),
    redundancies = integer(
      model_params$corda2_redundancies %||%
        model_params$corda_n %||% 3L,
      "corda2_redundancies", 1L
    ),
    support = integer(
      model_params$corda2_support %||%
        model_params$corda_support %||%
        model_params$corda_p %||% 5L,
      "corda2_support", 1L
    ),
    flux_tolerance = scalar(
      model_params$corda2_flux_tolerance %||%
        model_params$corda_flux_tolerance %||% 1e-8,
      "corda2_flux_tolerance", .Machine$double.eps, Inf
    ),
    upper_bound = 1e6,
    medium_confidence_threshold = scalar(
      model_params$corda_medium_confidence_threshold %||% 0.75,
      "corda_medium_confidence_threshold", 0, 1
    ),
    negative_confidence_threshold = scalar(
      model_params$corda_negative_confidence_threshold %||% 0.10,
      "corda_negative_confidence_threshold", 0, 1
    ),
    regulatory_weight = scalar(
      model_params$corda_regulatory_weight %||% 0.20,
      "corda_regulatory_weight", 0, 1
    ),
    include_evidence_outside_modules = flag(
      model_params$corda_include_evidence_outside_modules %||% TRUE,
      "corda_include_evidence_outside_modules"
    ),
    max_medium_confidence_reactions = max_mc,
    evidence_definition = paste(
      "(1 - regulatory_weight) * max(within-cell-type RNA percentile,",
      "multiome percentile) + regulatory_weight * regulatory support"
    ),
    algorithm =
      "resendislab_python_CORDA2_corrected_redundant_path_assessment",
    python_reference_repository = "resendislab/corda",
    python_reference_commit =
      "c02e06d50606bf93f23d8f2e6d6ade0e996ca70e",
    python_defaults = c(
      redundancies = 3,
      penalty_factor = 100,
      support = 5,
      cost_increase = 1.01,
      target_flux = 1,
      upper_bound = 1e6
    ),
    intentional_corrections = c(
      "maximize remaining medium-confidence directional flux",
      "block the opposite direction of a reversible target",
      "apply penalties from each directional variable's own confidence"
    ),
    deterministic = TRUE,
    random_noise = FALSE
  )
  if (answer$cost_increase <= 1) {
    stop("`corda2_cost_increase` must be greater than 1.", call. = FALSE)
  }
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
