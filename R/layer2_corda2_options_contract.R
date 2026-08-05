# Public parameter contract for exact Python CORDA2 execution semantics.

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

  non_source <- intersect(names(model_params), c(
    "corda2_cost_increase", "corda_cost_increase", "corda_kappa",
    "corda2_target_flux", "corda_tflux", "corda_epsilon",
    "corda2_flux_tolerance", "corda_flux_tolerance", "corda_seed",
    "corda_gamma", "corda_n", "corda_p",
    "corda_other_penalty", "corda_negative_penalty"
  ))
  if (length(non_source)) {
    stop(
      "Exact Python CORDA2 does not expose parameter(s): ",
      paste(non_source, collapse = ", "),
      ". The source fixes CI=1.01, tflux=1 and UPPER=1e6; only n, ",
      "penalty_factor and support are constructor controls.",
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
    redundancies = integer(
      model_params$corda2_redundancies %||% 3L,
      "corda2_redundancies", 0L
    ),
    penalty_factor = scalar(
      model_params$corda2_penalty_factor %||%
        model_params$corda_penalty_factor %||% 100,
      "corda2_penalty_factor", -Inf, Inf
    ),
    support = integer(
      model_params$corda2_support %||%
        model_params$corda_support %||% 5L,
      "corda2_support", 0L
    ),
    cost_increase = 1.01,
    target_flux = 1,
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
    algorithm = "resendislab_python_CORDA2_c02e06d_exact_semantics",
    python_reference_repository = "resendislab/corda",
    python_reference_commit =
      "c02e06d50606bf93f23d8f2e6d6ade0e996ca70e",
    python_constructor_controls = c(
      n = "corda2_redundancies",
      penalty_factor = "corda2_penalty_factor",
      support = "corda2_support"
    ),
    python_fixed_constants = c(
      CI = 1.01, tflux = 1, UPPER = 1e6
    ),
    supported_scope = "met_prod = NULL",
    source_semantics = "exact",
    random_noise = FALSE
  )
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
