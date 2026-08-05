# Public parameter contract for exact Python CORDA2 execution semantics.

.rc_corda2_constructor_args <- function(model_params) {
  exact <- model_params$corda2_args %||% list()
  if (!is.list(exact)) {
    stop("`corda2_args` must be a named list.", call. = FALSE)
  }
  if (length(exact) && (is.null(names(exact)) || any(!nzchar(names(exact))))) {
    stop("`corda2_args` must be a named list.", call. = FALSE)
  }
  allowed <- c("met_prod", "n", "penalty_factor", "support")
  unknown <- setdiff(names(exact), allowed)
  if (length(unknown)) {
    stop(
      "Unknown Python CORDA constructor parameter(s) in `corda2_args`: ",
      paste(unknown, collapse = ", "),
      ". Allowed names are met_prod, n, penalty_factor and support.",
      call. = FALSE
    )
  }

  same_value <- function(left, right) {
    isTRUE(all.equal(left, right, check.attributes = FALSE))
  }
  resolve <- function(name, legacy, default) {
    has_exact <- name %in% names(exact)
    has_legacy <- !is.null(legacy)
    if (has_exact && has_legacy && !same_value(exact[[name]], legacy)) {
      stop(
        "Conflicting CORDA2 values for Python parameter `", name,
        "` and its legacy RegCompass alias.",
        call. = FALSE
      )
    }
    if (has_exact) return(exact[[name]])
    if (has_legacy) return(legacy)
    default
  }

  met_prod <- if ("met_prod" %in% names(exact)) exact$met_prod else NULL
  if (!is.null(met_prod)) {
    stop(
      "RegCompass exact CORDA2 currently supports the Python constructor only ",
      "for `met_prod = NULL`; non-NULL mock metabolic targets are not silently ",
      "approximated.",
      call. = FALSE
    )
  }

  n <- .rc_corda_integer(
    resolve(
      "n",
      model_params$corda2_redundancies %||% NULL,
      3L
    ),
    "corda2_args$n",
    0L
  )
  penalty_factor <- .rc_corda_scalar(
    resolve(
      "penalty_factor",
      model_params$corda2_penalty_factor %||%
        model_params$corda_penalty_factor %||% NULL,
      100
    ),
    "corda2_args$penalty_factor",
    -Inf,
    Inf
  )
  support <- .rc_corda_integer(
    resolve(
      "support",
      model_params$corda2_support %||%
        model_params$corda_support %||% NULL,
      5L
    ),
    "corda2_args$support",
    0L
  )

  list(
    met_prod = NULL,
    n = n,
    penalty_factor = penalty_factor,
    support = support
  )
}

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
      ". The source fixes CI=1.01, tflux=1 and UPPER=1e6; its constructor ",
      "controls are met_prod, n, penalty_factor and support.",
      call. = FALSE
    )
  }

  constructor <- .rc_corda2_constructor_args(model_params)
  scalar <- .rc_corda_scalar
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
    met_prod = constructor$met_prod,
    n = constructor$n,
    # Internal compatibility alias. The canonical source parameter is `n`.
    redundancies = constructor$n,
    penalty_factor = constructor$penalty_factor,
    support = constructor$support,
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
    python_constructor_signature = c(
      "model", "confidence", "met_prod", "n", "penalty_factor", "support"
    ),
    python_constructor_controls = c(
      met_prod = "corda2_args$met_prod",
      n = "corda2_args$n",
      penalty_factor = "corda2_args$penalty_factor",
      support = "corda2_args$support"
    ),
    python_constructor_defaults = list(
      met_prod = NULL,
      n = 3L,
      penalty_factor = 100,
      support = 5L
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
