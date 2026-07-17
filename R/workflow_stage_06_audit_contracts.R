# Workflow stage 6: harden effective contracts identified by repository audit.
#
# This stage keeps the existing staged architecture while ensuring that the
# parameters recorded in stratum artifacts and final results are the parameters
# used during global reaction-capacity calculation. It also restores the
# advertised Q95 diagnostic arguments without allowing Q95 scaling to replace
# absolute LP evidence.

.rc_required_previous_run_regcompass_stratum <- .rc_run_regcompass_stratum
.rc_required_previous_run_regcompass_audit <- rc_run_regcompass
.rc_required_previous_metacell_logcpm_audit <- .rc_metacell_logcpm

.rc_resolve_workflow_or_method <- function(
    layer1_args = list(), strict_biological_defaults = NULL) {
  if (!is.list(layer1_args)) {
    stop("`layer1_args` must be a list.", call. = FALSE)
  }
  if (!is.null(layer1_args$or_method)) {
    return(match.arg(
      as.character(layer1_args$or_method),
      c("max", "sum_sqrtK", "prob_or", "sum")
    ))
  }
  if (isTRUE(strict_biological_defaults)) return("max")
  strict_defaults <- identical(
    layer1_args$promiscuity_mode %||% "sqrt",
    "none"
  ) && identical(
    layer1_args$and_method %||% "boltzmann",
    "min"
  )
  if (strict_defaults) "max" else "sum_sqrtK"
}

.rc_workflow_or_method_source <- function(
    layer1_args = list(), strict_biological_defaults = NULL) {
  if (!is.null(layer1_args$or_method)) return("explicit_layer1_args")
  if (isTRUE(strict_biological_defaults)) {
    return("strict_biological_default")
  }
  strict_defaults <- identical(
    layer1_args$promiscuity_mode %||% "sqrt",
    "none"
  ) && identical(
    layer1_args$and_method %||% "boltzmann",
    "min"
  )
  if (strict_defaults) {
    "strict_biological_default"
  } else {
    "legacy_sensitivity_default"
  }
}

.rc_finalize_stratum_capacity_params <- function(artifact_file,
                                                   layer1_args = list()) {
  if (!is.character(artifact_file) || length(artifact_file) != 1L ||
      is.na(artifact_file) || !nzchar(artifact_file) ||
      !file.exists(artifact_file)) {
    stop("A valid stratum artifact file is required.", call. = FALSE)
  }
  artifact <- readRDS(artifact_file)
  if (!is.list(artifact$capacity_params)) {
    stop("The stratum artifact lacks `capacity_params`.", call. = FALSE)
  }
  or_method <- .rc_resolve_workflow_or_method(layer1_args)
  artifact$capacity_params$or_method <- or_method
  artifact$capacity_params$or_method_source <-
    .rc_workflow_or_method_source(layer1_args)
  saveRDS(artifact, artifact_file)
  invisible(or_method)
}

# `global_workflow.R` historically rejected `layer1_args$or_method` and wrote
# `sum_sqrtK` into every stratum artifact. The global merge subsequently used
# that stored value even when strict biological defaults were requested. Accept
# the argument here, remove it before legacy validation, and correct the artifact
# before the global barrier reads it.
.rc_run_regcompass_stratum <- function(...) {
  args <- list(...)
  layer1_args <- args$layer1_args %||% list()
  if (!is.list(layer1_args)) {
    stop("`layer1_args` must be a list.", call. = FALSE)
  }
  args$layer1_args <- layer1_args
  args$layer1_args$or_method <- NULL
  answer <- do.call(.rc_required_previous_run_regcompass_stratum, args)
  if (identical(answer$status, "ok")) {
    .rc_finalize_stratum_capacity_params(
      answer$artifact_file,
      layer1_args = layer1_args
    )
  }
  answer
}

# Earlier result wrappers rewrote strict-mode metadata after computation. Apply
# the resolved value once more at the public boundary and resave the canonical
# result files so in-memory and serialized provenance agree.
rc_run_regcompass <- function(
    object, gem, outdir, pfm, genome,
    fragment_files = NULL,
    sample_col = "sample_id",
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
    upstream_workers = NULL,
    layer2_workers = NULL,
    parallel_backend = c("auto", "serial", "snow", "multicore"),
    strict_biological_defaults = TRUE,
    inference_unit = c("sample_celltype", "metacell")) {
  answer <- .rc_required_previous_run_regcompass_audit(
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
    parallel_backend = parallel_backend,
    strict_biological_defaults = strict_biological_defaults,
    inference_unit = inference_unit
  )

  or_method <- .rc_resolve_workflow_or_method(
    layer1_args,
    strict_biological_defaults = strict_biological_defaults
  )
  or_method_source <- .rc_workflow_or_method_source(
    layer1_args,
    strict_biological_defaults = strict_biological_defaults
  )
  if (is.list(answer$layer1$capacity_params)) {
    answer$layer1$capacity_params$or_method <- or_method
    answer$layer1$capacity_params$or_method_source <- or_method_source
  }
  answer$params$gpr_or_method <- or_method
  answer$params$gpr_or_method_source <- or_method_source

  if (is.character(outdir) && length(outdir) == 1L &&
      !is.na(outdir) && nzchar(outdir) && dir.exists(outdir)) {
    saveRDS(answer, file.path(outdir, "regcompass_global_metacell_result.rds"))
    saveRDS(answer, file.path(outdir, "regcompass_result.rds"))
  }
  answer
}

# Sparse matrix multiplication can drop column names. Preserve the complete
# metacell identifiers because subsequent alignment and diagnostics depend on
# exact dimnames.
.rc_metacell_logcpm <- function(counts, scale_factor = 1e6,
                                 library_size = NULL) {
  input_dimnames <- dimnames(counts)
  answer <- .rc_required_previous_metacell_logcpm_audit(
    counts = counts,
    scale_factor = scale_factor,
    library_size = library_size
  )
  dimnames(answer) <- input_dimnames
  answer
}

# Q95 remains a within-reaction diagnostic. The LP continues to use bounded
# absolute evidence in C_rel/C_abs, while n0, stratification and bootstrap
# arguments now affect the diagnostic output as their signature states.
rc_q95_calibrate <- function(C_raw, eps = 1e-6, bootstrap = TRUE,
                             B = 500, BPPARAM = NULL, n0 = 80,
                             unit_meta = NULL, stratum_col = NULL) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) {
    stop("`C_raw` must have reaction rownames and unit colnames.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps <= 0) {
    stop("`eps` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(n0) || length(n0) != 1L ||
      !is.finite(n0) || n0 < 0) {
    stop("`n0` must be one finite non-negative number.", call. = FALSE)
  }
  if (!is.logical(bootstrap) || length(bootstrap) != 1L || is.na(bootstrap)) {
    stop("`bootstrap` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(bootstrap) &&
      (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1)) {
    stop("`B` must be one positive finite number.", call. = FALSE)
  }

  diagnostic <- rc_q95_shrink(
    C_raw,
    unit_meta = unit_meta,
    stratum_col = stratum_col,
    q = 0.95,
    n0 = n0,
    eps = eps
  )
  relative <- .rc_clamp01(diagnostic$C_rel)
  finite_range <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    diff(range(x))
  })
  noninformative <- !is.finite(finite_range) | finite_range <= eps
  if (any(noninformative)) relative[noninformative, ] <- NA_real_

  absolute <- .rc_clamp01(C_raw)
  all_missing <- rowSums(is.finite(C_raw)) == 0L
  if (any(all_missing)) absolute[all_missing, ] <- NA_real_
  all_zero <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) > 0L && max(x) <= eps
  })
  if (any(all_zero)) absolute[all_zero, ] <- 0

  Q <- diagnostic$Q
  names(Q)[names(Q) == "q_shrink"] <- "q_value"
  reaction_index <- match(Q$reaction_id, rownames(C_raw))
  n_finite_global <- rowSums(is.finite(C_raw))
  minimum <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  })
  maximum <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  })
  informative <- rowSums(is.finite(relative)) > 0L

  Q$quantile_used <- 0.95
  Q$n_finite <- Q$n
  Q$n_finite_global <- as.integer(n_finite_global[reaction_index])
  Q$low_n_flag <- Q$n_finite < 20L
  Q$all_zero_reaction_flag <-
    is.finite(maximum[reaction_index]) & maximum[reaction_index] <= eps
  Q$constant_reaction_flag <-
    is.finite(minimum[reaction_index]) &
    is.finite(maximum[reaction_index]) &
    abs(maximum[reaction_index] - minimum[reaction_index]) <= eps
  Q$raw_out_of_unit_interval_flag <- vapply(
    reaction_index,
    function(i) any(is.finite(C_raw[i, ]) &
                      (C_raw[i, ] < -eps | C_raw[i, ] > 1 + eps)),
    logical(1)
  )
  Q$relative_capacity_informative <- informative[reaction_index]
  Q$sample_balanced <- FALSE
  Q$calibration_role <- "diagnostic_only_not_lp_capacity"

  if (isTRUE(bootstrap)) {
    boot <- rc_q95_bootstrap_diagnostics(
      C_raw,
      Q,
      unit_meta = unit_meta,
      stratum_col = stratum_col,
      B = as.integer(B),
      BPPARAM = BPPARAM
    )
    Q$q95_bootstrap <- boot[, "q95"]
    Q$q95_ci_low <- boot[, "ci_low"]
    Q$q95_ci_high <- boot[, "ci_high"]
    Q$q95_ci_width <- boot[, "width"]
    Q$q95_unstable_flag <- is.finite(Q$q95_ci_width) &
      (Q$q95_ci_width / pmax(Q$q_value, eps)) > 0.5
  }

  list(
    C_rel = absolute,
    C_abs = absolute,
    C_within_reaction_relative = relative,
    Q = Q
  )
}
