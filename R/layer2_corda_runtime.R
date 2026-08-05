# Native runtime dispatch for exact pinned Python CORDA2 Layer 2 completion.

.rc_layer2_completion_context <- new.env(parent = emptyenv())
.rc_layer2_completion_context$active <- FALSE
.rc_layer2_completion_context$model_completion <- "fastcore"
.rc_layer2_completion_context$reaction_evidence <- NULL
.rc_layer2_completion_context$corda_options <- NULL

.rc_regcompass_step_layer2_completion_base <- rc_regcompass_step_layer2
.rc_build_celltype_medium_union_gem_cache_fastcore <-
  .rc_build_celltype_medium_union_gem_cache

.rc_is_corda2_options <- function(options) {
  is.list(options) && identical(
    as.character(options$algorithm %||% ""),
    "resendislab_python_CORDA2_c02e06d_exact_semantics"
  )
}

# Persistent HiGHS updates preserve one solver instance across the Python-order
# target loop. One-shot solving remains a runtime fallback only.
.rc_corda_engine_solve <- function(
    engine, objective, lower, upper) {
  engine$n_solves <- engine$n_solves + 1L
  if (!identical(engine$type, "highs_persistent_cpp")) {
    return(list(
      answer = .rc_corda_one_shot_solve(engine, objective, lower, upper),
      engine = engine
    ))
  }
  persistent <- tryCatch({
    index <- seq_along(objective) - 1L
    .rc_corda_highs_call(
      "hi_solver_set_objective", engine$pointer,
      index = index, coeff = as.numeric(objective)
    )
    .rc_corda_highs_call(
      "hi_solver_set_variable_bounds", engine$pointer,
      index = index,
      lower = as.numeric(lower),
      upper = as.numeric(upper)
    )
    .rc_corda_highs_call("hi_solver_run", engine$pointer)
    status_message <- .rc_corda_highs_call(
      "hi_solver_status_message", engine$pointer
    )
    status <- .rc_lp_status(status_message)
    solution <- if (identical(status, "optimal")) {
      as.numeric(.rc_corda_highs_call(
        "hi_solver_get_solution", engine$pointer
      )$col_value)
    } else {
      numeric()
    }
    info <- tryCatch(
      .rc_corda_highs_call("hi_solver_info", engine$pointer),
      error = function(e) list()
    )
    list(
      status = status,
      solution = solution,
      objective = as.numeric(
        info$objective_function_value %||% NA_real_
      ),
      backend = "highs_persistent_cpp_basis_reuse",
      solver_message = as.character(status_message)
    )
  }, error = function(e) e)
  if (!inherits(persistent, "error")) {
    return(list(answer = persistent, engine = engine))
  }
  failure_message <- conditionMessage(persistent)
  engine$n_fallback <- engine$n_fallback + 1L
  fallback <- .rc_corda_one_shot_solve(engine, objective, lower, upper)
  fallback$backend <- "highs_persistent_failed_one_shot_fallback"
  fallback$solver_message <- paste(
    failure_message, fallback$solver_message
  )
  list(answer = fallback, engine = engine)
}

.rc_corda_pool_workers <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE)) return(1L)
  if (!is.null(BPPARAM) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(BPPARAM, "BiocParallelParam")) {
    return(max(1L, BiocParallel::bpnworkers(BPPARAM)))
  }
  config <- rc_parallel_config(workers = NULL, backend = "auto")
  max(1L, config$workers)
}
