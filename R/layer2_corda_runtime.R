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

# Retain the exact constructor from layer2_corda_lp.R and add only runtime
# bookkeeping. The mathematical model, bounds, tolerance and solver options are
# unchanged.
if (!exists(".rc_corda_new_lp_engine_runtime_base", inherits = FALSE)) {
  .rc_corda_new_lp_engine_runtime_base <- .rc_corda_new_lp_engine
}

.rc_corda_new_lp_engine <- function(split, solver, time_limit) {
  engine <- .rc_corda_new_lp_engine_runtime_base(
    split = split,
    solver = solver,
    time_limit = time_limit
  )
  n_variables <- ncol(split$S)
  engine$current_objective <- rep(0, n_variables)
  engine$current_lower <- as.numeric(split$lb)
  engine$current_upper <- as.numeric(split$ub)
  engine$n_objective_coeff_updates <- 0L
  engine$n_bound_index_updates <- 0L
  engine$n_sparse_update_calls <- 0L
  engine$n_full_vector_numeric_values <- 0
  engine$n_transmitted_numeric_values <- 0
  engine$n_full_vector_numeric_values_avoided <- 0
  engine$persistent_disabled <- FALSE
  engine$released <- FALSE
  engine
}

.rc_corda_changed_indices <- function(current, requested, label) {
  current <- as.numeric(current)
  requested <- as.numeric(requested)
  if (length(current) != length(requested)) {
    stop("CORDA2 ", label, " length changed during persistent solving.",
         call. = FALSE)
  }
  which(
    (current != requested) |
      (is.na(current) & !is.na(requested)) |
      (!is.na(current) & is.na(requested))
  )
}

# Explicitly clear the native HiGHS model as soon as one independent CORDA2
# reconstruction finishes. The external pointer still has its own finalizer;
# this call releases the large model and solver state before worker reuse.
.rc_corda_release_lp_engine <- function(engine) {
  if (!is.list(engine) || isTRUE(engine$released)) return(engine)
  if (!is.null(engine$pointer) && requireNamespace("highs", quietly = TRUE)) {
    exports <- getNamespaceExports("highs")
    if ("hi_solver_clear" %in% exports) {
      try(.rc_corda_highs_call("hi_solver_clear", engine$pointer), silent = TRUE)
    } else if ("hi_solver_clear_model" %in% exports) {
      try(.rc_corda_highs_call(
        "hi_solver_clear_model", engine$pointer
      ), silent = TRUE)
    }
  }
  engine$pointer <- NULL
  engine$released <- TRUE
  engine
}

.rc_corda_execution_metrics <- function(engine) {
  n_variables <- if (is.list(engine$split) && !is.null(engine$split$S)) {
    ncol(engine$split$S)
  } else {
    0L
  }
  full_values <- as.numeric(engine$n_full_vector_numeric_values %||% 0)
  transmitted <- as.numeric(engine$n_transmitted_numeric_values %||% 0)
  list(
    n_variables = as.integer(n_variables),
    n_solves = as.integer(engine$n_solves %||% 0L),
    n_fallback = as.integer(engine$n_fallback %||% 0L),
    n_objective_coeff_updates = as.integer(
      engine$n_objective_coeff_updates %||% 0L
    ),
    n_bound_index_updates = as.integer(
      engine$n_bound_index_updates %||% 0L
    ),
    n_sparse_update_calls = as.integer(
      engine$n_sparse_update_calls %||% 0L
    ),
    n_full_vector_numeric_values = full_values,
    n_transmitted_numeric_values = transmitted,
    n_full_vector_numeric_values_avoided = as.numeric(
      engine$n_full_vector_numeric_values_avoided %||% 0
    ),
    transmitted_fraction_of_full = if (full_values > 0) {
      transmitted / full_values
    } else {
      NA_real_
    },
    persistent_solver = identical(engine$type, "highs_persistent_cpp"),
    persistent_disabled = isTRUE(engine$persistent_disabled),
    release_policy = "explicit_native_clear_on_reconstruction_exit"
  )
}

# Persistent HiGHS updates preserve one solver instance across the Python-order
# target loop. Only coefficients and bounds that differ from the current native
# state are transmitted. This is algebraically identical to complete-vector
# replacement and preserves the same simplex basis and serial mutation order.
.rc_corda_engine_solve <- function(
    engine, objective, lower, upper) {
  objective <- as.numeric(objective)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  n_variables <- length(objective)
  if (length(lower) != n_variables || length(upper) != n_variables) {
    stop("CORDA2 objective and bound vectors must have equal length.",
         call. = FALSE)
  }
  engine$n_solves <- engine$n_solves + 1L
  engine$n_full_vector_numeric_values <-
    engine$n_full_vector_numeric_values + 3 * n_variables

  if (!identical(engine$type, "highs_persistent_cpp") ||
      isTRUE(engine$persistent_disabled)) {
    engine$n_transmitted_numeric_values <-
      engine$n_transmitted_numeric_values + 3 * n_variables
    return(list(
      answer = .rc_corda_one_shot_solve(engine, objective, lower, upper),
      engine = engine
    ))
  }

  objective_index <- .rc_corda_changed_indices(
    engine$current_objective, objective, "objective"
  )
  bound_index <- union(
    .rc_corda_changed_indices(engine$current_lower, lower, "lower bound"),
    .rc_corda_changed_indices(engine$current_upper, upper, "upper bound")
  )
  transmitted <- length(objective_index) + 2L * length(bound_index)

  persistent <- tryCatch({
    if (length(objective_index)) {
      .rc_corda_highs_call(
        "hi_solver_set_objective", engine$pointer,
        index = as.integer(objective_index - 1L),
        coeff = objective[objective_index]
      )
    }
    if (length(bound_index)) {
      .rc_corda_highs_call(
        "hi_solver_set_variable_bounds", engine$pointer,
        index = as.integer(bound_index - 1L),
        lower = lower[bound_index],
        upper = upper[bound_index]
      )
    }
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
      backend = "highs_persistent_cpp_sparse_delta_basis_reuse",
      solver_message = as.character(status_message)
    )
  }, error = function(e) e)

  if (!inherits(persistent, "error")) {
    engine$current_objective <- objective
    engine$current_lower <- lower
    engine$current_upper <- upper
    engine$n_objective_coeff_updates <-
      engine$n_objective_coeff_updates + length(objective_index)
    engine$n_bound_index_updates <-
      engine$n_bound_index_updates + length(bound_index)
    engine$n_sparse_update_calls <- engine$n_sparse_update_calls + 1L
    engine$n_transmitted_numeric_values <-
      engine$n_transmitted_numeric_values + transmitted
    engine$n_full_vector_numeric_values_avoided <-
      engine$n_full_vector_numeric_values_avoided +
      (3 * n_variables - transmitted)
    return(list(answer = persistent, engine = engine))
  }

  failure_message <- conditionMessage(persistent)
  engine$n_fallback <- engine$n_fallback + 1L
  engine <- .rc_corda_release_lp_engine(engine)
  engine$type <- "one_shot"
  engine$persistent_disabled <- TRUE
  engine$current_objective <- objective
  engine$current_lower <- lower
  engine$current_upper <- upper
  engine$n_transmitted_numeric_values <-
    engine$n_transmitted_numeric_values + 3 * n_variables
  fallback <- .rc_corda_one_shot_solve(engine, objective, lower, upper)
  fallback$backend <- "highs_persistent_failed_one_shot_fallback"
  fallback$solver_message <- paste(
    failure_message, fallback$solver_message
  )
  list(answer = fallback, engine = engine)
}

.rc_corda_task_bpparam <- function() {
  .rc_layer2_task_bpparam()
}

.rc_corda_worker_count <- function(BPPARAM, n_tasks) {
  # Exact CORDA2 reconstruction is serial within one model. Parallelism is
  # permitted only across independent cell-type x medium model instances.
  1L
}

.rc_corda_tune_task_bpparam <- function(BPPARAM, n_tasks) {
  n_tasks <- max(1L, as.integer(n_tasks[[1L]]))
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM) || n_tasks <= 1L) {
    return(BPPARAM)
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) {
    return(BPPARAM)
  }
  setter <- get0(
    "bptasks<-", envir = asNamespace("BiocParallel"),
    mode = "function", inherits = FALSE
  )
  if (is.function(setter)) {
    tuned <- tryCatch(
      setter(BPPARAM, n_tasks),
      error = function(e) BPPARAM
    )
  } else {
    tuned <- BPPARAM
  }
  attr(tuned, "regcompass_corda2_dynamic_tasks") <- n_tasks
  tuned
}

.rc_corda_should_outer_parallel <- function(n_tasks, pool_workers) {
  as.integer(n_tasks) > 1L && as.integer(pool_workers) > 1L
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
