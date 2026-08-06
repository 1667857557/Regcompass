.rc_flatten_microcompass_vmax_cache <- function(grouped, expected_names) {
  expected_names <- as.character(expected_names)
  answer <- unlist(
    unname(grouped),
    recursive = FALSE,
    use.names = TRUE
  )
  observed_names <- names(answer)
  if (is.null(observed_names) ||
      length(answer) != length(expected_names) ||
      anyNA(observed_names) ||
      any(!nzchar(observed_names)) ||
      anyDuplicated(observed_names) ||
      !setequal(observed_names, expected_names)) {
    stop("The shared directional vmax cache is incomplete.", call. = FALSE)
  }
  answer[expected_names]
}

.rc_microcompass_worker_count <- function(parallel, BPPARAM, n_tasks) {
  if (!isTRUE(parallel) || n_tasks <= 1L) return(1L)
  if (!is.null(BPPARAM) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(BPPARAM, "BiocParallelParam")) {
    return(max(1L, min(
      as.integer(n_tasks),
      as.integer(BiocParallel::bpnworkers(BPPARAM))
    )))
  }
  available <- get0(
    "rc_available_workers", mode = "function", inherits = TRUE
  )
  workers <- if (is.function(available)) {
    available(default = 1L)
  } else {
    detected <- suppressWarnings(as.integer(
      parallel::detectCores(logical = TRUE)[[1L]]
    ))
    if (is.finite(detected) && detected > 1L) detected - 1L else 1L
  }
  max(1L, min(as.integer(n_tasks), as.integer(workers)))
}

.rc_microcompass_vmax_tasks <- function(model_keys, workers) {
  if (is.null(names(model_keys)) || any(!nzchar(names(model_keys)))) {
    stop("Directional vmax model keys require named target rows.",
         call. = FALSE)
  }
  unique_keys <- unique(unname(model_keys))
  workers <- max(1L, as.integer(workers[[1L]]))
  batches_per_model <- max(1L, ceiling(workers / length(unique_keys)))
  tasks <- list()
  cursor <- 0L

  for (model_index in seq_along(unique_keys)) {
    model_key <- unique_keys[[model_index]]
    selected_rows <- names(model_keys)[model_keys == model_key]
    n_batches <- min(length(selected_rows), batches_per_model)
    batch_id <- ceiling(seq_along(selected_rows) * n_batches /
                          length(selected_rows))
    groups <- split(selected_rows, batch_id)
    for (batch_index in seq_along(groups)) {
      cursor <- cursor + 1L
      tasks[[cursor]] <- list(
        model_key = model_key,
        row_ids = as.character(groups[[batch_index]])
      )
      names(tasks)[[cursor]] <- paste0(
        "model_", model_index, "__batch_", batch_index
      )
    }
  }
  tasks
}

.rc_build_microcompass_vmax_cache_core <- function(
    model_cache, mode, model_keys, solver, flux_threshold,
    parallel = TRUE, BPPARAM = NULL) {
  workers <- .rc_microcompass_worker_count(
    parallel = parallel,
    BPPARAM = BPPARAM,
    n_tasks = length(model_cache)
  )
  tasks <- .rc_microcompass_vmax_tasks(model_keys, workers)
  grouped <- rc_parallel_lapply(
    tasks,
    function(task) {
      selected_rows <- as.character(task$row_ids)
      first_entry <- model_cache[[selected_rows[[1L]]]]
      model <- .rc_load_microcompass_model(first_entry, mode)
      values <- lapply(selected_rows, function(row_id) {
        entry <- model_cache[[row_id]]
        rc_compass_vmax_directional(
          S = model$S,
          lb = model$lb,
          ub = model$ub,
          target_reaction = entry$reaction_id,
          direction = entry$target_direction,
          solver = solver,
          flux_threshold = flux_threshold
        )
      })
      names(values) <- selected_rows
      values
    },
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  answer <- .rc_flatten_microcompass_vmax_cache(
    grouped, names(model_cache)
  )
  attr(answer, "parallel_tasks") <- length(tasks)
  attr(answer, "parallel_workers") <- workers
  attr(answer, "parallel_scope") <-
    "directional_target_batches_within_shared_models"
  answer
}

.rc_compass_step2_align_penalties <- function(reactions, penalties) {
  if (!is.null(names(penalties))) {
    missing_penalties <- setdiff(reactions, names(penalties))
    if (length(missing_penalties)) {
      stop(
        "Reaction penalties are missing for: ",
        paste(utils::head(missing_penalties, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    penalties <- as.numeric(penalties[reactions])
  } else {
    penalties <- as.numeric(penalties)
  }
  if (length(penalties) != length(reactions) ||
      any(!is.finite(penalties)) || any(penalties < 0)) {
    stop("`penalties` must provide one finite non-negative value per reaction.",
         call. = FALSE)
  }
  penalties
}

.rc_compass_step2_prepare <- function(
    S, lb, ub, target_reaction, vmax_result,
    target_direction = c("forward", "reverse"),
    omega = 0.95, flux_threshold = 1e-8) {
  target_direction <- match.arg(target_direction)
  if (!is.numeric(omega) || length(omega) != 1L ||
      !is.finite(omega) || omega <= 0 || omega > 1) {
    stop("`omega` must be in (0, 1].", call. = FALSE)
  }
  if (!is.list(vmax_result) ||
      !all(c("feasible", "vmax", "status") %in% names(vmax_result))) {
    stop("`vmax_result` must be a directional Step 1 result.", call. = FALSE)
  }
  reactions <- colnames(S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions)) ||
      anyDuplicated(reactions)) {
    stop("`S` must have unique non-empty reaction IDs in colnames().",
         call. = FALSE)
  }
  if (!target_reaction %in% reactions) {
    stop("`target_reaction` is missing from the stoichiometric matrix.",
         call. = FALSE)
  }
  lb <- rc_align_bound(lb, reactions, default = -1000, name = "lb")
  ub <- rc_align_bound(ub, reactions, default = 1000, name = "ub")
  if (any(lb > ub)) {
    stop("Reaction lower bounds cannot exceed upper bounds.", call. = FALSE)
  }
  early <- if (!isTRUE(vmax_result$feasible)) {
    list(
      feasible = FALSE,
      penalty = NA_real_,
      vmax = as.numeric(vmax_result$vmax),
      solver_status = as.character(vmax_result$status),
      step1_status = as.character(vmax_result$status),
      step2_status = "not_run",
      solver_backend = "not_run",
      flux = numeric()
    )
  } else {
    NULL
  }
  if (!is.null(early)) {
    return(list(runnable = FALSE, reactions = reactions, result = early))
  }

  vmax <- as.numeric(vmax_result$vmax)
  if (length(vmax) != 1L || !is.finite(vmax) || vmax < flux_threshold) {
    stop("Cached directional vmax is not a positive feasible value.",
         call. = FALSE)
  }

  S <- .rc_as_dgCMatrix(S)
  n_reactions <- ncol(S)
  zero <- Matrix::Matrix(
    0, nrow = nrow(S), ncol = n_reactions, sparse = TRUE
  )
  mass_balance <- cbind(S, zero)
  positive <- Matrix::Matrix(
    0, nrow = n_reactions, ncol = 2L * n_reactions, sparse = TRUE
  )
  negative <- positive
  positive[cbind(seq_len(n_reactions), seq_len(n_reactions))] <- 1
  positive[cbind(
    seq_len(n_reactions), n_reactions + seq_len(n_reactions)
  )] <- -1
  negative[cbind(seq_len(n_reactions), seq_len(n_reactions))] <- -1
  negative[cbind(
    seq_len(n_reactions), n_reactions + seq_len(n_reactions)
  )] <- -1
  target <- Matrix::Matrix(
    0, nrow = 1, ncol = 2L * n_reactions, sparse = TRUE
  )
  target_index <- match(target_reaction, reactions)
  target[1, target_index] <- if (
    identical(target_direction, "forward")
  ) 1 else -1
  A <- rbind(mass_balance, positive, negative, target)
  lhs <- c(
    rep(0, nrow(S)),
    rep(-Inf, 2L * n_reactions),
    omega * vmax
  )
  rhs <- c(
    rep(0, nrow(S)),
    rep(0, 2L * n_reactions),
    Inf
  )
  auxiliary_upper <- pmax(abs(lb), abs(ub))

  list(
    runnable = TRUE,
    reactions = reactions,
    template = list(
      A = A,
      lhs = lhs,
      rhs = rhs,
      lb = c(lb, rep(0, n_reactions)),
      ub = c(ub, auxiliary_upper),
      n_reactions = n_reactions,
      reactions = reactions,
      vmax = vmax,
      step1_status = as.character(vmax_result$status),
      target_reaction = target_reaction,
      target_direction = target_direction
    )
  )
}

.rc_microcompass_highs_api_available <- function() {
  if (!requireNamespace("highs", quietly = TRUE)) return(FALSE)
  required <- c(
    "highs_model", "hi_new_solver", "hi_solver_set_objective",
    "hi_solver_run", "hi_solver_status_message",
    "hi_solver_get_solution", "hi_solver_info", "hi_solver_set_option"
  )
  all(required %in% getNamespaceExports("highs"))
}

.rc_microcompass_highs_call <- function(name, ...) {
  do.call(getExportedValue("highs", name), list(...))
}

.rc_compass_step2_release_engine <- function(engine) {
  if (!is.list(engine) || is.null(engine$pointer)) return(invisible(engine))
  if (requireNamespace("highs", quietly = TRUE)) {
    exports <- getNamespaceExports("highs")
    if ("hi_solver_clear" %in% exports) {
      try(.rc_microcompass_highs_call(
        "hi_solver_clear", engine$pointer
      ), silent = TRUE)
    } else if ("hi_solver_clear_model" %in% exports) {
      try(.rc_microcompass_highs_call(
        "hi_solver_clear_model", engine$pointer
      ), silent = TRUE)
    }
  }
  engine$pointer <- NULL
  invisible(engine)
}

.rc_compass_step2_new_engine <- function(template, solver) {
  solver <- match.arg(solver, c("highs", "gurobi", "glpk"))
  n_reactions <- as.integer(template$n_reactions)
  engine <- list(
    type = "one_shot",
    pointer = NULL,
    template = template,
    solver = solver,
    current_penalties = rep(0, n_reactions),
    n_solves = 0L,
    n_objective_updates = 0L,
    n_fallback = 0L,
    persistent_disabled = FALSE
  )
  if (!identical(solver, "highs") ||
      !.rc_microcompass_highs_api_available()) {
    return(engine)
  }

  persistent <- tryCatch({
    model <- .rc_microcompass_highs_call(
      "highs_model",
      L = rep(0, 2L * n_reactions),
      lower = template$lb,
      upper = template$ub,
      A = template$A,
      lhs = template$lhs,
      rhs = template$rhs,
      maximum = FALSE
    )
    pointer <- .rc_microcompass_highs_call("hi_new_solver", model)
    .rc_microcompass_highs_call(
      "hi_solver_set_option", pointer, "output_flag", FALSE
    )
    .rc_microcompass_highs_call(
      "hi_solver_set_option", pointer, "threads", 1L
    )
    .rc_microcompass_highs_call(
      "hi_solver_set_option", pointer, "solver", "simplex"
    )
    pointer
  }, error = function(e) e)

  if (inherits(persistent, "error")) {
    engine$persistent_disabled <- TRUE
    engine$persistent_message <- conditionMessage(persistent)
    return(engine)
  }
  engine$type <- "highs_persistent_cpp"
  engine$pointer <- persistent
  engine
}

.rc_compass_step2_one_shot <- function(engine, penalties) {
  template <- engine$template
  n_reactions <- template$n_reactions
  answer <- rc_solve_lp(
    obj = c(rep(0, n_reactions), penalties),
    A = template$A,
    lhs = template$lhs,
    rhs = template$rhs,
    lb = template$lb,
    ub = template$ub,
    solver = engine$solver
  )
  answer$backend <- paste0("one_shot_", engine$solver)
  answer
}

.rc_compass_step2_engine_solve <- function(engine, penalties) {
  penalties <- .rc_compass_step2_align_penalties(
    engine$template$reactions, penalties
  )
  engine$n_solves <- engine$n_solves + 1L
  if (!identical(engine$type, "highs_persistent_cpp") ||
      isTRUE(engine$persistent_disabled)) {
    return(list(
      engine = engine,
      answer = .rc_compass_step2_one_shot(engine, penalties)
    ))
  }

  changed <- which(engine$current_penalties != penalties)
  persistent <- tryCatch({
    if (length(changed)) {
      .rc_microcompass_highs_call(
        "hi_solver_set_objective", engine$pointer,
        index = as.integer(engine$template$n_reactions + changed - 1L),
        coeff = penalties[changed]
      )
    }
    .rc_microcompass_highs_call("hi_solver_run", engine$pointer)
    status_message <- .rc_microcompass_highs_call(
      "hi_solver_status_message", engine$pointer
    )
    status <- .rc_lp_status(status_message)
    solution <- if (identical(status, "optimal")) {
      as.numeric(.rc_microcompass_highs_call(
        "hi_solver_get_solution", engine$pointer
      )$col_value)
    } else {
      numeric()
    }
    info <- tryCatch(
      .rc_microcompass_highs_call("hi_solver_info", engine$pointer),
      error = function(e) list()
    )
    objective <- as.numeric(
      info$objective_function_value %||% NA_real_
    )
    if (!is.finite(objective) &&
        length(solution) == 2L * engine$template$n_reactions) {
      objective <- sum(
        penalties * solution[
          engine$template$n_reactions + seq_len(
            engine$template$n_reactions
          )
        ]
      )
    }
    list(
      status = status,
      solution = solution,
      objective = objective,
      solver = "highs",
      backend = "highs_persistent_cpp_basis_reuse",
      solver_message = as.character(status_message)
    )
  }, error = function(e) e)

  if (!inherits(persistent, "error")) {
    engine$current_penalties <- penalties
    engine$n_objective_updates <-
      engine$n_objective_updates + length(changed)
    return(list(engine = engine, answer = persistent))
  }

  engine$n_fallback <- engine$n_fallback + 1L
  engine$persistent_message <- conditionMessage(persistent)
  engine <- .rc_compass_step2_release_engine(engine)
  engine$type <- "one_shot"
  engine$persistent_disabled <- TRUE
  fallback <- .rc_compass_step2_one_shot(engine, penalties)
  fallback$backend <- "highs_persistent_failed_one_shot_fallback"
  fallback$solver_message <- paste(
    engine$persistent_message,
    fallback$solver_message %||% ""
  )
  list(engine = engine, answer = fallback)
}

.rc_compass_step2_engine_metrics <- function(engine) {
  if (!is.list(engine)) {
    return(list(
      engine = "not_run", n_solves = 0L,
      n_objective_updates = 0L, n_fallback = 0L
    ))
  }
  list(
    engine = as.character(engine$type %||% "one_shot"),
    n_solves = as.integer(engine$n_solves %||% 0L),
    n_objective_updates = as.integer(
      engine$n_objective_updates %||% 0L
    ),
    n_fallback = as.integer(engine$n_fallback %||% 0L)
  )
}

.rc_compass_step2_result <- function(template, answer) {
  n_reactions <- template$n_reactions
  if (!identical(answer$status, "optimal") ||
      length(answer$solution) != 2L * n_reactions) {
    return(list(
      feasible = FALSE,
      penalty = NA_real_,
      vmax = template$vmax,
      solver_status = answer$status,
      step1_status = template$step1_status,
      step2_status = answer$status,
      solver_backend = answer$backend %||% "unknown",
      flux = numeric()
    ))
  }
  flux <- answer$solution[seq_len(n_reactions)]
  names(flux) <- template$reactions
  list(
    feasible = TRUE,
    penalty = max(0, as.numeric(answer$objective)),
    vmax = template$vmax,
    solver_status = answer$status,
    step1_status = template$step1_status,
    step2_status = answer$status,
    solver_backend = answer$backend %||% "unknown",
    flux = flux
  )
}

.rc_compass_step2_from_vmax_directional <- function(
    S, lb, ub, target_reaction, penalties, vmax_result,
    target_direction = c("forward", "reverse"),
    omega = 0.95,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8) {
  target_direction <- match.arg(target_direction)
  solver <- match.arg(solver)
  prepared <- .rc_compass_step2_prepare(
    S = S,
    lb = lb,
    ub = ub,
    target_reaction = target_reaction,
    vmax_result = vmax_result,
    target_direction = target_direction,
    omega = omega,
    flux_threshold = flux_threshold
  )
  penalties <- .rc_compass_step2_align_penalties(
    prepared$reactions, penalties
  )
  if (!isTRUE(prepared$runnable)) return(prepared$result)

  engine <- .rc_compass_step2_new_engine(prepared$template, solver)
  on.exit(.rc_compass_step2_release_engine(engine), add = TRUE)
  solved <- .rc_compass_step2_engine_solve(engine, penalties)
  engine <- solved$engine
  .rc_compass_step2_result(prepared$template, solved$answer)
}

# Progress-aware entry point; the algorithm remains in the core above.
.rc_build_microcompass_vmax_cache <- function(...) {
  args <- list(...)
  contexts <- .rc_layer2_model_contexts(
    args$model_cache,
    mode = as.character(args$mode %||% "meta_module_gem")
  )
  parts_dir <- .rc_layer2_cache_progress_dir(args$model_cache)
  run_kind <- .rc_layer2_progress_state$run_kind %||% "primary"
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "directional_vmax_start", 1L, 4L,
      detail = paste0("directional_targets=", item$n_targets),
      scope = "scoring", run_kind = run_kind, parts_dir = parts_dir
    )
  }
  answer <- do.call(
    .rc_build_microcompass_vmax_cache_core,
    args
  )
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "directional_vmax_complete", 2L, 4L,
      detail = paste0("directional_targets=", item$n_targets),
      scope = "scoring", run_kind = run_kind,
      status = "complete", parts_dir = parts_dir
    )
    .rc_layer2_task_event(
      item$context, "penalty_step2_start", 3L, 4L,
      "scoring matching metacells with the cached directional vmax",
      scope = "scoring", run_kind = run_kind, parts_dir = parts_dir
    )
  }
  if (identical(run_kind, "primary")) {
    .rc_layer2_overall_event(
      "primary_vmax_complete", 4L,
      detail = paste0(
        "directional target batches completed; tasks=",
        attr(answer, "parallel_tasks") %||% length(answer)
      )
    )
  }
  answer
}
