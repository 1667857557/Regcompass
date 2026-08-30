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

.rc_microcompass_vmax_cache_contract <- function(
    model_cache, mode, solver, flux_threshold) {
  if (!is.list(model_cache) || !length(model_cache) ||
      is.null(names(model_cache)) || anyNA(names(model_cache)) ||
      any(!nzchar(names(model_cache))) || anyDuplicated(names(model_cache))) {
    stop(
      "Directional vmax cache contracts require a uniquely named model cache.",
      call. = FALSE
    )
  }
  mode <- match.arg(
    as.character(mode), c("full_gem", "meta_module_gem")
  )
  solver <- match.arg(
    as.character(solver), c("highs", "gurobi", "glpk")
  )
  if (!is.numeric(flux_threshold) || length(flux_threshold) != 1L ||
      !is.finite(flux_threshold) || flux_threshold <= 0) {
    stop("`flux_threshold` must be one positive finite number.",
         call. = FALSE)
  }
  scalar_text <- function(x, label, allow_empty = FALSE) {
    value <- as.character(x %||% "")
    if (length(value) != 1L || is.na(value) ||
        (!allow_empty && !nzchar(value))) {
      stop("Directional vmax cache entry has invalid `", label, "`.",
           call. = FALSE)
    }
    value
  }
  rows <- lapply(names(model_cache), function(row_id) {
    entry <- model_cache[[row_id]]
    model_file <- scalar_text(entry$file, "file", allow_empty = TRUE)
    checksum <- as.character(entry$file_checksum %||% NA_character_)
    if (length(checksum) != 1L || is.na(checksum) || !nzchar(checksum)) {
      checksum <- if (nzchar(model_file) && file.exists(model_file)) {
        unname(tools::md5sum(model_file)[[1L]])
      } else {
        NA_character_
      }
    }
    data.frame(
      row_id = row_id,
      model_file = model_file,
      model_file_checksum = checksum,
      cell_type = if (identical(mode, "meta_module_gem")) {
        scalar_text(entry$cell_type, "cell_type")
      } else {
        ""
      },
      reaction_id = scalar_text(entry$reaction_id, "reaction_id"),
      target_direction = scalar_text(
        entry$target_direction, "target_direction"
      ),
      medium_scenario = scalar_text(
        entry$medium_scenario, "medium_scenario"
      ),
      stringsAsFactors = FALSE
    )
  })
  rows <- do.call(rbind, rows)
  rownames(rows) <- NULL
  list(
    schema_version = "regcompass_directional_vmax_cache_contract_v1",
    mode = mode,
    solver = solver,
    flux_threshold = as.numeric(flux_threshold),
    rows = rows
  )
}

.rc_pack_microcompass_vmax_cache <- function(
    diagnostics, model_cache, mode, solver, flux_threshold) {
  tab <- as.data.frame(diagnostics)
  required <- c("row_id", "vmax", "feasible", "status")
  if (!all(required %in% colnames(tab))) {
    stop("Directional vmax diagnostics are incomplete.", call. = FALSE)
  }
  observed <- as.character(tab$row_id)
  expected <- names(model_cache)
  if (anyNA(observed) || any(!nzchar(observed)) || anyDuplicated(observed) ||
      !setequal(observed, expected)) {
    stop(
      "Directional vmax diagnostics do not match the structural model cache.",
      call. = FALSE
    )
  }
  tab <- tab[match(expected, observed), , drop = FALSE]
  values <- lapply(seq_along(expected), function(i) {
    vmax <- as.numeric(tab$vmax[[i]])
    feasible <- isTRUE(as.logical(tab$feasible[[i]]))
    status <- as.character(tab$status[[i]])
    if (length(vmax) != 1L ||
        (feasible && (!is.finite(vmax) || vmax < flux_threshold)) ||
        length(status) != 1L || is.na(status) || !nzchar(status)) {
      stop("Directional vmax diagnostics contain an invalid result.",
           call. = FALSE)
    }
    list(
      feasible = feasible,
      vmax = vmax,
      status = status,
      flux = numeric()
    )
  })
  names(values) <- expected
  attr(values, "vmax_cache_contract") <-
    .rc_microcompass_vmax_cache_contract(
      model_cache = model_cache,
      mode = mode,
      solver = solver,
      flux_threshold = flux_threshold
    )
  attr(values, "cache_source") <- "primary_directional_vmax"
  values
}

.rc_validate_microcompass_vmax_cache_override <- function(
    vmax_cache, model_cache, mode, solver, flux_threshold) {
  expected_names <- names(model_cache)
  if (!is.list(vmax_cache) || !length(vmax_cache) ||
      is.null(names(vmax_cache)) || anyNA(names(vmax_cache)) ||
      any(!nzchar(names(vmax_cache))) || anyDuplicated(names(vmax_cache)) ||
      !setequal(names(vmax_cache), expected_names)) {
    stop("The reused directional vmax cache is incomplete.", call. = FALSE)
  }
  expected_contract <- .rc_microcompass_vmax_cache_contract(
    model_cache = model_cache,
    mode = mode,
    solver = solver,
    flux_threshold = flux_threshold
  )
  observed_contract <- attr(
    vmax_cache, "vmax_cache_contract", exact = TRUE
  )
  if (!identical(observed_contract, expected_contract)) {
    stop(
      "The reused directional vmax cache contract does not match the current ",
      "structural models, solver, or flux threshold.",
      call. = FALSE
    )
  }
  answer <- vmax_cache[expected_names]
  for (row_id in expected_names) {
    value <- answer[[row_id]]
    if (!is.list(value) ||
        !all(c("feasible", "vmax", "status") %in% names(value))) {
      stop("The reused directional vmax cache contains a malformed result.",
           call. = FALSE)
    }
    vmax <- as.numeric(value$vmax)
    status <- as.character(value$status)
    if (length(vmax) != 1L ||
        (isTRUE(value$feasible) &&
         (!is.finite(vmax) || vmax < flux_threshold)) ||
        length(status) != 1L || is.na(status) || !nzchar(status)) {
      stop("The reused directional vmax cache contains an invalid result.",
           call. = FALSE)
    }
  }
  attr(answer, "vmax_cache_contract") <- expected_contract
  attr(answer, "cache_source") <- "primary_directional_vmax"
  attr(answer, "parallel_tasks") <- 0L
  attr(answer, "parallel_workers") <- 0L
  attr(answer, "parallel_scope") <-
    "reused_primary_directional_vmax_cache"
  answer
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
  override <- attr(
    model_cache, "directional_vmax_cache_override", exact = TRUE
  )
  if (!is.null(override)) {
    return(.rc_validate_microcompass_vmax_cache_override(
      vmax_cache = override,
      model_cache = model_cache,
      mode = mode,
      solver = solver,
      flux_threshold = flux_threshold
    ))
  }
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
        value <- rc_compass_vmax_directional(
          S = model$S,
          lb = model$lb,
          ub = model$ub,
          target_reaction = entry$reaction_id,
          direction = entry$target_direction,
          solver = solver,
          flux_threshold = flux_threshold
        )
        list(
          feasible = isTRUE(value$feasible),
          vmax = as.numeric(value$vmax),
          status = as.character(value$status),
          flux = numeric()
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
  attr(answer, "vmax_cache_contract") <-
    .rc_microcompass_vmax_cache_contract(
      model_cache = model_cache,
      mode = mode,
      solver = solver,
      flux_threshold = flux_threshold
    )
  attr(answer, "cache_source") <- "computed"
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
  n_metabolites <- nrow(S)
  n_reactions <- ncol(S)
  target_index <- match(target_reaction, reactions)
  reaction_index <- seq_len(n_reactions)
  s_nnz_per_col <- diff(S@p)
  s_j <- rep.int(reaction_index, s_nnz_per_col)
  s_i <- S@i + 1L

  A <- Matrix::sparseMatrix(
    i = c(
      s_i,
      n_metabolites + reaction_index,
      n_metabolites + reaction_index,
      n_metabolites + n_reactions + reaction_index,
      n_metabolites + n_reactions + reaction_index,
      n_metabolites + 2L * n_reactions + 1L
    ),
    j = c(
      s_j,
      reaction_index,
      n_reactions + reaction_index,
      reaction_index,
      n_reactions + reaction_index,
      target_index
    ),
    x = c(
      S@x,
      rep.int(1, n_reactions),
      rep.int(-1, n_reactions),
      rep.int(-1, n_reactions),
      rep.int(-1, n_reactions),
      if (identical(target_direction, "forward")) 1 else -1
    ),
    dims = c(
      n_metabolites + 2L * n_reactions + 1L,
      2L * n_reactions
    ),
    giveCsparse = TRUE
  )
  lhs <- c(
    rep(0, n_metabolites),
    rep(-Inf, 2L * n_reactions),
    omega * vmax
  )
  rhs <- c(
    rep(0, n_metabolites),
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
  progress_state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  args <- list(...)
  contexts <- .rc_layer2_model_contexts(
    args$model_cache,
    mode = as.character(args$mode %||% "meta_module_gem")
  )
  parts_dir <- .rc_layer2_cache_progress_dir(args$model_cache)
  run_kind <- progress_state$run_kind %||% "primary"
  reused <- !is.null(attr(
    args$model_cache, "directional_vmax_cache_override", exact = TRUE
  ))
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "directional_vmax_start", 1L, 4L,
      detail = paste0(
        "directional_targets=", item$n_targets,
        if (reused) "; source=primary_cache; solves=0" else "; source=lp"
      ),
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
      detail = paste0(
        "directional_targets=", item$n_targets,
        if (reused) "; reused_from=primary; solves=0" else "; solves=completed"
      ),
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
