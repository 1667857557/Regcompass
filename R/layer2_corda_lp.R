# Original CORDA directional dependency assessment and accelerated LP runtime.

.rc_corda_split_model <- function(gem, tolerance = 1e-8) {
  validated <- rc_validate_gem(gem)
  reactions <- validated$reactions
  columns <- list()
  rows <- list()
  lower <- upper <- numeric()
  variable_id <- character()
  index <- 0L
  for (i in seq_along(reactions)) {
    reaction <- reactions[[i]]
    lb <- validated$lb[[i]]
    ub <- validated$ub[[i]]
    if (ub > tolerance) {
      index <- index + 1L
      variable <- paste0(reaction, "::forward")
      variable_id[[index]] <- variable
      columns[[index]] <- validated$S[, i, drop = FALSE]
      lower[[index]] <- max(0, lb)
      upper[[index]] <- ub
      rows[[index]] <- data.frame(
        variable_id = variable,
        reaction_id = reaction,
        direction = "forward",
        original_index = i,
        stringsAsFactors = FALSE
      )
    }
    if (lb < -tolerance) {
      index <- index + 1L
      variable <- paste0(reaction, "::reverse")
      variable_id[[index]] <- variable
      columns[[index]] <- -validated$S[, i, drop = FALSE]
      lower[[index]] <- max(0, -ub)
      upper[[index]] <- -lb
      rows[[index]] <- data.frame(
        variable_id = variable,
        reaction_id = reaction,
        direction = "reverse",
        original_index = i,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(columns)) {
    stop("CORDA parent model has no allowed directional variables.",
         call. = FALSE)
  }
  S <- do.call(cbind, columns)
  colnames(S) <- variable_id
  lower <- stats::setNames(as.numeric(lower), variable_id)
  upper <- stats::setNames(as.numeric(upper), variable_id)
  direction_table <- do.call(rbind, rows)
  rownames(direction_table) <- NULL
  list(
    S = .rc_as_dgCMatrix(S),
    lb = lower,
    ub = upper,
    direction_table = direction_table,
    variable_to_reaction = stats::setNames(
      direction_table$reaction_id, direction_table$variable_id
    ),
    variable_to_direction = stats::setNames(
      direction_table$direction, direction_table$variable_id
    ),
    algorithm = "single_direction_split_before_CORDA",
    tolerance = tolerance
  )
}

.rc_corda_block_reactions <- function(split, reactions) {
  reactions <- unique(as.character(reactions))
  blocked <- split$direction_table$variable_id[
    split$direction_table$reaction_id %in% reactions
  ]
  split$ub[blocked] <- 0
  split$lb[blocked] <- pmin(split$lb[blocked], split$ub[blocked])
  split
}

.rc_corda_target_table <- function(split, reactions) {
  reactions <- unique(as.character(reactions))
  answer <- split$direction_table[
    split$direction_table$reaction_id %in% reactions &
      split$ub[split$direction_table$variable_id] > split$tolerance,
    c("variable_id", "reaction_id", "direction"),
    drop = FALSE
  ]
  rownames(answer) <- NULL
  answer
}

.rc_corda_string_hash <- function(value) {
  raw <- as.integer(charToRaw(enc2utf8(paste(value, collapse = "|"))))
  hash <- 0
  if (length(raw)) {
    for (x in raw) hash <- (hash * 131 + x) %% 2147483629
  }
  as.integer(hash)
}

.rc_corda_noise <- function(n, seed, key, kappa) {
  if (n <= 0L || kappa <= 0) return(rep(0, n))
  derived <- (as.double(seed) + .rc_corda_string_hash(key)) %% 2147483646
  derived <- as.integer(derived + 1)
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (existed) {
      assign(".Random.seed", old, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(derived)
  stats::runif(n, min = 0, max = kappa)
}

.rc_corda_base_cost <- function(split, confidence, stage, gamma) {
  reaction_confidence <- as.character(
    confidence[split$direction_table$reaction_id]
  )
  cost <- rep(0, nrow(split$direction_table))
  if (identical(stage, "stage1_hc_dependencies")) {
    cost[grepl("^MC", reaction_confidence)] <- 1
    cost[reaction_confidence == "NC"] <- gamma
  } else if (identical(stage, "stage2_mc_nc_support")) {
    cost[reaction_confidence == "NC"] <- gamma
  } else if (identical(stage, "stage3_re_ot_dependencies")) {
    cost[reaction_confidence == "OT"] <- 1
  } else {
    stop("Unknown CORDA dependency stage: ", stage, call. = FALSE)
  }
  stats::setNames(cost, split$direction_table$variable_id)
}

.rc_corda_highs_api_available <- function() {
  if (!requireNamespace("highs", quietly = TRUE)) return(FALSE)
  required <- c(
    "highs_model", "hi_new_solver", "hi_solver_set_objective",
    "hi_solver_set_variable_bounds", "hi_solver_run",
    "hi_solver_status_message", "hi_solver_get_solution",
    "hi_solver_info", "hi_solver_set_option"
  )
  all(required %in% getNamespaceExports("highs"))
}

.rc_corda_highs_call <- function(name, ...) {
  do.call(getExportedValue("highs", name), list(...))
}

.rc_corda_new_lp_engine <- function(split, solver, time_limit) {
  if (identical(solver, "highs") && .rc_corda_highs_api_available()) {
    persistent <- tryCatch({
      model <- .rc_corda_highs_call(
        "highs_model",
        L = rep(0, ncol(split$S)),
        lower = as.numeric(split$lb),
        upper = as.numeric(split$ub),
        A = split$S,
        lhs = rep(0, nrow(split$S)),
        rhs = rep(0, nrow(split$S)),
        maximum = FALSE
      )
      pointer <- .rc_corda_highs_call("hi_new_solver", model)
      try(.rc_corda_highs_call(
        "hi_solver_set_option", pointer, "output_flag", "off"
      ), silent = TRUE)
      try(.rc_corda_highs_call(
        "hi_solver_set_option", pointer, "threads", "1"
      ), silent = TRUE)
      try(.rc_corda_highs_call(
        "hi_solver_set_option", pointer, "solver", "simplex"
      ), silent = TRUE)
      if (is.finite(time_limit)) {
        try(.rc_corda_highs_call(
          "hi_solver_set_option", pointer, "time_limit",
          as.character(time_limit)
        ), silent = TRUE)
      }
      list(
        type = "highs_persistent_cpp",
        pointer = pointer,
        split = split,
        solver = solver,
        time_limit = time_limit,
        n_solves = 0L,
        n_fallback = 0L
      )
    }, error = function(e) NULL)
    if (!is.null(persistent)) return(persistent)
  }
  list(
    type = "one_shot",
    pointer = NULL,
    split = split,
    solver = solver,
    time_limit = time_limit,
    n_solves = 0L,
    n_fallback = 0L
  )
}

.rc_corda_one_shot_solve <- function(
    engine, objective, lower, upper) {
  answer <- rc_solve_lp(
    obj = objective,
    A = engine$split$S,
    lhs = rep(0, nrow(engine$split$S)),
    rhs = rep(0, nrow(engine$split$S)),
    lb = lower,
    ub = upper,
    solver = engine$solver,
    time_limit = engine$time_limit
  )
  list(
    status = answer$status,
    solution = answer$solution,
    objective = answer$objective,
    backend = "one_shot",
    solver_message = answer$solver_message %||% ""
  )
}

.rc_corda_engine_solve <- function(
    engine, objective, lower, upper) {
  engine$n_solves <- engine$n_solves + 1L
  if (!identical(engine$type, "highs_persistent_cpp")) {
    return(list(
      answer = .rc_corda_one_shot_solve(engine, objective, lower, upper),
      engine = engine
    ))
  }
  answer <- tryCatch({
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
  if (inherits(answer, "error")) {
    engine$n_fallback <- engine$n_fallback + 1L
    answer <- .rc_corda_one_shot_solve(engine, objective, lower, upper)
    answer$backend <- "highs_persistent_failed_one_shot_fallback"
    answer$solver_message <- paste(
      conditionMessage(answer), answer$solver_message
    )
  }
  list(answer = answer, engine = engine)
}

.rc_corda_dependency_task <- function(
    engine, task, confidence, options) {
  split <- engine$split
  target <- as.character(task$variable_id[[1L]])
  target_index <- match(target, colnames(split$S))
  stage <- as.character(task$stage[[1L]])
  replicate <- as.integer(task$replicate[[1L]])
  base_cost <- .rc_corda_base_cost(
    split, confidence, stage, options$gamma
  )
  objective <- as.numeric(base_cost) + .rc_corda_noise(
    length(base_cost), options$seed,
    c(stage, target, replicate), options$kappa
  )
  lower <- split$lb
  upper <- split$ub
  if (is.na(target_index) || upper[[target_index]] < options$epsilon) {
    return(list(
      result = list(
        task = task, status = "target_blocked", associated = character(),
        target_flux = 0, objective = NA_real_, backend = engine$type
      ),
      engine = engine
    ))
  }
  lower[[target_index]] <- max(lower[[target_index]], options$epsilon)
  solved <- .rc_corda_engine_solve(engine, objective, lower, upper)
  engine <- solved$engine
  answer <- solved$answer
  flux <- answer$solution
  associated <- character()
  target_flux <- NA_real_
  if (identical(answer$status, "optimal") &&
      length(flux) == ncol(split$S)) {
    names(flux) <- colnames(split$S)
    target_flux <- flux[[target]]
    penalized <- names(base_cost)[base_cost > 0]
    active <- penalized[flux[penalized] > options$flux_tolerance]
    associated <- unique(as.character(
      split$variable_to_reaction[active]
    ))
    associated <- setdiff(
      associated, as.character(task$reaction_id[[1L]])
    )
  }
  list(
    result = list(
      task = task,
      status = answer$status,
      associated = sort(associated),
      target_flux = target_flux,
      objective = answer$objective,
      backend = answer$backend,
      solver_message = answer$solver_message
    ),
    engine = engine
  )
}

.rc_corda_feasibility_task <- function(engine, task, options) {
  split <- engine$split
  target <- as.character(task$variable_id[[1L]])
  target_index <- match(target, colnames(split$S))
  objective <- rep(0, ncol(split$S))
  lower <- split$lb
  upper <- split$ub
  if (is.na(target_index) || upper[[target_index]] < options$epsilon) {
    return(list(
      result = list(
        task = task, status = "target_blocked", associated = character(),
        target_flux = 0, objective = NA_real_, backend = engine$type
      ),
      engine = engine
    ))
  }
  objective[[target_index]] <- -1
  solved <- .rc_corda_engine_solve(engine, objective, lower, upper)
  engine <- solved$engine
  answer <- solved$answer
  target_flux <- NA_real_
  status <- answer$status
  if (identical(status, "optimal") &&
      length(answer$solution) == ncol(split$S)) {
    target_flux <- as.numeric(answer$solution[[target_index]])
    if (!is.finite(target_flux) || target_flux < options$epsilon) {
      status <- "blocked"
    }
  }
  list(
    result = list(
      task = task,
      status = status,
      associated = character(),
      target_flux = target_flux,
      objective = answer$objective,
      backend = answer$backend,
      solver_message = answer$solver_message
    ),
    engine = engine
  )
}

.rc_corda_task_bpparam <- function() {
  .rc_layer2_task_bpparam()
}

.rc_corda_worker_count <- function(BPPARAM, n_tasks) {
  if (identical(BPPARAM, FALSE) || n_tasks < 2L) return(1L)
  if (!is.null(BPPARAM) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(BPPARAM, "BiocParallelParam")) {
    return(max(1L, min(n_tasks, BiocParallel::bpnworkers(BPPARAM))))
  }
  config <- rc_parallel_config(workers = NULL, backend = "auto")
  max(1L, min(n_tasks, config$workers))
}

.rc_corda_run_tasks <- function(
    split, tasks, confidence, options, solver, time_limit,
    BPPARAM = .rc_corda_task_bpparam()) {
  if (!is.data.frame(tasks) || !nrow(tasks)) {
    return(list(
      results = list(),
      execution = list(
        n_tasks = 0L, n_chunks = 0L, workers = 1L,
        task_granularity = "direction_x_replicate",
        solver_runtime = "not_run"
      )
    ))
  }
  workers <- .rc_corda_worker_count(BPPARAM, nrow(tasks))
  n_chunks <- min(workers, nrow(tasks))
  chunk_id <- rep(seq_len(n_chunks), length.out = nrow(tasks))
  chunks <- split(tasks, chunk_id)
  run_chunk <- function(chunk) {
    engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
    rows <- vector("list", nrow(chunk))
    for (i in seq_len(nrow(chunk))) {
      task <- chunk[i, , drop = FALSE]
      solved <- if (identical(as.character(task$kind[[1L]]), "dependency")) {
        .rc_corda_dependency_task(engine, task, confidence, options)
      } else {
        .rc_corda_feasibility_task(engine, task, options)
      }
      engine <- solved$engine
      rows[[i]] <- solved$result
    }
    list(
      results = rows,
      engine = list(
        type = engine$type,
        n_solves = engine$n_solves,
        n_fallback = engine$n_fallback
      )
    )
  }
  parts <- rc_parallel_lapply(
    chunks, run_chunk,
    BPPARAM = if (n_chunks > 1L) BPPARAM else FALSE
  )
  results <- unlist(lapply(parts, `[[`, "results"), recursive = FALSE)
  engine_rows <- lapply(parts, function(part) {
    data.frame(
      type = part$engine$type,
      n_solves = part$engine$n_solves,
      n_fallback = part$engine$n_fallback,
      stringsAsFactors = FALSE
    )
  })
  engines <- do.call(rbind, engine_rows)
  list(
    results = results,
    execution = list(
      n_tasks = nrow(tasks),
      n_chunks = n_chunks,
      workers = workers,
      task_granularity = "direction_x_replicate",
      stage_barrier = TRUE,
      persistent_solver = any(
        engines$type == "highs_persistent_cpp"
      ),
      solver_runtime = paste(unique(engines$type), collapse = ";"),
      n_solves = sum(engines$n_solves),
      n_fallback = sum(engines$n_fallback)
    )
  )
}

.rc_corda_make_tasks <- function(
    split, reactions, stage, n = 1L,
    kind = c("dependency", "feasibility")) {
  kind <- match.arg(kind)
  targets <- .rc_corda_target_table(split, reactions)
  if (!nrow(targets)) return(data.frame())
  if (identical(kind, "dependency")) {
    targets <- targets[rep(seq_len(nrow(targets)), each = n), , drop = FALSE]
    targets$replicate <- rep(seq_len(n), times = nrow(targets) / n)
  } else {
    targets$replicate <- 1L
  }
  targets$stage <- stage
  targets$kind <- kind
  rownames(targets) <- NULL
  targets
}

.rc_corda_results_table <- function(results) {
  if (!length(results)) return(data.frame())
  rows <- lapply(results, function(x) {
    task <- x$task
    data.frame(
      variable_id = as.character(task$variable_id[[1L]]),
      reaction_id = as.character(task$reaction_id[[1L]]),
      direction = as.character(task$direction[[1L]]),
      stage = as.character(task$stage[[1L]]),
      replicate = as.integer(task$replicate[[1L]]),
      kind = as.character(task$kind[[1L]]),
      status = as.character(x$status),
      target_flux = as.numeric(x$target_flux),
      objective = as.numeric(x$objective),
      backend = as.character(x$backend),
      n_associated = length(x$associated),
      associated = paste(x$associated, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.rc_corda_associated <- function(results, allowed = NULL) {
  value <- unique(unlist(lapply(results, `[[`, "associated"), use.names = FALSE))
  value <- value[!is.na(value) & nzchar(value)]
  if (!is.null(allowed)) value <- intersect(value, allowed)
  sort(value)
}
