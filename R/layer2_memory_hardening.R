# Memory hardening for Layer 2 directional LP scoring.
#
# These helpers preserve the exact COMPASS mathematics and result semantics.
# They only remove avoidable controller/worker payloads and transient sparse
# matrix copies from the Step 5 execution path.

.rc_step5_compact_vmax_result <- function(value) {
  list(
    feasible = isTRUE(value$feasible),
    vmax = as.numeric(value$vmax),
    status = as.character(value$status),
    flux = numeric()
  )
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
  task_specs <- .rc_microcompass_vmax_tasks(model_keys, workers)

  tasks <- lapply(task_specs, function(task) {
    row_ids <- as.character(task$row_ids)
    entries <- model_cache[row_ids]
    list(
      model_key = as.character(task$model_key),
      row_ids = row_ids,
      entries = entries
    )
  })
  names(tasks) <- names(task_specs)

  grouped <- rc_parallel_lapply(
    tasks,
    function(task) {
      row_ids <- as.character(task$row_ids)
      entries <- task$entries
      first_entry <- entries[[1L]]
      model <- .rc_load_microcompass_model(first_entry, mode)
      on.exit({
        rm(model)
        invisible(gc(verbose = FALSE, full = FALSE))
      }, add = TRUE)

      values <- lapply(row_ids, function(row_id) {
        entry <- entries[[row_id]]
        value <- rc_compass_vmax_directional(
          S = model$S,
          lb = model$lb,
          ub = model$ub,
          target_reaction = entry$reaction_id,
          direction = entry$target_direction,
          solver = solver,
          flux_threshold = flux_threshold
        )
        .rc_step5_compact_vmax_result(value)
      })
      names(values) <- row_ids
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
  attr(answer, "memory_contract") <-
    "worker_compacts_directional_flux_before_controller_return"
  answer
}

.rc_step5_abs_flux_constraint_matrix <- function(S, target_index, target_sign) {
  S <- .rc_as_dgCMatrix(S)
  n_metabolites <- nrow(S)
  n_reactions <- ncol(S)
  nnz <- length(S@x)
  mass_cols <- if (nnz) {
    rep.int(seq_len(n_reactions), diff(S@p))
  } else {
    integer()
  }
  reaction_index <- seq_len(n_reactions)
  positive_rows <- n_metabolites + reaction_index
  negative_rows <- n_metabolites + n_reactions + reaction_index
  target_row <- n_metabolites + 2L * n_reactions + 1L

  Matrix::sparseMatrix(
    i = c(
      S@i + 1L,
      positive_rows,
      positive_rows,
      negative_rows,
      negative_rows,
      target_row
    ),
    j = c(
      mass_cols,
      reaction_index,
      n_reactions + reaction_index,
      reaction_index,
      n_reactions + reaction_index,
      target_index
    ),
    x = c(
      S@x,
      rep(1, n_reactions),
      rep(-1, n_reactions),
      rep(-1, n_reactions),
      rep(-1, n_reactions),
      target_sign
    ),
    dims = c(
      n_metabolites + 2L * n_reactions + 1L,
      2L * n_reactions
    ),
    repr = "C"
  )
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
  target_index <- match(target_reaction, reactions)
  target_sign <- if (identical(target_direction, "forward")) 1 else -1
  A <- .rc_step5_abs_flux_constraint_matrix(
    S = S,
    target_index = target_index,
    target_sign = target_sign
  )
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
