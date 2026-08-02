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

.rc_build_microcompass_vmax_cache <- function(
    model_cache, mode, model_keys, solver, flux_threshold,
    parallel = TRUE, BPPARAM = NULL) {
  unique_model_keys <- unique(unname(model_keys))
  tasks <- stats::setNames(as.list(unique_model_keys), unique_model_keys)
  grouped <- rc_parallel_lapply(
    tasks,
    function(model_key) {
      selected_rows <- names(model_keys)[model_keys == model_key]
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
  .rc_flatten_microcompass_vmax_cache(grouped, names(model_cache))
}

.rc_compass_step2_from_vmax_directional <- function(
    S, lb, ub, target_reaction, penalties, vmax_result,
    target_direction = c("forward", "reverse"),
    omega = 0.95,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8) {
  target_direction <- match.arg(target_direction)
  solver <- match.arg(solver)
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
  if (!isTRUE(vmax_result$feasible)) {
    return(list(
      feasible = FALSE,
      penalty = NA_real_,
      vmax = as.numeric(vmax_result$vmax),
      solver_status = as.character(vmax_result$status),
      step1_status = as.character(vmax_result$status),
      step2_status = "not_run",
      flux = numeric()
    ))
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
  step2 <- rc_solve_lp(
    obj = c(rep(0, n_reactions), penalties),
    A = A,
    lhs = lhs,
    rhs = rhs,
    lb = c(lb, rep(0, n_reactions)),
    ub = c(ub, auxiliary_upper),
    solver = solver
  )
  if (!identical(step2$status, "optimal") ||
      length(step2$solution) != 2L * n_reactions) {
    return(list(
      feasible = FALSE,
      penalty = NA_real_,
      vmax = vmax,
      solver_status = step2$status,
      step1_status = as.character(vmax_result$status),
      step2_status = step2$status,
      flux = numeric()
    ))
  }
  flux <- step2$solution[seq_len(n_reactions)]
  names(flux) <- reactions
  list(
    feasible = TRUE,
    penalty = max(0, as.numeric(step2$objective)),
    vmax = vmax,
    solver_status = step2$status,
    step1_status = as.character(vmax_result$status),
    step2_status = step2$status,
    flux = flux
  )
}
