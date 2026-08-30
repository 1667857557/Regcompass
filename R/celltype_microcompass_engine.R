# Directional LP scoring on cell-type-specific structural models.

.rc_celltype_model_contract <- function(model_cache) {
  if (!is.list(model_cache) || !length(model_cache)) {
    stop("The cell-type model cache is empty.", call. = FALSE)
  }
  files <- vapply(model_cache, function(entry) as.character(entry$file), character(1))
  representative <- match(unique(files), files)
  rows <- lapply(representative, function(i) {
    entry <- model_cache[[i]]
    model <- .rc_read_celltype_union_gem(
      entry$file, entry$cell_type, entry$medium_scenario,
      entry$file_checksum
    )
    validated <- rc_validate_gem(model)
    data.frame(
      cell_type = as.character(entry$cell_type),
      medium_scenario = as.character(entry$medium_scenario),
      model_file = as.character(entry$file),
      model_file_checksum = unname(tools::md5sum(entry$file)[[1L]]),
      n_reactions = ncol(validated$S),
      n_metabolites = nrow(validated$S),
      reaction_order_checksum =
        .rc_microcompass_object_checksum(colnames(validated$S)),
      metabolite_order_checksum =
        .rc_microcompass_object_checksum(rownames(validated$S)),
      stoichiometry_bounds_checksum =
        .rc_microcompass_object_checksum(list(
          S = validated$S, lb = validated$lb, ub = validated$ub
        )),
      shared_across_conditions = TRUE,
      shared_across_cell_types = FALSE,
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer[order(answer$cell_type, answer$medium_scenario), , drop = FALSE]
}

# Build the identical Step-2 LP constraint matrix directly in compressed sparse
# form. This replaces the previous zero/mass_balance/positive/negative/target
# intermediates, whose simultaneous residency multiplied peak RSS across workers.
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
    dims = c(n_metabolites + 2L * n_reactions + 1L,
             2L * n_reactions),
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

# Directional Vmax Step 1 needs only feasibility, scalar Vmax and status in the
# subsequent Step-2 objective. Drop the full primal flux vector inside each
# worker before BiocParallel returns it to the controller.
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
