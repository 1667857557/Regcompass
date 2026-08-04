# Reaction-granular Layer 2 parallel runtime and single-thread LP backends.

.rc_layer2_parallel_context <- new.env(parent = emptyenv())
.rc_layer2_parallel_context$active <- FALSE
.rc_layer2_parallel_context$parallel <- TRUE
.rc_layer2_parallel_context$BPPARAM <- NULL
.rc_layer2_parallel_context$nested_serial <- FALSE

.rc_regcompass_step_layer2_runtime_base <- rc_regcompass_step_layer2

rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  previous <- as.list(.rc_layer2_parallel_context)
  .rc_layer2_parallel_context$active <- TRUE
  .rc_layer2_parallel_context$parallel <- isTRUE(parallel)
  .rc_layer2_parallel_context$BPPARAM <- BPPARAM
  on.exit({
    rm(list = ls(.rc_layer2_parallel_context, all.names = TRUE),
       envir = .rc_layer2_parallel_context)
    list2env(previous, envir = .rc_layer2_parallel_context)
  }, add = TRUE)
  .rc_regcompass_step_layer2_runtime_base(
    layer1 = layer1,
    meta_modules = meta_modules,
    gem = gem,
    medium_scenarios = medium_scenarios,
    outdir = outdir,
    model_mode = model_mode,
    layer2_args = layer2_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
}

.rc_solve_lp_gurobi_base <- rc_solve_lp

rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = c("highs", "gurobi", "glpk"),
                        time_limit = Inf) {
  solver <- match.arg(solver)
  if (!identical(solver, "gurobi")) {
    return(.rc_solve_lp_gurobi_base(
      obj = obj, A = A, lhs = lhs, rhs = rhs, lb = lb, ub = ub,
      solver = solver, time_limit = time_limit
    ))
  }

  obj <- as.numeric(obj)
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  lhs <- as.numeric(lhs)
  rhs <- as.numeric(rhs)
  if (is.null(dim(A)) || length(dim(A)) != 2L) {
    stop("`A` must be a two-dimensional constraint matrix.", call. = FALSE)
  }
  if (ncol(A) != length(obj) || length(lb) != length(obj) ||
      length(ub) != length(obj)) {
    stop("LP objective, bounds and constraint columns are misaligned.",
         call. = FALSE)
  }
  if (nrow(A) != length(lhs) || nrow(A) != length(rhs)) {
    stop("LP constraint rows, lhs and rhs are misaligned.", call. = FALSE)
  }
  if (anyNA(obj) || any(!is.finite(obj)) || anyNA(lb) || anyNA(ub) ||
      anyNA(lhs) || anyNA(rhs) || any(lb > ub)) {
    stop("LP coefficients and bounds are invalid.", call. = FALSE)
  }
  if (!is.numeric(time_limit) || length(time_limit) != 1L ||
      is.na(time_limit) || time_limit <= 0) {
    stop("Internal solver `time_limit` must be one positive number or Inf.",
         call. = FALSE)
  }
  failure <- function(message) {
    list(
      status = "error", solution = numeric(), objective = NA_real_,
      solver = solver, solver_message = as.character(message)
    )
  }
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    return(failure("Package 'gurobi' is not installed."))
  }
  A <- .rc_as_dgCMatrix(A)
  expanded <- .rc_expand_ranged_constraints(A, lhs, rhs)
  model <- list(
    A = expanded$A,
    obj = obj,
    lb = lb,
    ub = ub,
    sense = ifelse(
      expanded$direction == "==", "=",
      ifelse(expanded$direction == ">=", ">", "<")
    ),
    rhs = expanded$rhs,
    modelsense = "min"
  )
  parameters <- list(
    OutputFlag = 0,
    Threads = 1
  )
  if (is.finite(time_limit)) parameters$TimeLimit <- time_limit
  answer <- tryCatch(
    gurobi::gurobi(model, params = parameters),
    error = function(e) e
  )
  if (inherits(answer, "error")) return(failure(conditionMessage(answer)))
  list(
    status = .rc_lp_status(answer$status %||% ""),
    solution = as.numeric(answer$x %||% numeric()),
    objective = as.numeric(answer$objval %||% NA_real_),
    solver = solver,
    solver_message = answer$status %||% ""
  )
}

.rc_atomic_save_rds <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(file), ".tmp_"),
    tmpdir = dirname(file)
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, file)) {
    if (!file.copy(temporary, file, overwrite = TRUE)) {
      stop("Cannot persist Layer 2 task result: ", file, call. = FALSE)
    }
    unlink(temporary, force = TRUE)
  }
  invisible(file)
}

.rc_layer2_task_bpparam <- function() {
  if (!isTRUE(.rc_layer2_parallel_context$active) ||
      !isTRUE(.rc_layer2_parallel_context$parallel) ||
      isTRUE(.rc_layer2_parallel_context$nested_serial)) {
    return(FALSE)
  }
  .rc_layer2_parallel_context$BPPARAM
}

.rc_directional_feasibility <- function(
    gem, targets, solver = "highs", time_limit = 60,
    flux_threshold = 1e-8) {
  required <- c("reaction_id", "target_direction")
  if (!is.data.frame(targets) || !all(required %in% colnames(targets))) {
    stop("`targets` must contain reaction_id and target_direction.",
         call. = FALSE)
  }
  if (!nrow(targets)) {
    return(data.frame(
      reaction_id = character(), target_direction = character(),
      feasible = logical(), vmax = numeric(), solver_status = character(),
      stringsAsFactors = FALSE
    ))
  }
  validated <- rc_validate_gem(gem)
  tasks <- split(targets, seq_len(nrow(targets)))
  run_one <- function(task) {
    on.exit(invisible(gc(verbose = FALSE, full = TRUE)), add = TRUE)
    reaction <- as.character(task$reaction_id[[1L]])
    direction <- as.character(task$target_direction[[1L]])
    if (!reaction %in% validated$reactions) {
      return(data.frame(
        reaction_id = reaction, target_direction = direction,
        feasible = FALSE, vmax = NA_real_, solver_status = "reaction_missing",
        stringsAsFactors = FALSE
      ))
    }
    if (!direction %in% c("forward", "reverse")) {
      return(data.frame(
        reaction_id = reaction, target_direction = direction,
        feasible = FALSE, vmax = 0,
        solver_status = "no_allowed_direction",
        stringsAsFactors = FALSE
      ))
    }
    answer <- rc_compass_vmax_directional(
      S = validated$S,
      lb = validated$lb,
      ub = validated$ub,
      target_reaction = reaction,
      direction = direction,
      solver = solver,
      time_limit = time_limit,
      flux_threshold = flux_threshold
    )
    data.frame(
      reaction_id = reaction,
      target_direction = direction,
      feasible = isTRUE(answer$feasible),
      vmax = answer$vmax,
      solver_status = answer$status,
      stringsAsFactors = FALSE
    )
  }
  rows <- rc_parallel_lapply(
    tasks, run_one, BPPARAM = .rc_layer2_task_bpparam()
  )
  do.call(rbind, rows)
}

.rc_build_celltype_medium_union_gem_cache <- function(
    gem, reaction_membership, core_reactions,
    target_reactions = NULL, medium_scenarios = NULL,
    celltype_col = "cell_type",
    cache_dir = tempfile("RegCompassR_celltype_medium_union_gem_"),
    target_direction = c("both", "forward", "reverse"),
    solver = "highs", time_limit = 300,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE) {
  target_direction <- match.arg(target_direction)
  scoped <- .rc_validate_celltype_reaction_scope(
    reaction_membership, core_reactions, celltype_col
  )
  reaction_membership <- scoped$reaction_membership
  core_reactions <- scoped$core_reactions
  if (!is.numeric(time_limit) || length(time_limit) != 1L ||
      !is.finite(time_limit) || time_limit <= 0) {
    stop("Cell-type union-GEM FASTCORE time limit must be positive.",
         call. = FALSE)
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  scenarios <- scenarios[!is.na(scenarios) & nzchar(scenarios)]
  if (!length(scenarios)) {
    stop("No medium scenarios are available for union-GEM construction.",
         call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  task_grid <- expand.grid(
    cell_type = scoped$cell_types,
    medium_scenario = scenarios,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  tasks <- split(task_grid, seq_len(nrow(task_grid)))

  run_one <- function(task, suppress_nested = FALSE) {
    previous_nested <- .rc_layer2_parallel_context$nested_serial
    if (isTRUE(suppress_nested)) {
      .rc_layer2_parallel_context$nested_serial <- TRUE
    }
    on.exit({
      .rc_layer2_parallel_context$nested_serial <- previous_nested
      invisible(gc(verbose = FALSE, full = TRUE))
    }, add = TRUE)
    cell_type <- as.character(task$cell_type[[1L]])
    scenario <- as.character(task$medium_scenario[[1L]])
    membership <- reaction_membership[
      reaction_membership[[celltype_col]] == cell_type, , drop = FALSE
    ]
    core <- core_reactions[
      core_reactions[[celltype_col]] == cell_type, , drop = FALSE
    ]
    selected_targets <- .rc_celltype_target_reactions(
      target_reactions, cell_type, celltype_col
    )
    if (!is.null(selected_targets)) {
      core <- core[core$reaction_id %in% selected_targets, , drop = FALSE]
    }
    if (!nrow(core)) {
      stop("No core targets remain for cell type `", cell_type, "`.",
           call. = FALSE)
    }
    medium <- medium_scenarios[
      as.character(medium_scenarios$medium_scenario_id) == scenario,
      , drop = FALSE
    ]
    if (!nrow(medium) ||
        (".no_constraints" %in% colnames(medium) &&
         all(medium$.no_constraints))) {
      medium <- NULL
    }
    model <- .rc_complete_celltype_medium_union_gem(
      gem = gem,
      reaction_membership = membership,
      core_reactions = core,
      cell_type = cell_type,
      medium_table = medium,
      target_direction = target_direction,
      solver = solver,
      time_limit = time_limit,
      fastcore_epsilon = fastcore_epsilon,
      max_support_reactions = max_support_reactions,
      strict = strict
    )
    model$shared_across_conditions <- TRUE
    model$shared_across_cell_types <- FALSE
    model$is_union_gem <- TRUE
    model$union_gem_cell_type <- cell_type
    model$union_gem_medium_scenario <- scenario
    model$union_gem_scope <-
      "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"

    file <- file.path(
      cache_dir,
      paste0(
        "union_gem__celltype_", .rc_safe_cache_token(cell_type),
        "__medium_", .rc_safe_cache_token(scenario), ".rds"
      )
    )
    .rc_atomic_save_rds(model, file)
    checksum <- unname(tools::md5sum(file))
    summary <- data.frame(
      cell_type = cell_type,
      medium_scenario = scenario,
      file = file,
      file_checksum = checksum,
      n_reactions = ncol(model$S),
      n_metabolites = nrow(model$S),
      n_celltype_biological_reactions =
        model$build_params$n_celltype_biological_reactions,
      n_celltype_fastcore_support_reactions =
        model$build_params$n_celltype_fastcore_support_reactions,
      target_status = model$target_status,
      build_strategy = "celltype_medium_union_gem",
      completion_stage =
        "parallel_celltype_specific_fastcore_after_condition_module_union",
      completion_time_limit = time_limit,
      stringsAsFactors = FALSE
    )

    cache <- list()
    if (nrow(model$target_directions)) {
      for (i in seq_len(nrow(model$target_directions))) {
        reaction <- as.character(model$target_directions$reaction_id[[i]])
        direction <- as.character(
          model$target_directions$target_direction[[i]]
        )
        key <- paste0(
          "celltype=", utils::URLencode(cell_type, reserved = TRUE),
          "::reaction=", utils::URLencode(reaction, reserved = TRUE),
          "::direction=", direction,
          "::medium=", utils::URLencode(scenario, reserved = TRUE)
        )
        cache[[key]] <- list(
          module_id = "CELLTYPE_MEDIUM_UNION_GEM",
          cell_type = cell_type,
          reaction_id = reaction,
          target_direction = direction,
          medium_scenario = scenario,
          condition = "all",
          file = file,
          file_checksum = checksum,
          build_strategy = "celltype_medium_union_gem"
        )
      }
    }
    rm(model)
    list(cache = cache, summary = summary)
  }

  parts <- if (length(tasks) > 1L) {
    rc_parallel_lapply(
      tasks,
      function(task) run_one(
        task[1, , drop = FALSE], suppress_nested = TRUE
      ),
      BPPARAM = .rc_layer2_task_bpparam()
    )
  } else {
    lapply(tasks, function(task) {
      run_one(task[1, , drop = FALSE], suppress_nested = FALSE)
    })
  }
  cache <- list()
  summaries <- vector("list", length(parts))
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    summaries[[i]] <- part$summary
    if (length(part$cache)) {
      duplicated <- intersect(names(cache), names(part$cache))
      if (length(duplicated)) {
        stop("Parallel FASTCORE tasks produced duplicate cache keys.",
             call. = FALSE)
      }
      cache[names(part$cache)] <- part$cache
    }
  }
  attr(cache, "summary") <- .rc_bind_frames_fill(summaries)
  attr(cache, "celltype_col") <- celltype_col
  attr(cache, "structural_scope") <- "cell_type_x_medium"
  attr(cache, "fastcore_parallel_task") <- "cell_type_x_medium"
  cache
}

.rc_build_microcompass_vmax_cache <- function(
    model_cache, mode, model_keys, solver, flux_threshold,
    parallel = TRUE, BPPARAM = NULL) {
  row_ids <- names(model_cache)
  if (is.null(row_ids) || anyNA(row_ids) || any(!nzchar(row_ids)) ||
      anyDuplicated(row_ids)) {
    stop("The directional model cache must have unique row IDs.",
         call. = FALSE)
  }
  if (!identical(names(model_keys), row_ids)) {
    model_keys <- model_keys[row_ids]
  }
  checkpoint_dir <- file.path(
    dirname(model_cache[[row_ids[[1L]]]]$file), "vmax_reaction_cache"
  )
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  tasks <- stats::setNames(as.list(row_ids), row_ids)

  run_one <- function(row_id) {
    on.exit(invisible(gc(verbose = FALSE, full = TRUE)), add = TRUE)
    entry <- model_cache[[row_id]]
    if (!identical(as.character(model_keys[[row_id]]),
                   as.character(entry$file))) {
      stop("Directional vmax task and model key are inconsistent.",
           call. = FALSE)
    }
    token <- substr(.rc_microcompass_object_checksum(list(
      row_id = row_id,
      file_checksum = entry$file_checksum,
      solver = solver,
      flux_threshold = flux_threshold
    )), 1L, 24L)
    checkpoint <- file.path(
      checkpoint_dir, paste0("vmax__", token, ".rds")
    )
    if (file.exists(checkpoint)) {
      cached <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
      if (is.list(cached) && identical(cached$row_id, row_id) &&
          identical(cached$file_checksum, entry$file_checksum) &&
          identical(cached$solver, solver) &&
          isTRUE(all.equal(cached$flux_threshold, flux_threshold)) &&
          is.list(cached$result)) {
        return(checkpoint)
      }
    }
    model <- .rc_load_microcompass_model(entry, mode)
    result <- rc_compass_vmax_directional(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_reaction = entry$reaction_id,
      direction = entry$target_direction,
      solver = solver,
      flux_threshold = flux_threshold
    )
    .rc_atomic_save_rds(list(
      row_id = row_id,
      file_checksum = entry$file_checksum,
      solver = solver,
      flux_threshold = flux_threshold,
      result = result
    ), checkpoint)
    rm(model, result)
    checkpoint
  }

  checkpoints <- rc_parallel_lapply(
    tasks, run_one,
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  answer <- lapply(checkpoints, function(file) readRDS(file)$result)
  names(answer) <- row_ids
  if (length(answer) != length(row_ids) ||
      !identical(names(answer), row_ids)) {
    stop("The shared directional vmax cache is incomplete.", call. = FALSE)
  }
  answer
}
