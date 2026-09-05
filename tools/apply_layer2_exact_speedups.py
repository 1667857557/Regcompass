from pathlib import Path
import re


def replace_function(path, name, new_text):
    p = Path(path)
    text = p.read_text()
    m = re.search(rf'(?m)^{re.escape(name)}\s*<-\s*function\b', text)
    if not m:
        raise RuntimeError(f'Cannot find {name} in {path}')
    brace = text.find('{', m.end())
    if brace < 0:
        raise RuntimeError(f'Cannot find opening brace for {name}')
    depth = 0
    quote = None
    escape = False
    comment = False
    end = None
    for i in range(brace, len(text)):
        ch = text[i]
        if comment:
            if ch == '\n':
                comment = False
            continue
        if quote is not None:
            if escape:
                escape = False
            elif ch == '\\':
                escape = True
            elif ch == quote:
                quote = None
            continue
        if ch in ('"', "'", '`'):
            quote = ch
            continue
        if ch == '#':
            comment = True
            continue
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    if end is None:
        raise RuntimeError(f'Cannot find closing brace for {name}')
    p.write_text(text[:m.start()] + new_text.rstrip() + text[end:])


def insert_before(path, marker, block):
    p = Path(path)
    text = p.read_text()
    if block.strip() in text:
        return
    idx = text.find(marker)
    if idx < 0:
        raise RuntimeError(f'Cannot find marker {marker!r} in {path}')
    p.write_text(text[:idx] + block.rstrip() + '\n\n' + text[idx:])


def replace_once(path, old, new):
    p = Path(path)
    text = p.read_text()
    if text.count(old) != 1:
        raise RuntimeError(
            f'Expected exactly one occurrence in {path}, found {text.count(old)}: {old[:80]!r}'
        )
    p.write_text(text.replace(old, new, 1))


# ---- Solver/native hot path -------------------------------------------------
solver_path = 'R/microcompass_vmax_cache.R'
replace_function(solver_path, '.rc_microcompass_highs_api_available', r'''.rc_microcompass_highs_api_available <- function() {
  if (!requireNamespace("highs", quietly = TRUE)) return(FALSE)
  required <- c(
    "highs_model", "hi_new_solver", "hi_solver_set_objective",
    "hi_solver_run", "hi_solver_status_message",
    "hi_solver_get_solution", "hi_solver_info", "hi_solver_set_option"
  )
  all(required %in% getNamespaceExports("highs"))
}''')

replace_function(solver_path, '.rc_microcompass_highs_call', r'''.rc_microcompass_highs_call <- function(name, ...) {
  cache <- get0(
    ".rc_microcompass_highs_api_cache",
    envir = asNamespace("RegCompassR"), inherits = FALSE
  )
  if (!is.environment(cache)) {
    cache <- new.env(parent = emptyenv())
    assign(
      ".rc_microcompass_highs_api_cache", cache,
      envir = asNamespace("RegCompassR")
    )
  }
  fun <- get0(name, envir = cache, mode = "function", inherits = FALSE)
  if (!is.function(fun)) {
    fun <- getExportedValue("highs", name)
    assign(name, fun, envir = cache)
  }
  fun(...)
}''')

insert_before(
    solver_path,
    '.rc_microcompass_highs_api_available <- function()',
    '.rc_microcompass_highs_api_cache <- new.env(parent = emptyenv())'
)

insert_before(
    solver_path,
    '.rc_compass_step2_prepare <- function(',
    r'''.rc_step2_penalty_evidence_stats <- function(penalty) {
  if (!is.numeric(penalty) || is.null(dim(penalty)) || nrow(penalty) < 1L) {
    stop("Step 2 penalty evidence requires a non-empty numeric matrix.",
         call. = FALSE)
  }
  unavailable <- vapply(
    seq_len(ncol(penalty)),
    function(i) sum(!is.finite(penalty[, i])),
    integer(1)
  )
  list(
    unavailable = unavailable,
    fraction = 1 - unavailable / nrow(penalty),
    all_finite = unavailable == 0L
  )
}'''
)

replace_function(solver_path, '.rc_compass_step2_new_engine', r'''.rc_compass_step2_new_engine <- function(
    template, solver, persistent_required = FALSE) {
  solver <- match.arg(solver, c("highs", "gurobi", "glpk"))
  n_reactions <- as.integer(template$n_reactions)
  engine <- list(
    type = "one_shot",
    pointer = NULL,
    template = template,
    solver = solver,
    current_penalties = rep(0, n_reactions),
    objective_index = as.integer(n_reactions + seq_len(n_reactions) - 1L),
    n_solves = 0L,
    n_objective_updates = 0L,
    n_fallback = 0L,
    persistent_disabled = FALSE,
    persistent_required = isTRUE(persistent_required)
  )
  if (!identical(solver, "highs")) return(engine)
  if (!.rc_microcompass_highs_api_available()) {
    if (isTRUE(persistent_required)) {
      stop(
        "Target-parallel Layer 2 requires the persistent HiGHS solver API. ",
        "Install a current `highs` package (>= 1.12.0-1).",
        call. = FALSE
      )
    }
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
    if (isTRUE(persistent_required)) {
      stop(
        "Persistent HiGHS initialization failed for target-parallel Layer 2: ",
        conditionMessage(persistent),
        call. = FALSE
      )
    }
    engine$persistent_disabled <- TRUE
    engine$persistent_message <- conditionMessage(persistent)
    return(engine)
  }
  engine$type <- "highs_persistent_cpp"
  engine$pointer <- persistent
  engine
}''')

replace_function(solver_path, '.rc_compass_step2_engine_solve', r'''.rc_compass_step2_engine_solve <- function(
    engine, penalties, return_solution = TRUE, trusted_aligned = FALSE) {
  if (isTRUE(trusted_aligned)) {
    penalties <- as.numeric(penalties)
    if (length(penalties) != engine$template$n_reactions ||
        any(!is.finite(penalties)) || any(penalties < 0)) {
      stop("Trusted Step 2 penalties are not finite non-negative aligned values.",
           call. = FALSE)
    }
  } else {
    penalties <- .rc_compass_step2_align_penalties(
      engine$template$reactions, penalties
    )
  }
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
        index = engine$objective_index[changed],
        coeff = penalties[changed]
      )
    }
    .rc_microcompass_highs_call("hi_solver_run", engine$pointer)
    status_message <- .rc_microcompass_highs_call(
      "hi_solver_status_message", engine$pointer
    )
    status <- .rc_lp_status(status_message)
    info <- tryCatch(
      .rc_microcompass_highs_call("hi_solver_info", engine$pointer),
      error = function(e) list()
    )
    objective <- as.numeric(
      info$objective_function_value %||% NA_real_
    )
    solution <- numeric()
    if (identical(status, "optimal") &&
        (isTRUE(return_solution) || !is.finite(objective))) {
      solution <- as.numeric(.rc_microcompass_highs_call(
        "hi_solver_get_solution", engine$pointer
      )$col_value)
    }
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
      backend = if (isTRUE(return_solution)) {
        "highs_persistent_cpp_basis_reuse"
      } else {
        "highs_persistent_cpp_basis_reuse_objective_only"
      },
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
  if (isTRUE(engine$persistent_required)) {
    .rc_compass_step2_release_engine(engine)
    stop(
      "Persistent HiGHS failed inside target-parallel Layer 2 instead of ",
      "silently rebuilding millions of one-shot LPs: ",
      engine$persistent_message,
      call. = FALSE
    )
  }
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
}''')

replace_function(solver_path, '.rc_compass_step2_result', r'''.rc_compass_step2_result <- function(
    template, answer, require_solution = TRUE) {
  n_reactions <- template$n_reactions
  optimal <- identical(answer$status, "optimal")
  objective <- as.numeric(answer$objective %||% NA_real_)
  solution_ok <- length(answer$solution) == 2L * n_reactions
  if (!optimal || !is.finite(objective) ||
      (isTRUE(require_solution) && !solution_ok)) {
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
  flux <- if (isTRUE(require_solution)) {
    value <- answer$solution[seq_len(n_reactions)]
    names(value) <- template$reactions
    value
  } else {
    numeric()
  }
  list(
    feasible = TRUE,
    penalty = max(0, objective),
    vmax = template$vmax,
    solver_status = answer$status,
    step1_status = template$step1_status,
    step2_status = answer$status,
    solver_backend = answer$backend %||% "unknown",
    flux = flux
  )
}''')

insert_before(
    solver_path,
    '.rc_build_microcompass_vmax_cache_core <- function(',
    r'''.rc_compass_vmax_batch_highs <- function(
    S, lb, ub, target_reaction, direction, flux_threshold = 1e-8) {
  if (!.rc_microcompass_highs_api_available()) return(NULL)
  S <- .rc_as_dgCMatrix(S)
  reactions <- colnames(S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions)) ||
      anyDuplicated(reactions)) {
    stop("`S` must have unique non-empty reaction IDs in colnames().",
         call. = FALSE)
  }
  target_reaction <- as.character(target_reaction)
  direction <- as.character(direction)
  if (!length(target_reaction) || length(direction) != length(target_reaction) ||
      any(!direction %in% c("forward", "reverse")) ||
      any(!target_reaction %in% reactions)) {
    stop("Persistent Vmax targets or directions are invalid.", call. = FALSE)
  }
  lb <- rc_align_bound(lb, reactions, default = -1000, name = "lb")
  ub <- rc_align_bound(ub, reactions, default = 1000, name = "ub")
  if (any(lb > ub)) {
    stop("Reaction lower bounds cannot exceed upper bounds.", call. = FALSE)
  }

  pointer <- NULL
  on.exit({
    if (!is.null(pointer)) {
      .rc_compass_step2_release_engine(list(pointer = pointer))
    }
  }, add = TRUE)
  answer <- tryCatch({
    model <- .rc_microcompass_highs_call(
      "highs_model",
      L = rep(0, length(reactions)),
      lower = lb,
      upper = ub,
      A = S,
      lhs = rep(0, nrow(S)),
      rhs = rep(0, nrow(S)),
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

    values <- vector("list", length(target_reaction))
    previous_index <- NA_integer_
    for (i in seq_along(target_reaction)) {
      target_index <- match(target_reaction[[i]], reactions)
      allowed <- if (identical(direction[[i]], "forward")) {
        ub[[target_index]] > flux_threshold
      } else {
        lb[[target_index]] < -flux_threshold
      }
      if (!isTRUE(allowed)) {
        values[[i]] <- list(
          feasible = FALSE, vmax = 0, status = "no_allowed_direction",
          flux = numeric()
        )
        next
      }

      update_index <- unique(c(previous_index, target_index))
      update_index <- update_index[is.finite(update_index)]
      coeff <- numeric(length(update_index))
      coeff[match(target_index, update_index)] <-
        if (identical(direction[[i]], "forward")) -1 else 1
      .rc_microcompass_highs_call(
        "hi_solver_set_objective", pointer,
        index = as.integer(update_index - 1L), coeff = coeff
      )
      previous_index <- target_index
      .rc_microcompass_highs_call("hi_solver_run", pointer)
      status_message <- .rc_microcompass_highs_call(
        "hi_solver_status_message", pointer
      )
      status <- .rc_lp_status(status_message)
      if (!identical(status, "optimal")) {
        values[[i]] <- list(
          feasible = FALSE, vmax = NA_real_, status = status,
          flux = numeric()
        )
        next
      }
      info <- tryCatch(
        .rc_microcompass_highs_call("hi_solver_info", pointer),
        error = function(e) list()
      )
      objective <- as.numeric(info$objective_function_value %||% NA_real_)
      vmax <- if (is.finite(objective)) {
        max(0, -objective)
      } else {
        solution <- as.numeric(.rc_microcompass_highs_call(
          "hi_solver_get_solution", pointer
        )$col_value)
        if (length(solution) != length(reactions)) NA_real_ else if (
          identical(direction[[i]], "forward")
        ) {
          solution[[target_index]]
        } else {
          -solution[[target_index]]
        }
      }
      feasible <- is.finite(vmax) && vmax >= flux_threshold
      values[[i]] <- list(
        feasible = feasible,
        vmax = if (is.finite(vmax)) max(0, vmax) else NA_real_,
        status = if (feasible) "optimal" else "blocked",
        flux = numeric()
      )
    }
    values
  }, error = function(e) e)
  if (inherits(answer, "error")) return(NULL)
  answer
}'''
)

old_vmax_callback = r'''      model <- .rc_load_microcompass_model(first_entry, mode)
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
      names(values) <- selected_rows'''
new_vmax_callback = r'''      model <- .rc_load_microcompass_model(first_entry, mode)
      entries <- lapply(selected_rows, function(row_id) model_cache[[row_id]])
      values <- if (identical(solver, "highs")) {
        .rc_compass_vmax_batch_highs(
          S = model$S,
          lb = model$lb,
          ub = model$ub,
          target_reaction = vapply(
            entries, function(entry) as.character(entry$reaction_id),
            character(1)
          ),
          direction = vapply(
            entries, function(entry) as.character(entry$target_direction),
            character(1)
          ),
          flux_threshold = flux_threshold
        )
      } else {
        NULL
      }
      if (is.null(values)) {
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
      }
      names(values) <- selected_rows'''
replace_once(solver_path, old_vmax_callback, new_vmax_callback)

# ---- Full-GEM target-parallel worker ----------------------------------------
engine_path = 'R/microcompass_engine.R'
replace_function(engine_path, '.rc_full_gem_step2_model_payload', r'''.rc_full_gem_step2_model_payload <- function(
    model_key, row_ids, model_cache, penalties, vmax_cache,
    omega, solver, flux_threshold, payload_dir) {
  row_ids <- as.character(row_ids)
  if (!length(row_ids)) {
    stop("A full-GEM Step 2 payload requires at least one target row.",
         call. = FALSE)
  }
  first_entry <- model_cache[[row_ids[[1L]]]]
  if (is.null(first_entry)) {
    stop("A full-GEM Step 2 payload references an unknown target row.",
         call. = FALSE)
  }
  model <- .rc_load_microcompass_model(first_entry, "full_gem")
  reactions <- colnames(model$S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions)) ||
      anyDuplicated(reactions)) {
    stop("A full-GEM Step 2 payload has invalid reaction identifiers.",
         call. = FALSE)
  }
  units <- colnames(penalties$penalty)
  if (is.null(units) || anyNA(units) || any(!nzchar(units)) ||
      anyDuplicated(units)) {
    stop("A full-GEM Step 2 payload has invalid unit identifiers.",
         call. = FALSE)
  }
  missing_penalty_rows <- setdiff(reactions, rownames(penalties$penalty))
  if (length(missing_penalty_rows)) {
    stop(
      "Full-GEM Step 2 penalties are missing model reactions: ",
      paste(utils::head(missing_penalty_rows, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  entries <- lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    if (is.null(entry) ||
        !identical(as.character(entry$file), as.character(model_key))) {
      stop("A full-GEM Step 2 payload mixes different model files.",
           call. = FALSE)
    }
    list(
      reaction_id = as.character(entry$reaction_id),
      target_direction = as.character(entry$target_direction),
      medium_scenario = as.character(entry$medium_scenario),
      condition = as.character(entry$condition %||% "all")
    )
  })
  names(entries) <- row_ids
  vmax_values <- lapply(row_ids, function(row_id) {
    value <- vmax_cache[[row_id]]
    list(
      feasible = isTRUE(value$feasible),
      vmax = as.numeric(value$vmax),
      status = as.character(value$status),
      flux = numeric()
    )
  })
  names(vmax_values) <- row_ids
  model_checksum <- if (!is.null(first_entry$file) &&
                        file.exists(first_entry$file)) {
    unname(tools::md5sum(first_entry$file)[[1L]])
  } else {
    NA_character_
  }
  penalty_matrix <- penalties$penalty[reactions, units, drop = FALSE]
  penalty_evidence <- .rc_step2_penalty_evidence_stats(penalty_matrix)
  payload <- list(
    schema_version = "regcompass_full_gem_step2_compact_payload_v1",
    model = list(
      S = .rc_as_dgCMatrix(model$S),
      lb = model$lb,
      ub = model$ub,
      target_status = model$target_status %||% "not_prechecked",
      file_checksum = model_checksum,
      medium_scenario = as.character(first_entry$medium_scenario),
      condition = as.character(first_entry$condition %||% "all")
    ),
    reactions = reactions,
    units = units,
    penalty = penalty_matrix,
    penalty_evidence = penalty_evidence,
    entries = entries,
    vmax = vmax_values,
    omega = as.numeric(omega),
    solver = as.character(solver),
    flux_threshold = as.numeric(flux_threshold)
  )
  token <- substr(.rc_microcompass_object_checksum(list(
    file = as.character(first_entry$file %||% model_key),
    checksum = model_checksum,
    units = units,
    row_ids = row_ids
  )), 1L, 24L)
  file <- file.path(payload_dir, paste0("payload__", token, ".rds"))
  .rc_atomic_save_rds(payload, file)
  rm(model, payload, entries, vmax_values, penalty_matrix, penalty_evidence)
  invisible(gc(verbose = FALSE, full = FALSE))
  file
}''')

replace_function(engine_path, '.rc_full_gem_step2_reaction_batch_worker', r'''.rc_full_gem_step2_reaction_batch_worker <- function(task) {
  if (!is.list(task) ||
      !all(c("payload_file", "row_ids", "checkpoint_dir") %in% names(task))) {
    stop("Malformed full-GEM Step 2 reaction-batch task.", call. = FALSE)
  }
  payload <- readRDS(task$payload_file)
  if (!is.list(payload) ||
      !identical(
        payload$schema_version,
        "regcompass_full_gem_step2_compact_payload_v1"
      )) {
    stop("Malformed full-GEM Step 2 compact payload.", call. = FALSE)
  }
  row_ids <- as.character(task$row_ids)
  if (!length(row_ids) ||
      !all(row_ids %in% names(payload$entries)) ||
      !all(row_ids %in% names(payload$vmax))) {
    stop("Full-GEM Step 2 reaction-batch rows are absent from the payload.",
         call. = FALSE)
  }
  model <- payload$model
  if (!is.list(model) || is.null(model$S) || is.null(model$lb) ||
      is.null(model$ub)) {
    stop("Full-GEM Step 2 compact payload lacks the required LP model state.",
         call. = FALSE)
  }
  if (!identical(colnames(model$S), as.character(payload$reactions))) {
    stop("Full-GEM Step 2 payload reaction order differs from its model.",
         call. = FALSE)
  }
  if (!identical(colnames(payload$penalty), as.character(payload$units)) ||
      !identical(rownames(payload$penalty), as.character(payload$reactions))) {
    stop("Full-GEM Step 2 compact payload penalties are not aligned.",
         call. = FALSE)
  }

  checkpoint_files <- character(length(row_ids))
  step2_engine <- NULL
  on.exit({
    .rc_compass_step2_release_engine(step2_engine)
    rm(model, payload)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  units <- as.character(payload$units)
  evidence <- payload$penalty_evidence %||%
    .rc_step2_penalty_evidence_stats(payload$penalty)
  if (length(evidence$all_finite) != length(units) ||
      length(evidence$fraction) != length(units) ||
      length(evidence$unavailable) != length(units)) {
    stop("Full-GEM Step 2 penalty evidence summary is malformed.",
         call. = FALSE)
  }

  for (j in seq_along(row_ids)) {
    row_id <- row_ids[[j]]
    entry <- payload$entries[[row_id]]
    target_index <- match(entry$reaction_id, payload$reactions)
    if (is.na(target_index)) {
      stop("A full-GEM target reaction is absent from its shared model.",
           call. = FALSE)
    }
    prepared <- .rc_compass_step2_prepare(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_reaction = entry$reaction_id,
      vmax_result = payload$vmax[[row_id]],
      target_direction = entry$target_direction,
      omega = payload$omega,
      flux_threshold = payload$flux_threshold
    )
    step2_engine <- if (isTRUE(prepared$runnable)) {
      .rc_compass_step2_new_engine(
        prepared$template, payload$solver,
        persistent_required = identical(payload$solver, "highs")
      )
    } else {
      NULL
    }

    n_units <- length(units)
    task_penalty <- rep(NA_real_, n_units)
    task_vmax <- rep(NA_real_, n_units)
    task_feasible <- task_evaluated <- rep(FALSE, n_units)
    names(task_penalty) <- names(task_vmax) <-
      names(task_feasible) <- names(task_evaluated) <- units
    solver_status <- step1_status <- step2_status <-
      solver_backend <- rep(NA_character_, n_units)
    target_available <- is.finite(payload$penalty[target_index, ])

    for (i in seq_along(units)) {
      unit_penalty <- payload$penalty[, i]
      solver_penalty <- if (isTRUE(evidence$all_finite[[i]])) {
        unit_penalty
      } else {
        value <- unit_penalty
        value[!is.finite(value)] <- 0
        value
      }

      if (isTRUE(prepared$runnable)) {
        solved <- .rc_compass_step2_engine_solve(
          step2_engine, solver_penalty,
          return_solution = FALSE,
          trusted_aligned = TRUE
        )
        step2_engine <- solved$engine
        fit <- .rc_compass_step2_result(
          prepared$template, solved$answer,
          require_solution = FALSE
        )
      } else {
        fit <- prepared$result
        solved <- NULL
      }

      task_penalty[[i]] <- if (target_available[[i]]) {
        fit$penalty
      } else {
        NA_real_
      }
      task_vmax[[i]] <- fit$vmax
      task_feasible[[i]] <- isTRUE(fit$feasible)
      task_evaluated[[i]] <- isTRUE(fit$feasible) && target_available[[i]]
      solver_status[[i]] <- as.character(fit$solver_status)
      solver_backend[[i]] <- as.character(fit$solver_backend %||% "unknown")
      step1_status[[i]] <- as.character(fit$step1_status)
      step2_status[[i]] <- as.character(fit$step2_status)
      rm(unit_penalty, solver_penalty, fit, solved)
    }

    diagnostics <- data.frame(
      row_id = rep(row_id, n_units),
      unit_id = units,
      module_id = rep(NA_character_, n_units),
      reaction_id = rep(entry$reaction_id, n_units),
      target_direction = rep(entry$target_direction, n_units),
      medium_scenario = rep(entry$medium_scenario, n_units),
      condition = rep(entry$condition, n_units),
      strict_feasible = unname(task_feasible),
      solver_status = solver_status,
      solver_backend = solver_backend,
      step1_status = step1_status,
      step2_status = step2_status,
      target_status = ifelse(
        task_feasible, "ok", "medium_directionally_infeasible"
      ),
      objective_value = unname(task_penalty),
      vmax = unname(task_vmax),
      vmax_reused_from_shared_cache = rep(TRUE, n_units),
      step2_model_reused_across_metacells = rep(TRUE, n_units),
      target_expression_available = target_available,
      objective_evidence_fraction = evidence$fraction,
      unavailable_objective_terms = evidence$unavailable,
      parallel_task = rep(
        "directional_reaction_x_all_metacells", n_units
      ),
      stringsAsFactors = FALSE
    )

    engine_metrics <- .rc_compass_step2_engine_metrics(step2_engine)
    token <- substr(.rc_microcompass_object_checksum(list(
      row_id = row_id,
      file_checksum = model$file_checksum,
      units = units,
      omega = payload$omega,
      solver = payload$solver,
      flux_threshold = payload$flux_threshold
    )), 1L, 24L)
    checkpoint <- file.path(
      task$checkpoint_dir, paste0("step2__", token, ".rds")
    )
    .rc_atomic_save_rds(list(
      row_id = row_id,
      units = units,
      penalty = task_penalty,
      vmax = task_vmax,
      feasible = task_feasible,
      evaluated = task_evaluated,
      diagnostics = diagnostics,
      engine_metrics = engine_metrics
    ), checkpoint)
    checkpoint_files[[j]] <- checkpoint

    .rc_compass_step2_release_engine(step2_engine)
    step2_engine <- NULL
    rm(
      diagnostics, task_penalty, task_vmax, task_feasible,
      task_evaluated, prepared, engine_metrics, solver_status,
      solver_backend, step1_status, step2_status, target_available
    )
    invisible(gc(verbose = FALSE, full = FALSE))
  }
  checkpoint_files
}''')

# ---- Cell-type union-GEM target worker: same exact fast path -----------------
cell_path = 'R/celltype_microcompass_reaction_parallel.R'
replace_function(cell_path, '.rc_step2_model_payload', r'''.rc_step2_model_payload <- function(
    model_key, row_ids, model_cache, unit_celltype, penalties, vmax_cache,
    omega, solver, flux_threshold, payload_dir) {
  row_ids <- as.character(row_ids)
  first_entry <- model_cache[[row_ids[[1L]]]]
  eligible <- names(unit_celltype)[unit_celltype == first_entry$cell_type]
  if (!length(eligible)) {
    stop("No Layer 1 units match cell type `", first_entry$cell_type, "`.",
         call. = FALSE)
  }
  model <- .rc_read_celltype_union_gem(
    first_entry$file, first_entry$cell_type,
    first_entry$medium_scenario, first_entry$file_checksum
  )
  reactions <- colnames(model$S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions))) {
    stop("A cell-type union GEM has invalid reaction identifiers.",
         call. = FALSE)
  }
  entries <- lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    if (!identical(as.character(entry$file), as.character(model_key))) {
      stop("A Step 2 payload mixes different union GEM files.", call. = FALSE)
    }
    list(
      reaction_id = as.character(entry$reaction_id),
      target_direction = as.character(entry$target_direction),
      cell_type = as.character(entry$cell_type),
      medium_scenario = as.character(entry$medium_scenario)
    )
  })
  names(entries) <- row_ids
  vmax_values <- lapply(row_ids, function(row_id) {
    .rc_step2_compact_vmax_value(vmax_cache[[row_id]])
  })
  names(vmax_values) <- row_ids
  penalty_matrix <- penalties$penalty[reactions, eligible, drop = FALSE]
  penalty_evidence <- .rc_step2_penalty_evidence_stats(penalty_matrix)
  payload <- list(
    schema_version = "regcompass_step2_compact_payload_v1",
    model = list(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_status = model$target_status %||% NA_character_,
      file_checksum = as.character(first_entry$file_checksum),
      cell_type = as.character(first_entry$cell_type),
      medium_scenario = as.character(first_entry$medium_scenario)
    ),
    reactions = reactions,
    units = eligible,
    penalty = penalty_matrix,
    penalty_evidence = penalty_evidence,
    entries = entries,
    vmax = vmax_values,
    omega = as.numeric(omega),
    solver = as.character(solver),
    flux_threshold = as.numeric(flux_threshold)
  )
  token <- substr(.rc_microcompass_object_checksum(list(
    file = first_entry$file,
    checksum = first_entry$file_checksum,
    units = eligible,
    row_ids = row_ids
  )), 1L, 24L)
  file <- file.path(payload_dir, paste0("payload__", token, ".rds"))
  .rc_atomic_save_rds(payload, file)
  rm(model, payload, entries, vmax_values, penalty_matrix, penalty_evidence)
  invisible(gc(verbose = FALSE, full = FALSE))
  file
}''')

replace_function(cell_path, '.rc_step2_reaction_batch_worker', r'''.rc_step2_reaction_batch_worker <- function(task) {
  if (!is.list(task) ||
      !all(c("payload_file", "row_ids", "checkpoint_dir") %in% names(task))) {
    stop("Malformed Step 2 reaction-batch task.", call. = FALSE)
  }
  payload <- readRDS(task$payload_file)
  if (!is.list(payload) ||
      !identical(payload$schema_version, "regcompass_step2_compact_payload_v1")) {
    stop("Malformed Step 2 compact payload.", call. = FALSE)
  }
  row_ids <- as.character(task$row_ids)
  if (!length(row_ids) || !all(row_ids %in% names(payload$entries)) ||
      !all(row_ids %in% names(payload$vmax))) {
    stop("Step 2 reaction-batch rows are absent from the compact payload.",
         call. = FALSE)
  }
  model <- payload$model
  if (!is.list(model) || is.null(model$S) || is.null(model$lb) ||
      is.null(model$ub)) {
    stop("Step 2 compact payload lacks the required LP model state.",
         call. = FALSE)
  }
  if (!identical(colnames(model$S), as.character(payload$reactions))) {
    stop("Step 2 compact payload reaction order differs from its union GEM.",
         call. = FALSE)
  }
  if (!identical(colnames(payload$penalty), as.character(payload$units)) ||
      !identical(rownames(payload$penalty), as.character(payload$reactions))) {
    stop("Step 2 compact payload penalties are not aligned.", call. = FALSE)
  }
  checkpoint_files <- character(length(row_ids))
  step2_engine <- NULL
  on.exit({
    .rc_compass_step2_release_engine(step2_engine)
    rm(model, payload)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  units <- as.character(payload$units)
  evidence <- payload$penalty_evidence %||%
    .rc_step2_penalty_evidence_stats(payload$penalty)
  if (length(evidence$all_finite) != length(units) ||
      length(evidence$fraction) != length(units) ||
      length(evidence$unavailable) != length(units)) {
    stop("Step 2 penalty evidence summary is malformed.", call. = FALSE)
  }

  for (j in seq_along(row_ids)) {
    row_id <- row_ids[[j]]
    entry <- payload$entries[[row_id]]
    target_index <- match(entry$reaction_id, payload$reactions)
    if (is.na(target_index)) {
      stop("A target reaction is absent from its cell-type union GEM.",
           call. = FALSE)
    }
    prepared <- .rc_compass_step2_prepare(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_reaction = entry$reaction_id,
      vmax_result = payload$vmax[[row_id]],
      target_direction = entry$target_direction,
      omega = payload$omega,
      flux_threshold = payload$flux_threshold
    )
    step2_engine <- if (isTRUE(prepared$runnable)) {
      .rc_compass_step2_new_engine(
        prepared$template, payload$solver,
        persistent_required = identical(payload$solver, "highs")
      )
    } else {
      NULL
    }

    n_units <- length(units)
    task_penalty <- rep(NA_real_, n_units)
    task_vmax <- rep(NA_real_, n_units)
    task_feasible <- task_evaluated <- rep(FALSE, n_units)
    names(task_penalty) <- names(task_vmax) <-
      names(task_feasible) <- names(task_evaluated) <- units
    solver_status <- step1_status <- step2_status <-
      solver_backend <- rep(NA_character_, n_units)
    target_available <- is.finite(payload$penalty[target_index, ])

    for (i in seq_along(units)) {
      unit_penalty <- payload$penalty[, i]
      solver_penalty <- if (isTRUE(evidence$all_finite[[i]])) {
        unit_penalty
      } else {
        value <- unit_penalty
        value[!is.finite(value)] <- 0
        value
      }

      if (isTRUE(prepared$runnable)) {
        solved <- .rc_compass_step2_engine_solve(
          step2_engine, solver_penalty,
          return_solution = FALSE,
          trusted_aligned = TRUE
        )
        step2_engine <- solved$engine
        fit <- .rc_compass_step2_result(
          prepared$template, solved$answer,
          require_solution = FALSE
        )
      } else {
        .rc_compass_step2_align_penalties(
          prepared$reactions, solver_penalty
        )
        fit <- prepared$result
        solved <- NULL
      }

      task_penalty[[i]] <- if (target_available[[i]]) fit$penalty else NA_real_
      task_vmax[[i]] <- fit$vmax
      task_feasible[[i]] <- isTRUE(fit$feasible)
      task_evaluated[[i]] <- isTRUE(fit$feasible) && target_available[[i]]
      solver_status[[i]] <- as.character(fit$solver_status)
      solver_backend[[i]] <- as.character(fit$solver_backend %||% "unknown")
      step1_status[[i]] <- as.character(fit$step1_status)
      step2_status[[i]] <- as.character(fit$step2_status)
      rm(unit_penalty, solver_penalty, fit, solved)
    }

    target_status <- if (!is.null(model$target_status)) {
      rep(as.character(model$target_status), n_units)
    } else {
      ifelse(task_feasible, "ok", "structurally_infeasible")
    }
    diagnostics <- data.frame(
      row_id = rep(row_id, n_units),
      unit_id = units,
      module_id = rep("CELLTYPE_MEDIUM_UNION_GEM", n_units),
      cell_type = rep(entry$cell_type, n_units),
      reaction_id = rep(entry$reaction_id, n_units),
      target_direction = rep(entry$target_direction, n_units),
      medium_scenario = rep(entry$medium_scenario, n_units),
      condition = rep("all", n_units),
      strict_feasible = unname(task_feasible),
      solver_status = solver_status,
      solver_backend = solver_backend,
      step1_status = step1_status,
      step2_status = step2_status,
      target_status = target_status,
      objective_value = unname(task_penalty),
      vmax = unname(task_vmax),
      vmax_reused_from_celltype_cache = rep(TRUE, n_units),
      step2_model_reused_across_metacells = rep(TRUE, n_units),
      target_expression_available = target_available,
      objective_evidence_fraction = evidence$fraction,
      unavailable_objective_terms = evidence$unavailable,
      parallel_task = rep(
        "directional_reaction_x_matching_metacells", n_units
      ),
      stringsAsFactors = FALSE
    )

    engine_metrics <- .rc_compass_step2_engine_metrics(step2_engine)
    token <- substr(.rc_microcompass_object_checksum(list(
      row_id = row_id,
      file_checksum = model$file_checksum,
      units = units,
      omega = payload$omega,
      solver = payload$solver,
      flux_threshold = payload$flux_threshold
    )), 1L, 24L)
    checkpoint <- file.path(
      task$checkpoint_dir, paste0("step2__", token, ".rds")
    )
    .rc_atomic_save_rds(list(
      row_id = row_id,
      units = units,
      penalty = task_penalty,
      vmax = task_vmax,
      feasible = task_feasible,
      evaluated = task_evaluated,
      diagnostics = diagnostics,
      engine_metrics = engine_metrics
    ), checkpoint)
    checkpoint_files[[j]] <- checkpoint

    .rc_compass_step2_release_engine(step2_engine)
    step2_engine <- NULL
    rm(
      diagnostics, task_penalty, task_vmax, task_feasible,
      task_evaluated, prepared, engine_metrics, solver_status,
      solver_backend, step1_status, step2_status, target_available,
      target_status
    )
    invisible(gc(verbose = FALSE, full = FALSE))
  }
  checkpoint_files
}''')

# ---- Package contract / CI guards ------------------------------------------
desc = Path('DESCRIPTION')
text = desc.read_text()
text = text.replace('Version: 2.4.22', 'Version: 2.4.23', 1)
text = text.replace('    highs,\n', '    highs (>= 1.12.0-1),\n', 1)
desc.write_text(text)

workflow = Path('.github/workflows/full-gem-step2-checks.yaml')
wtext = workflow.read_text()
needle = "              'step2_parallel_workers = as.integer(step2_workers)',\n"
addition = (
    needle
    + "              'return_solution = FALSE',\n"
    + "              'require_solution = FALSE',\n"
    + "              'persistent_required = identical(payload$solver, \"highs\")',\n"
)
if "'return_solution = FALSE'" not in wtext:
    if needle not in wtext:
        raise RuntimeError('Cannot extend full-GEM fast-path audit')
    wtext = wtext.replace(needle, addition, 1)
workflow.write_text(wtext)

# ---- Numerical regressions --------------------------------------------------
test_path = Path('tests/layer2-native-acceleration-check.R')
ttext = test_path.read_text()
marker = '  metrics <- .rc_compass_step2_engine_metrics(engine)\n'
score_check = r'''  score_only_solved <- .rc_compass_step2_engine_solve(
    engine, penalty_sets[[2L]],
    return_solution = FALSE,
    trusted_aligned = TRUE
  )
  engine <- score_only_solved$engine
  score_only <- .rc_compass_step2_result(
    prepared$template, score_only_solved$answer,
    require_solution = FALSE
  )
  score_reference <- rc_compass_two_step_lp_directional(
    S, lb, ub, "TARGET", penalty_sets[[2L]],
    target_direction = "forward", omega = 0.95, solver = "highs"
  )
  stopifnot(
    length(score_only_solved$answer$solution) == 0L,
    identical(score_only$feasible, score_reference$feasible),
    isTRUE(all.equal(
      score_only$penalty, score_reference$penalty, tolerance = 1e-9
    )),
    length(score_only$flux) == 0L
  )

'''
if 'score_only_solved <-' not in ttext:
    if marker not in ttext:
        raise RuntimeError('Cannot insert objective-only Step 2 regression')
    ttext = ttext.replace(marker, score_check + marker, 1)

vmax_marker = 'run_persistent_step2_check <- function() {\n'
vmax_check = r'''run_persistent_vmax_check <- function() {
  if (!.rc_microcompass_highs_api_available()) return(invisible(NULL))
  S <- Matrix::Matrix(
    matrix(
      c(1, -1, -1), nrow = 1,
      dimnames = list("M", c("UP", "T1", "T2"))
    ),
    sparse = TRUE
  )
  lb <- c(UP = 0, T1 = 0, T2 = 0)
  ub <- c(UP = 10, T1 = 10, T2 = 10)
  batch <- .rc_compass_vmax_batch_highs(
    S, lb, ub,
    target_reaction = c("T1", "T2"),
    direction = c("forward", "forward"),
    flux_threshold = 1e-8
  )
  stopifnot(length(batch) == 2L)
  for (i in seq_along(batch)) {
    target <- c("T1", "T2")[[i]]
    reference <- rc_compass_vmax_directional(
      S, lb, ub, target, direction = "forward", solver = "highs"
    )
    stopifnot(
      identical(batch[[i]]$feasible, reference$feasible),
      identical(batch[[i]]$status, reference$status),
      isTRUE(all.equal(batch[[i]]$vmax, reference$vmax, tolerance = 1e-10)),
      length(batch[[i]]$flux) == 0L
    )
  }
  invisible(batch)
}

'''
if 'run_persistent_vmax_check <- function()' not in ttext:
    if vmax_marker not in ttext:
        raise RuntimeError('Cannot insert persistent Vmax regression')
    ttext = ttext.replace(vmax_marker, vmax_check + vmax_marker, 1)

call_marker = 'run_persistent_step2_check()\nrun_compact_step2_worker_check()\n'
if 'run_persistent_vmax_check()\nrun_persistent_step2_check()' not in ttext:
    if call_marker not in ttext:
        raise RuntimeError('Cannot insert persistent Vmax test call')
    ttext = ttext.replace(
        call_marker,
        'run_persistent_vmax_check()\n' + call_marker,
        1
    )
test_path.write_text(ttext)

# Testthat regression for package-level CI.
testthat_path = Path('tests/testthat/test-layer2-native-acceleration.R')
tt = testthat_path.read_text()
append = r'''

test_that("objective-only persistent Step 2 preserves the canonical penalty", {
  skip_if_not_installed("highs")
  S <- Matrix::Matrix(
    matrix(c(1, -1), nrow = 1,
           dimnames = list("M", c("UP", "TARGET"))),
    sparse = TRUE
  )
  lb <- c(UP = 0, TARGET = 0)
  ub <- c(UP = 10, TARGET = 10)
  penalties <- c(UP = 0.75, TARGET = 0.10)
  vmax <- rc_compass_vmax_directional(
    S, lb, ub, "TARGET", direction = "forward", solver = "highs"
  )
  prepared <- RegCompassR:::.rc_compass_step2_prepare(
    S, lb, ub, "TARGET", vmax,
    target_direction = "forward", omega = 0.95
  )
  engine <- RegCompassR:::.rc_compass_step2_new_engine(
    prepared$template, "highs", persistent_required = TRUE
  )
  on.exit(RegCompassR:::.rc_compass_step2_release_engine(engine), add = TRUE)
  solved <- RegCompassR:::.rc_compass_step2_engine_solve(
    engine, penalties, return_solution = FALSE, trusted_aligned = TRUE
  )
  engine <- solved$engine
  observed <- RegCompassR:::.rc_compass_step2_result(
    prepared$template, solved$answer, require_solution = FALSE
  )
  reference <- rc_compass_two_step_lp_directional(
    S, lb, ub, "TARGET", penalties,
    target_direction = "forward", omega = 0.95, solver = "highs"
  )
  expect_length(solved$answer$solution, 0L)
  expect_identical(observed$feasible, reference$feasible)
  expect_equal(observed$penalty, reference$penalty, tolerance = 1e-9)
  expect_length(observed$flux, 0L)
})

test_that("persistent batched Vmax matches independent directional LPs", {
  skip_if_not_installed("highs")
  S <- Matrix::Matrix(
    matrix(c(1, -1, -1), nrow = 1,
           dimnames = list("M", c("UP", "T1", "T2"))),
    sparse = TRUE
  )
  lb <- c(UP = 0, T1 = 0, T2 = 0)
  ub <- c(UP = 10, T1 = 10, T2 = 10)
  observed <- RegCompassR:::.rc_compass_vmax_batch_highs(
    S, lb, ub,
    target_reaction = c("T1", "T2"),
    direction = c("forward", "forward"),
    flux_threshold = 1e-8
  )
  skip_if(is.null(observed), "Persistent HiGHS API unavailable")
  for (i in seq_along(observed)) {
    target <- c("T1", "T2")[[i]]
    reference <- rc_compass_vmax_directional(
      S, lb, ub, target, direction = "forward", solver = "highs"
    )
    expect_identical(observed[[i]]$feasible, reference$feasible)
    expect_identical(observed[[i]]$status, reference$status)
    expect_equal(observed[[i]]$vmax, reference$vmax, tolerance = 1e-10)
    expect_length(observed[[i]]$flux, 0L)
  }
})
'''
if 'objective-only persistent Step 2 preserves the canonical penalty' not in tt:
    tt = tt.rstrip() + append + '\n'
testthat_path.write_text(tt)

print('Applied exact target-parallel Layer 2 speedups.')
