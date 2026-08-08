suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
  library(BiocParallel)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
.rc_bind_frames_fill <- function(values) {
  values <- values[vapply(values, is.data.frame, logical(1))]
  values <- values[vapply(values, nrow, integer(1)) > 0L]
  if (!length(values)) return(data.frame())
  columns <- unique(unlist(lapply(values, colnames), use.names = FALSE))
  values <- lapply(values, function(value) {
    missing <- setdiff(columns, colnames(value))
    for (name in missing) value[[name]] <- NA
    value[, columns, drop = FALSE]
  })
  answer <- do.call(rbind, values)
  rownames(answer) <- NULL
  answer
}
.rc_lp_status <- function(message = "", code = NA_integer_) {
  text <- tolower(paste(message, collapse = " "))
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("time|limit", text)) return("time_limit")
  if (grepl("optimal", text)) return("optimal")
  if (is.finite(code) && as.integer(code) == 0L) return("optimal")
  "error"
}
rc_validate_gem <- function(gem) {
  S <- .rc_as_dgCMatrix(gem$S)
  reactions <- colnames(S)
  lb <- as.numeric(gem$lb[reactions])
  ub <- as.numeric(gem$ub[reactions])
  names(lb) <- names(ub) <- reactions
  list(S = S, lb = lb, ub = ub, reactions = reactions)
}
rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = "highs", time_limit = Inf) {
  answer <- highs::highs_solve(
    L = as.numeric(obj), lower = as.numeric(lb), upper = as.numeric(ub),
    A = A, lhs = as.numeric(lhs), rhs = as.numeric(rhs), maximum = FALSE,
    control = highs::highs_control(
      log_to_console = FALSE, output_flag = FALSE,
      threads = 1L, solver = "simplex",
      primal_feasibility_tolerance = 1e-7,
      time_limit = as.numeric(time_limit)
    )
  )
  list(
    status = .rc_lp_status(answer$status_message, answer$status),
    solution = as.numeric(answer$primal_solution),
    objective = as.numeric(answer$objective_value),
    solver_message = answer$status_message
  )
}
.rc_progress_enabled <- function(x) FALSE
.rc_format_elapsed <- function(seconds) sprintf("%.3fs", seconds)
.rc_safe_cache_token <- function(value) {
  paste(sprintf("%02x", as.integer(charToRaw(as.character(value)))), collapse = "")
}
.rc_layer2_task_context <- function(
    cell_type = "ALL", medium_scenario = "base", route = "unknown") {
  list(cell_type = cell_type, medium_scenario = medium_scenario, route = route)
}
.rc_layer2_task_event <- function(...) invisible(NULL)
.rc_layer2_progress_dir_from_cache <- function(...) tempdir()
.rc_layer2_algorithm_once <- function(...) invisible(NULL)
.rc_layer2_progress_state <- new.env(parent = emptyenv())
.rc_layer2_progress_state$current_task <- NULL
.rc_layer2_progress_state$inside_dependency <- FALSE
.rc_layer2_progress_state$algorithm_flags <- new.env(parent = emptyenv())
.rc_layer2_parallel_context <- new.env(parent = emptyenv())
.rc_layer2_parallel_context$active <- FALSE
.rc_layer2_parallel_context$parallel <- FALSE
.rc_layer2_parallel_context$BPPARAM <- FALSE
.rc_layer2_parallel_context$nested_serial <- FALSE
.rc_layer2_enter_parallel_context <- function(parallel, BPPARAM) {
  previous <- as.list(.rc_layer2_parallel_context)
  .rc_layer2_parallel_context$active <- TRUE
  .rc_layer2_parallel_context$parallel <- isTRUE(parallel)
  .rc_layer2_parallel_context$BPPARAM <- BPPARAM
  .rc_layer2_parallel_context$nested_serial <- FALSE
  previous
}
.rc_layer2_restore_parallel_context <- function(previous) {
  rm(list = ls(.rc_layer2_parallel_context, all.names = TRUE),
     envir = .rc_layer2_parallel_context)
  list2env(previous, envir = .rc_layer2_parallel_context)
  invisible(NULL)
}
rc_available_workers <- function(default = 1L) 2L
rc_default_bpparam <- function(workers = NULL, backend = "auto") {
  workers <- as.integer(workers %||% 2L)
  BiocParallel::SnowParam(
    workers = workers, type = "SOCK", progressbar = FALSE,
    exportglobals = TRUE, exportvariables = TRUE
  )
}
.rc_layer2_tune_task_bpparam <- function(BPPARAM, n_tasks) BPPARAM
.rc_release_bpparam <- function(param) {
  if (!identical(param, FALSE) && !is.null(param) &&
      isTRUE(BiocParallel::bpisup(param))) {
    try(BiocParallel::bpstop(param), silent = TRUE)
  }
  invisible(NULL)
}
.rc_with_internal_single_thread <- function(FUN) FUN()
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) return(lapply(X, FUN, ...))
  was_started <- isTRUE(BiocParallel::bpisup(BPPARAM))
  if (!was_started) {
    BiocParallel::bpstart(BPPARAM)
    on.exit(.rc_release_bpparam(BPPARAM), add = TRUE)
  }
  names_all <- ls(.GlobalEnv, all.names = TRUE)
  runtime_functions <- mget(names_all, envir = .GlobalEnv, inherits = FALSE)
  runtime_functions <- runtime_functions[vapply(
    runtime_functions, is.function, logical(1)
  )]
  extra <- list(...)
  worker <- function(x, FUN, runtime_functions, extra) {
    list2env(runtime_functions, envir = .GlobalEnv)
    do.call(FUN, c(list(x), extra))
  }
  BiocParallel::bplapply(
    X, worker,
    FUN = FUN,
    runtime_functions = runtime_functions,
    extra = extra,
    BPPARAM = BPPARAM
  )
}

source("R/layer2_corda_evidence.R")
source("R/layer2_corda_lp.R")
source("R/layer2_corda_paper_contract.R")
source("R/layer2_corda_direction_contract.R")
source("R/layer2_corda2_algorithm.R")
source("R/layer2_corda_stage_parallel.R")
source("R/layer2_corda2_algorithm_build.R")
source("R/layer2_corda2_options_contract.R")
source("R/layer2_corda_runtime.R")

metabolites <- c("A1", "B1", "A2", "B2", "C1", "C2")
reactions <- c(
  "M1", "H1", "N1", "M2", "H2", "N2",
  "M3", "N3", "M4", "N4"
)
S <- Matrix::sparseMatrix(
  i = integer(), j = integer(), x = numeric(),
  dims = c(length(metabolites), length(reactions)),
  dimnames = list(metabolites, reactions), giveCsparse = TRUE
)
S["A1", "M1"] <- 1; S["A1", "H1"] <- -1
S["B1", "H1"] <- 1; S["B1", "N1"] <- -1
S["A2", "M2"] <- 1; S["A2", "H2"] <- -1
S["B2", "H2"] <- 1; S["B2", "N2"] <- -1
S["C1", "M3"] <- 1; S["C1", "N3"] <- -1
S["C2", "M4"] <- 1; S["C2", "N4"] <- -1
split <- .rc_corda2_split_original(list(
  S = S,
  lb = stats::setNames(rep(0, length(reactions)), reactions),
  ub = stats::setNames(rep(1000, length(reactions)), reactions)
))
initial <- stats::setNames(rep("OT", length(reactions)), reactions)
initial[c("H1", "H2")] <- "HC"
initial[c("M1", "M2", "M3", "M4")] <- "MC_module"
initial[c("N1", "N2", "N3", "N4")] <- "NC"
classes <- list(
  hc = c("H1", "H2"),
  mc_module = c("M1", "M2", "M3", "M4"),
  mc_evidence = character(),
  mc = c("M1", "M2", "M3", "M4"),
  nc = c("N1", "N2", "N3", "N4"),
  ot = character(),
  confidence = initial,
  initial_confidence = initial
)
options <- .rc_layer2_corda_options(list(model_completion = "corda2"))

previous <- .rc_layer2_enter_parallel_context(FALSE, FALSE)
serial <- tryCatch(
  .rc_corda_build_three_stage_core(
    split, classes, options, solver = "highs", time_limit = 30
  ),
  finally = .rc_layer2_restore_parallel_context(previous)
)

param <- BiocParallel::SnowParam(
  workers = 2L, type = "SOCK", progressbar = FALSE,
  exportglobals = TRUE, exportvariables = TRUE
)
previous <- .rc_layer2_enter_parallel_context(TRUE, param)
parallel <- tryCatch(
  .rc_corda_build_three_stage_core(
    split, classes, options, solver = "highs", time_limit = 30
  ),
  finally = .rc_layer2_restore_parallel_context(previous)
)

stopifnot(
  setequal(serial$included, parallel$included),
  setequal(serial$included_directional_variables,
           parallel$included_directional_variables),
  isTRUE(all.equal(serial$HCtoMC, parallel$HCtoMC)),
  isTRUE(all.equal(serial$HCtoNC, parallel$HCtoNC)),
  isTRUE(all.equal(serial$MCtoNC, parallel$MCtoNC)),
  setequal(serial$stage1_associated, parallel$stage1_associated),
  setequal(serial$stage2_promoted_nc, parallel$stage2_promoted_nc),
  setequal(serial$stage2_promoted_mc, parallel$stage2_promoted_mc),
  setequal(serial$stage3_associated_ot, parallel$stage3_associated_ot),
  identical(serial$stage_update_policy, "original_matlab_directional_order"),
  identical(parallel$stage_update_policy, "original_matlab_directional_order"),
  identical(
    parallel$parallel_execution_policy,
    "stage_barrier_parallel_targets_deterministic_ordered_reduce"
  ),
  isTRUE(parallel$solver_performance$stage_barrier),
  identical(parallel$solver_performance$target_parallelism,
            "within_corda2_stage"),
  !BiocParallel::bpisup(param)
)

cat("CORDA2 stage-parallel serial-equivalence check passed\n")
