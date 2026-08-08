suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
  library(BiocParallel)
})

e <- new.env(parent = baseenv())
e$`%||%` <- function(x, y) if (is.null(x)) y else x
e$.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
e$.rc_bind_frames_fill <- function(values) {
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
e$.rc_lp_status <- function(message = "", code = NA_integer_) {
  text <- tolower(paste(message, collapse = " "))
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("time|limit", text)) return("time_limit")
  if (grepl("optimal", text)) return("optimal")
  if (is.finite(code) && as.integer(code) == 0L) return("optimal")
  "error"
}
e$rc_validate_gem <- function(gem) {
  S <- .rc_as_dgCMatrix(gem$S)
  reactions <- colnames(S)
  lb <- as.numeric(gem$lb[reactions])
  ub <- as.numeric(gem$ub[reactions])
  names(lb) <- names(ub) <- reactions
  if (anyNA(lb) || anyNA(ub) || any(lb > ub)) stop("invalid GEM bounds")
  list(S = S, lb = lb, ub = ub, reactions = reactions)
}
e$rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
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
e$.rc_progress_enabled <- function(x) FALSE
e$.rc_release_bpparam <- function(param) {
  if (isTRUE(BiocParallel::bpisup(param))) BiocParallel::bpstop(param)
  invisible(NULL)
}
e$rc_available_workers <- function(default = 1L) 2L
e$rc_parallel_config <- function(...) list(workers = 2L)
e$rc_default_bpparam <- function(workers = NULL, backend = "auto") {
  workers <- as.integer(workers %||% 2L)
  if (workers <= 1L) return(NULL)
  BiocParallel::SnowParam(workers = workers, type = "SOCK", progressbar = FALSE)
}
e$.rc_layer2_tune_task_bpparam <- function(BPPARAM, n_tasks) BPPARAM
e$.rc_layer2_pool_workers <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) return(1L)
  as.integer(BiocParallel::bpnworkers(BPPARAM))
}
e$rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) {
    return(lapply(X, FUN, ...))
  }
  was_started <- isTRUE(BiocParallel::bpisup(BPPARAM))
  if (!was_started) {
    BiocParallel::bpstart(BPPARAM)
    on.exit({
      BiocParallel::bpstop(BPPARAM)
      invisible(gc(verbose = FALSE, full = TRUE))
    }, add = TRUE)
  }
  BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)
}
e$.rc_layer2_parallel_context <- new.env(parent = emptyenv())
e$.rc_layer2_parallel_context$active <- FALSE
e$.rc_layer2_parallel_context$parallel <- FALSE
e$.rc_layer2_parallel_context$BPPARAM <- FALSE
e$.rc_layer2_parallel_context$nested_serial <- FALSE
e$.rc_layer2_enter_parallel_context <- function(parallel, BPPARAM) {
  previous <- as.list(.rc_layer2_parallel_context)
  .rc_layer2_parallel_context$active <- TRUE
  .rc_layer2_parallel_context$parallel <- isTRUE(parallel)
  .rc_layer2_parallel_context$BPPARAM <- BPPARAM
  .rc_layer2_parallel_context$nested_serial <- FALSE
  previous
}
e$.rc_layer2_restore_parallel_context <- function(previous) {
  rm(list = ls(.rc_layer2_parallel_context, all.names = TRUE),
     envir = .rc_layer2_parallel_context)
  list2env(previous, envir = .rc_layer2_parallel_context)
  invisible(NULL)
}

e$e <- e
for (name in ls(e, all.names = TRUE)) {
  if (is.function(e[[name]])) environment(e[[name]]) <- e
}

for (file in c(
  "R/layer2_corda_lp.R",
  "R/layer2_corda_paper_contract.R",
  "R/layer2_corda_direction_contract.R",
  "R/layer2_corda2_algorithm.R",
  "R/layer2_corda_serial_core.R",
  "R/layer2_corda_stage_parallel.R",
  "R/layer2_corda2_algorithm_build.R",
  "R/layer2_corda_dispatch.R",
  "R/layer2_corda2_options_contract.R",
  "R/layer2_corda_runtime.R"
)) source(file, local = e)

metabolites <- c("A1", "B1", "A2", "B2", "C1", "C2")
reactions <- c("M1", "H1", "N1", "M2", "H2", "N2", "M3", "N3", "M4", "N4")
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

gem <- list(
  S = S,
  lb = stats::setNames(rep(0, length(reactions)), reactions),
  ub = stats::setNames(rep(1000, length(reactions)), reactions)
)
split <- e$.rc_corda2_split_original(gem)
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
  ot = character(), confidence = initial, initial_confidence = initial
)
corda_options <- e$.rc_layer2_corda_options(list(model_completion = "corda2"))

serial <- e$.rc_corda_build_three_stage_serial_core(
  split = split, classes = classes, options = corda_options,
  solver = "highs", time_limit = 30
)

param <- BiocParallel::SnowParam(workers = 2L, type = "SOCK", progressbar = FALSE)
previous <- e$.rc_layer2_enter_parallel_context(TRUE, param)
parallel_result <- e$.rc_corda_build_three_stage_core(
  split = split, classes = classes, options = corda_options,
  solver = "highs", time_limit = 30
)
e$.rc_layer2_restore_parallel_context(previous)

stopifnot(
  setequal(parallel_result$included, serial$included),
  setequal(
    parallel_result$included_directional_variables,
    serial$included_directional_variables
  ),
  identical(parallel_result$HCtoMC, serial$HCtoMC),
  identical(parallel_result$HCtoNC, serial$HCtoNC),
  identical(parallel_result$MCtoNC, serial$MCtoNC),
  setequal(parallel_result$stage1_associated, serial$stage1_associated),
  setequal(parallel_result$stage2_promoted_nc, serial$stage2_promoted_nc),
  setequal(parallel_result$stage2_promoted_mc, serial$stage2_promoted_mc),
  setequal(parallel_result$stage3_associated_ot, serial$stage3_associated_ot),
  identical(serial$parallel_execution_policy,
            "serial_original_persistent_engine"),
  identical(
    parallel_result$stage_update_policy,
    "original_matlab_directional_order"
  ),
  identical(
    parallel_result$parallel_execution_policy,
    "stage_barrier_parallel_targets_deterministic_ordered_reduce"
  ),
  isTRUE(parallel_result$solver_performance$stage_barrier),
  identical(
    parallel_result$solver_performance$target_parallelism,
    "within_corda2_stage"
  ),
  !isTRUE(BiocParallel::bpisup(param))
)

previous <- e$.rc_layer2_enter_parallel_context(FALSE, FALSE)
dispatched_serial <- e$.rc_corda_build_three_stage_dispatch(
  split = split, classes = classes, options = corda_options,
  solver = "highs", time_limit = 30
)
e$.rc_layer2_restore_parallel_context(previous)
stopifnot(identical(
  dispatched_serial$parallel_execution_policy,
  "serial_original_persistent_engine"
))

cat("CORDA2 stage-parallel serial-equivalence check passed\n")
