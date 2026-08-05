suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
.rc_bind_frames_fill <- function(values) data.frame()
.rc_lp_status <- function(message = "", code = NA_integer_) {
  text <- tolower(paste(message, collapse = " "))
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("optimal", text)) return("optimal")
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
rc_parallel_config <- function(...) list(workers = 1L)
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) lapply(X, FUN, ...)
.rc_layer2_task_bpparam <- function() FALSE

source("R/layer2_corda_evidence.R")
source("R/layer2_corda_lp.R")
source("R/layer2_corda_output_contract.R")
source("R/layer2_corda_direction_contract.R")
source("R/layer2_corda2_algorithm.R")
source("R/layer2_corda2_algorithm_build.R")
source("R/layer2_corda2_options_contract.R")

S <- Matrix::Matrix(
  matrix(c(-1, 1), nrow = 2), sparse = TRUE,
  dimnames = list(c("A", "B"), "REV")
)
gem <- list(S = S, lb = c(REV = -10), ub = c(REV = 10))
split <- .rc_corda_split_model(
  gem,
  tolerance = .rc_corda2_solver_feasibility_tolerance("highs")
)
stopifnot(
  identical(
    split$direction_table$variable_id,
    c("REV::forward", "REV::reverse")
  ),
  split$ub[["REV::forward"]] == 1e6,
  split$ub[["REV::reverse"]] == 1e6
)

bounds <- .rc_corda_target_bounds(split, "REV::forward", epsilon = 1)
stopifnot(
  identical(bounds$opposite_variables, "REV::reverse"),
  identical(bounds$opposite_direction_blocked, character()),
  bounds$lower[["REV::forward"]] == 1,
  bounds$upper[["REV::reverse"]] == 1e6
)
answer <- rc_solve_lp(
  obj = c(0, 0), A = split$S,
  lhs = rep(0, nrow(split$S)), rhs = rep(0, nrow(split$S)),
  lb = bounds$lower, ub = bounds$upper,
  solver = "highs", time_limit = Inf
)
stopifnot(
  identical(answer$status, "optimal"),
  answer$solution[[1L]] > split$tolerance,
  answer$solution[[2L]] > split$tolerance
)

small <- .rc_corda_split_model(
  list(S = S, lb = c(REV = 0), ub = c(REV = 0.5)),
  tolerance = .rc_corda2_solver_feasibility_tolerance("highs")
)
stopifnot(small$ub[["REV::forward"]] == 1e6)
small$ub[["REV::forward"]] <- 0.5
failed_order <- try(
  .rc_corda_target_bounds(small, "REV::forward", epsilon = 1),
  silent = TRUE
)
stopifnot(inherits(failed_order, "try-error"))

options <- .rc_layer2_corda_options(list(model_completion = "corda2"))
stopifnot(
  identical(options$met_prod, NULL),
  identical(options$n, 3L),
  identical(options$penalty_factor, 100),
  identical(options$support, 5L),
  identical(options$cost_increase, 1.01),
  identical(options$target_flux, 1),
  identical(options$upper_bound, 1e6),
  identical(options$source_semantics, "exact")
)
for (parameter in c(
  "corda2_cost_increase", "corda2_target_flux", "corda2_flux_tolerance"
)) {
  args <- list(model_completion = "corda2")
  args[[parameter]] <- 2
  stopifnot(inherits(try(.rc_layer2_corda_options(args), silent = TRUE),
                     "try-error"))
}

medium_code <- paste(
  deparse(body(.rc_corda2_minimize_medium_targets)), collapse = "\n"
)
stopifnot(
  grepl("objective[[target]] <- 1", medium_code, fixed = TRUE),
  grepl("answer$objective > options$target_flux", medium_code, fixed = TRUE),
  !grepl("objective[[target]] <- -1", medium_code, fixed = TRUE)
)

build_code <- paste(deparse(body(.rc_corda_build_three_stage)), collapse = "\n")
stopifnot(
  grepl("stage1_targets", build_code, fixed = TRUE),
  grepl("stage2_targets", build_code, fixed = TRUE),
  grepl("split_stage3$ub[[variable]] <- 0", build_code, fixed = TRUE),
  grepl("python_serial_mutation_order", build_code, fixed = TRUE)
)

cat("Direct exact Python CORDA2 target-direction semantics passed\n")
