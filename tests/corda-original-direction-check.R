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
  if (grepl("unbounded", text)) return("unbounded")
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
  control <- list(log_to_console = FALSE, threads = 1L, solver = "simplex")
  if (is.finite(time_limit)) control$time_limit <- time_limit
  answer <- highs::highs_solve(
    L = as.numeric(obj), lower = as.numeric(lb), upper = as.numeric(ub),
    A = A, lhs = as.numeric(lhs), rhs = as.numeric(rhs), maximum = FALSE,
    control = do.call(highs::highs_control, control)
  )
  list(
    status = .rc_lp_status(answer$status_message, answer$status),
    solution = as.numeric(answer$primal_solution),
    objective = as.numeric(answer$objective_value),
    solver_message = answer$status_message
  )
}
rc_parallel_config <- function(workers = NULL, backend = "auto") list(workers = 1L)
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) lapply(X, FUN, ...)
.rc_layer2_task_bpparam <- function() FALSE
.rc_subset_gem <- function(gem, reactions) gem
rc_prepare_directional_targets <- function(...) data.frame()
.rc_directional_feasibility <- function(...) data.frame()
rc_build_full_gem <- function(gem, ...) gem

source("R/layer2_corda_evidence.R")
source("R/layer2_corda_lp.R")
source("R/layer2_corda_paper_contract.R")
source("R/layer2_corda_direction_contract.R")
source("R/layer2_corda_model.R")
source("R/layer2_corda_output_contract.R")
source("R/layer2_corda_target_contract.R")
source("R/layer2_corda_parent_contract.R")
source("R/layer2_corda2_algorithm.R")

S <- Matrix::Matrix(
  matrix(c(-1, 1), nrow = 2), sparse = TRUE,
  dimnames = list(c("A", "B"), "REV")
)
gem <- list(S = S, lb = c(REV = -10), ub = c(REV = 10))
split <- .rc_corda_split_model(gem, tolerance = 1e-8)

raw_lower <- split$lb
raw_upper <- split$ub
raw_lower[["REV::forward"]] <- 1
raw <- rc_solve_lp(
  obj = c(0, 0), A = split$S,
  lhs = rep(0, nrow(split$S)), rhs = rep(0, nrow(split$S)),
  lb = raw_lower, ub = raw_upper,
  solver = "highs", time_limit = 60
)
stopifnot(identical(raw$status, "optimal"))

bounds <- .rc_corda_target_bounds(
  split, "REV::forward", epsilon = 1
)
closed <- rc_solve_lp(
  obj = c(0, 0), A = split$S,
  lhs = rep(0, nrow(split$S)), rhs = rep(0, nrow(split$S)),
  lb = bounds$lower, ub = bounds$upper,
  solver = "highs", time_limit = 60
)
stopifnot(
  identical(bounds$opposite_variables, "REV::reverse"),
  identical(closed$status, "infeasible")
)

options <- .rc_layer2_corda_options(list(
  model_completion = "corda2",
  corda2_target_flux = 1,
  corda2_penalty_factor = 100,
  corda2_cost_increase = 1.01,
  corda2_redundancies = 3L,
  corda2_support = 5L
))
engine <- .rc_corda_new_lp_engine(split, solver = "highs", time_limit = 60)
assessment <- .rc_corda2_associated_target(
  engine = engine,
  target = "REV::forward",
  directional_confidence = c(
    "REV::forward" = 3L,
    "REV::reverse" = 3L
  ),
  options = options,
  penalize_medium = TRUE,
  redundancies = TRUE,
  stage = "corda2_stage1_high_associations"
)
stopifnot(
  assessment$result$status %in% c("infeasible", "error"),
  identical(
    assessment$result$opposite_direction_blocked,
    "REV::reverse"
  )
)

cat("Corrected Python CORDA2 signed-direction check passed\n")
