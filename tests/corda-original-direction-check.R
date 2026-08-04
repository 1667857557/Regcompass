suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
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

source("R/layer2_corda_lp.R")
source("R/layer2_corda_paper_contract.R")
source("R/layer2_corda_direction_contract.R")

S <- Matrix::Matrix(
  matrix(c(-1, 1), nrow = 2), sparse = TRUE,
  dimnames = list(c("A", "B"), "REV")
)
gem <- list(S = S, lb = c(REV = -10), ub = c(REV = 10))
split <- .rc_corda_split_model(gem, tolerance = 1e-8)
stopifnot(setequal(
  split$direction_table$variable_id,
  c("REV::forward", "REV::reverse")
))

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
stopifnot(
  identical(bounds$opposite_variables, "REV::reverse"),
  identical(unname(bounds$lower[["REV::reverse"]]), 0),
  identical(unname(bounds$upper[["REV::reverse"]]), 0),
  identical(unname(bounds$lower[["REV::forward"]]), 1)
)
closed <- rc_solve_lp(
  obj = c(0, 0), A = split$S,
  lhs = rep(0, nrow(split$S)), rhs = rep(0, nrow(split$S)),
  lb = bounds$lower, ub = bounds$upper,
  solver = "highs", time_limit = 60
)
stopifnot(identical(closed$status, "infeasible"))

options <- list(
  gamma = 1e5, kappa = 0, epsilon = 1,
  seed = 1L, flux_tolerance = 1e-8
)
task <- data.frame(
  variable_id = "REV::forward",
  reaction_id = "REV",
  direction = "forward",
  stage = "stage1_hc_dependencies",
  replicate = 1L,
  kind = "dependency",
  stringsAsFactors = FALSE
)
engine <- .rc_corda_new_lp_engine(split, solver = "highs", time_limit = 60)
dependency <- .rc_corda_dependency_task(
  engine, task, confidence = c(REV = "HC"), options = options
)
stopifnot(
  dependency$result$status %in% c("infeasible", "error"),
  identical(
    dependency$result$opposite_direction_blocked,
    "REV::reverse"
  )
)
engine <- dependency$engine
feasibility <- .rc_corda_feasibility_task(engine, task, options)
stopifnot(
  feasibility$result$status %in% c("infeasible", "blocked", "error"),
  identical(
    feasibility$result$opposite_direction_blocked,
    "REV::reverse"
  )
)

cat("Original CORDA signed-direction check passed\n")
