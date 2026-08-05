suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
  library(BiocParallel)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
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
    L = as.numeric(obj),
    lower = as.numeric(lb),
    upper = as.numeric(ub),
    A = A,
    lhs = as.numeric(lhs),
    rhs = as.numeric(rhs),
    maximum = FALSE,
    control = highs::highs_control(
      log_to_console = FALSE,
      threads = 1L,
      solver = "simplex"
    )
  )
  list(
    status = .rc_lp_status(answer$status_message, answer$status),
    solution = as.numeric(answer$primal_solution),
    objective = as.numeric(answer$objective_value),
    solver_message = answer$status_message
  )
}
rc_parallel_config <- function(...) list(workers = 4L)
.rc_layer2_task_bpparam <- function() FALSE
rc_regcompass_step_layer2 <- function(...) NULL
.rc_build_celltype_medium_union_gem_cache <- function(...) NULL

source("R/layer2_corda_lp.R")
source("R/layer2_corda_runtime.R")

S <- Matrix::Matrix(
  matrix(
    c(
      1, 0,
      -1, 1,
      0, -1
    ),
    nrow = 2,
    dimnames = list(c("A", "B"), c("SRC", "R", "SINK"))
  ),
  sparse = TRUE
)
gem <- list(
  S = S,
  lb = c(SRC = 0, R = 0, SINK = 0),
  ub = c(SRC = 1000, R = 1000, SINK = 1000)
)
split <- .rc_corda_split_model(gem, tolerance = 1e-7)
engine <- .rc_corda_new_lp_engine(split, solver = "highs", time_limit = 60)
stopifnot(identical(engine$type, "highs_persistent_cpp"))

objective1 <- rep(0, ncol(split$S))
names(objective1) <- colnames(split$S)
objective1[["SRC::forward"]] <- 1
lower1 <- split$lb
upper1 <- split$ub
lower1[["R::forward"]] <- 1

first <- .rc_corda_engine_solve(
  engine,
  objective = objective1,
  lower = lower1,
  upper = upper1
)
engine <- first$engine
stopifnot(identical(first$answer$status, "optimal"))

objective2 <- rep(0, ncol(split$S))
names(objective2) <- colnames(split$S)
objective2[["SINK::forward"]] <- 1
lower2 <- split$lb
upper2 <- split$ub
lower2[["SINK::forward"]] <- 2

second <- .rc_corda_engine_solve(
  engine,
  objective = objective2,
  lower = lower2,
  upper = upper2
)
engine <- second$engine
oracle <- .rc_corda_one_shot_solve(
  engine,
  objective = objective2,
  lower = lower2,
  upper = upper2
)
stopifnot(
  identical(second$answer$status, "optimal"),
  identical(oracle$status, "optimal"),
  isTRUE(all.equal(second$answer$objective, oracle$objective,
                   tolerance = 1e-9)),
  isTRUE(all.equal(second$answer$solution, oracle$solution,
                   tolerance = 1e-8))
)

metrics <- .rc_corda_execution_metrics(engine)
stopifnot(
  metrics$n_solves == 2L,
  metrics$n_objective_coeff_updates == 3L,
  metrics$n_bound_index_updates == 3L,
  metrics$n_transmitted_numeric_values <
    metrics$n_full_vector_numeric_values,
  metrics$transmitted_fraction_of_full < 0.5,
  metrics$n_full_vector_numeric_values_avoided > 0
)

engine <- .rc_corda_release_lp_engine(engine)
stopifnot(is.null(engine$pointer), isTRUE(engine$released))

param <- BiocParallel::SnowParam(
  workers = 4L,
  type = "SOCK",
  tasks = 0L,
  progressbar = FALSE
)
tuned <- .rc_corda_tune_task_bpparam(param, 3L)
stopifnot(
  BiocParallel::bptasks(tuned) == 3L,
  .rc_corda_should_outer_parallel(2L, 8L),
  !.rc_corda_should_outer_parallel(1L, 8L),
  !.rc_corda_should_outer_parallel(8L, 1L)
)

cat(
  "CORDA2 sparse persistent runtime check passed; transmitted fraction = ",
  format(metrics$transmitted_fraction_of_full, digits = 4),
  "\n",
  sep = ""
)
