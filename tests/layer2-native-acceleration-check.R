`%||%` <- function(x, y) if (is.null(x)) y else x

library(Matrix)
library(Rcpp)

Rcpp::sourceCpp("src/layer2_native.cpp")
.rc_corda2_scan_flux_cpp <- rc_corda2_scan_flux_cpp
options(warn = 2)

source("R/00_utils.R", local = FALSE)

rc_align_bound <- function(x, reactions, default, name) {
  if (is.null(x)) {
    value <- rep(default, length(reactions))
    names(value) <- reactions
    return(value)
  }
  if (!is.null(names(x))) {
    missing <- setdiff(reactions, names(x))
    if (length(missing)) {
      stop(name, " is missing reactions: ", paste(missing, collapse = ", "))
    }
    x <- x[reactions]
  }
  x <- as.numeric(x)
  if (length(x) != length(reactions) || anyNA(x)) {
    stop(name, " does not align with reactions")
  }
  names(x) <- reactions
  x
}

source("R/lp_solver.R", local = FALSE)
source("R/microcompass_vmax_cache.R", local = FALSE)

flux <- c(0, 2, 3, 1e-9, 4, NA_real_)
class_code <- c(1L, 2L, 3L, 4L, 2L, 3L)
track_code <- c(2L, 3L)
scan <- .rc_corda2_scan_flux_cpp(flux, class_code, track_code, 1e-8)
expected_active <- which(is.finite(flux) & flux > 1e-8)
expected_used <- expected_active[class_code[expected_active] %in% track_code]
stopifnot(
  identical(as.integer(scan$active), as.integer(expected_active)),
  identical(as.integer(scan$used), as.integer(expected_used))
)

row_ids <- paste0("row_", seq_len(12L))
model_keys <- stats::setNames(
  rep("/tmp/one-shared-model.rds", length(row_ids)), row_ids
)
tasks <- .rc_microcompass_vmax_tasks(model_keys, workers = 4L)
observed <- unlist(lapply(tasks, `[[`, "row_ids"), use.names = FALSE)
stopifnot(
  length(tasks) == 4L,
  setequal(observed, row_ids),
  length(observed) == length(unique(observed))
)

S <- Matrix::Matrix(
  matrix(
    c(1, -1), nrow = 1,
    dimnames = list("M", c("UP", "TARGET"))
  ),
  sparse = TRUE
)
lb <- c(UP = 0, TARGET = 0)
ub <- c(UP = 10, TARGET = 10)
vmax <- rc_compass_vmax_directional(
  S, lb, ub, "TARGET", direction = "forward", solver = "highs"
)
stopifnot(isTRUE(vmax$feasible))
prepared <- .rc_compass_step2_prepare(
  S, lb, ub, "TARGET", vmax,
  target_direction = "forward", omega = 0.95
)
engine <- .rc_compass_step2_new_engine(prepared$template, "highs")
on.exit(.rc_compass_step2_release_engine(engine), add = TRUE)
if (!identical(engine$type, "highs_persistent_cpp")) {
  stop(
    "Persistent HiGHS Step 2 API was not activated: ",
    engine$persistent_message %||% "unknown reason"
  )
}

penalty_sets <- list(
  c(UP = 0.25, TARGET = 0.50),
  c(UP = 0.75, TARGET = 0.10),
  c(UP = 0.05, TARGET = 1.25)
)
for (penalties in penalty_sets) {
  solved <- .rc_compass_step2_engine_solve(engine, penalties)
  engine <- solved$engine
  persistent <- .rc_compass_step2_result(prepared$template, solved$answer)
  reference <- rc_compass_two_step_lp_directional(
    S, lb, ub, "TARGET", penalties,
    target_direction = "forward", omega = 0.95, solver = "highs"
  )
  stopifnot(
    identical(persistent$feasible, reference$feasible),
    isTRUE(all.equal(persistent$vmax, reference$vmax, tolerance = 1e-10)),
    isTRUE(all.equal(persistent$penalty, reference$penalty, tolerance = 1e-9))
  )
}
metrics <- .rc_compass_step2_engine_metrics(engine)
stopifnot(
  identical(metrics$engine, "highs_persistent_cpp"),
  metrics$n_solves == length(penalty_sets),
  metrics$n_objective_updates > 0L,
  metrics$n_fallback == 0L
)

cat("Layer 2 native acceleration regression passed.\n")
