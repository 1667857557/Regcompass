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
rc_parallel_config <- function(workers = NULL, backend = "auto") {
  list(workers = if (is.null(workers)) 2L else as.integer(workers))
}
.rc_with_internal_single_thread <- function(FUN) FUN()
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  extra <- list(...)
  worker_fun <- function(x) do.call(FUN, c(list(x), extra))
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) {
    return(lapply(X, worker_fun))
  }
  BiocParallel::bplapply(X, worker_fun, BPPARAM = BPPARAM)
}
.TEST_BPPARAM <- FALSE
.rc_layer2_task_bpparam <- function() .TEST_BPPARAM

source("R/layer2_corda_evidence.R")
source("R/layer2_corda_lp.R")
source("R/layer2_corda_paper_contract.R")
source("R/layer2_corda_direction_contract.R")
source("R/layer2_corda_model.R")
source("R/layer2_corda_output_contract.R")
source("R/layer2_corda_target_contract.R")
source("R/layer2_corda_parent_contract.R")
source("R/layer2_corda2_algorithm.R")
source("R/layer2_corda2_options_contract.R")

metabolites <- c("A", "B", "C", "D", "E", "F")
reactions <- c(
  "SRC_A", "M1", "SRC_B", "M2", "H1",
  "SRC_D", "N1", "M3", "M4",
  "SRC_F", "M5", "N2"
)
S <- Matrix::Matrix(
  0, nrow = length(metabolites), ncol = length(reactions), sparse = TRUE,
  dimnames = list(metabolites, reactions)
)
S["A", "SRC_A"] <- 1
S["A", "M1"] <- -1
S["C", "M1"] <- 1
S["B", "SRC_B"] <- 1
S["B", "M2"] <- -1
S["C", "M2"] <- 1
S["C", "H1"] <- -1
S["D", "SRC_D"] <- 1
S["D", "N1"] <- -1
S["E", "N1"] <- 1
S["E", "M3"] <- -1
S["E", "M4"] <- -1
S["F", "SRC_F"] <- 1
S["F", "M5"] <- -1
S["A", "N2"] <- -1

gem <- list(
  S = S,
  lb = stats::setNames(rep(0, length(reactions)), reactions),
  ub = stats::setNames(rep(100, length(reactions)), reactions)
)
split <- .rc_corda_split_model(gem, tolerance = 1e-8)
classes <- list(
  hc = "H1",
  mc_module = c("M1", "M2", "M3", "M4", "M5"),
  mc_evidence = character(),
  mc = c("M1", "M2", "M3", "M4", "M5"),
  nc = c("N1", "N2"),
  ot = c("SRC_A", "SRC_B", "SRC_D", "SRC_F")
)
classes$confidence <- stats::setNames(rep("OT", length(reactions)), reactions)
classes$confidence[classes$nc] <- "NC"
classes$confidence[classes$mc_module] <- "MC_module"
classes$confidence[classes$hc] <- "HC"
classes$initial_confidence <- classes$confidence

options <- .rc_layer2_corda_options(list(
  model_completion = "corda2",
  corda2_redundancies = 3L,
  corda2_penalty_factor = 100,
  corda2_support = 2L,
  corda2_cost_increase = 1.01,
  corda2_target_flux = 1,
  corda2_flux_tolerance = 1e-8
))
stopifnot(
  identical(options$model_completion, "corda2"),
  identical(options$requested_model_completion, "corda2"),
  identical(options$algorithm,
            "resendislab_python_CORDA2_corrected_redundant_path_assessment"),
  identical(options$redundancies, 3L),
  identical(options$support, 2L),
  identical(options$penalty_factor, 100),
  identical(options$cost_increase, 1.01),
  identical(options$random_noise, FALSE)
)

.TEST_BPPARAM <- FALSE
serial <- .rc_corda_build_three_stage(
  split, classes, options, solver = "highs", time_limit = 60
)
stopifnot(
  identical(serial$algorithm,
            "resendislab_python_CORDA2_corrected_redundant_path_assessment"),
  all(c("H1", "M1", "M2", "N1", "M3", "M4", "M5") %in%
        serial$included),
  !"N2" %in% serial$included,
  all(c("M1", "M2") %in% serial$stage1_associated),
  "N1" %in% serial$stage2_promoted_nc,
  all(c("M3", "M4", "M5") %in% serial$stage2_promoted_mc),
  all(c("SRC_A", "SRC_B", "SRC_D", "SRC_F") %in%
        serial$stage3_associated_ot),
  max(serial$redundancies, na.rm = TRUE) >= 2L
)

options_support3 <- options
options_support3$support <- 3L
support3 <- .rc_corda_build_three_stage(
  split, classes, options_support3, solver = "highs", time_limit = 60
)
stopifnot(
  !"N1" %in% support3$stage2_promoted_nc,
  !"M3" %in% support3$included,
  !"M4" %in% support3$included,
  "M5" %in% support3$included
)

if (.Platform$OS.type != "windows") {
  .TEST_BPPARAM <- BiocParallel::MulticoreParam(
    workers = 2L, progressbar = FALSE
  )
  parallel <- .rc_corda_build_three_stage(
    split, classes, options, solver = "highs", time_limit = 60
  )
  BiocParallel::bpstop(.TEST_BPPARAM)
  .TEST_BPPARAM <- FALSE
  stopifnot(
    setequal(serial$included, parallel$included),
    identical(serial$stage2_nc_support_count,
              parallel$stage2_nc_support_count),
    setequal(serial$stage3_associated_ot, parallel$stage3_associated_ot)
  )
}

if (.rc_corda_highs_api_available()) {
  stopifnot(any(vapply(serial$execution, function(x) {
    isTRUE(x$persistent_solver)
  }, logical(1))))
  original_api <- .rc_corda_highs_api_available
  .rc_corda_highs_api_available <- function() FALSE
  one_shot <- .rc_corda_build_three_stage(
    split, classes, options, solver = "highs", time_limit = 60
  )
  .rc_corda_highs_api_available <- original_api
  stopifnot(setequal(serial$included, one_shot$included))
}

cat("Corrected Python CORDA2 synthetic checks passed\n")
