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
  if (is.null(reactions) || anyDuplicated(reactions)) stop("invalid reactions")
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
    L = as.numeric(obj),
    lower = as.numeric(lb),
    upper = as.numeric(ub),
    A = A,
    lhs = as.numeric(lhs),
    rhs = as.numeric(rhs),
    maximum = FALSE,
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
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) {
    return(lapply(X, FUN, ...))
  }
  BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)
}
.TEST_BPPARAM <- FALSE
.rc_layer2_task_bpparam <- function() .TEST_BPPARAM

source("R/layer2_corda_evidence.R")
source("R/layer2_corda_lp.R")
source("R/layer2_corda_model.R")

metabolites <- c("A", "B", "X", "Y")
reactions <- c(
  "SRC_A", "NC1", "HC1", "NC_SHARED",
  "MC1", "MC2", "SRC_Y", "MC3"
)
S <- Matrix::Matrix(
  0, nrow = length(metabolites), ncol = length(reactions), sparse = TRUE,
  dimnames = list(metabolites, reactions)
)
S["A", "SRC_A"] <- 1
S["A", "NC1"] <- -1
S["B", "NC1"] <- 1
S["B", "HC1"] <- -1
S["X", "NC_SHARED"] <- 1
S["X", "MC1"] <- -1
S["X", "MC2"] <- -1
S["Y", "SRC_Y"] <- 1
S["Y", "MC3"] <- -1
gem <- list(
  S = S,
  lb = stats::setNames(rep(0, length(reactions)), reactions),
  ub = stats::setNames(rep(100, length(reactions)), reactions)
)
split <- .rc_corda_split_model(gem, tolerance = 1e-8)
classes <- list(
  hc = "HC1",
  mc_module = c("MC1", "MC2", "MC3"),
  mc_evidence = character(),
  mc = c("MC1", "MC2", "MC3"),
  nc = c("NC1", "NC_SHARED"),
  ot = c("SRC_A", "SRC_Y")
)
classes$confidence <- stats::setNames(rep("OT", length(reactions)), reactions)
classes$confidence[classes$nc] <- "NC"
classes$confidence[classes$mc] <- "MC_module"
classes$confidence[classes$hc] <- "HC"
classes$initial_confidence <- classes$confidence
options <- .rc_layer2_corda_options(list(
  model_completion = "corda",
  corda_gamma = 1e5,
  corda_kappa = 1e-2,
  corda_epsilon = 1,
  corda_n = 3L,
  corda_p = 2L,
  corda_seed = 19L
))

.TEST_BPPARAM <- FALSE
serial <- .rc_corda_build_three_stage(
  split, classes, options, solver = "highs", time_limit = 60
)
expected <- reactions
stopifnot(
  setequal(serial$included, expected),
  identical(serial$stage1_associated, "NC1"),
  identical(serial$stage2_promoted_nc, "NC_SHARED"),
  setequal(serial$stage2_promoted_mc, c("MC1", "MC2", "MC3")),
  setequal(serial$stage3_associated_ot, c("SRC_A", "SRC_Y")),
  identical(unname(serial$stage2_nc_support_count[["NC_SHARED"]]), 2L)
)

options_p3 <- options
options_p3$p <- 3L
p3 <- .rc_corda_build_three_stage(
  split, classes, options_p3, solver = "highs", time_limit = 60
)
stopifnot(
  !"NC_SHARED" %in% p3$included,
  !"MC1" %in% p3$included,
  !"MC2" %in% p3$included,
  "MC3" %in% p3$included
)

.TEST_BPPARAM <- BiocParallel::SnowParam(
  workers = 2L, type = "SOCK", progressbar = FALSE
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

persistent_available <- .rc_corda_highs_api_available()
if (persistent_available) {
  original_api <- .rc_corda_highs_api_available
  .rc_corda_highs_api_available <- function() FALSE
  one_shot <- .rc_corda_build_three_stage(
    split, classes, options, solver = "highs", time_limit = 60
  )
  .rc_corda_highs_api_available <- original_api
  stopifnot(setequal(serial$included, one_shot$included))
}

cat("CORDA synthetic checks passed\n")
