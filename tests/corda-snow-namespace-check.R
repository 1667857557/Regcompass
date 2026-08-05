suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
  library(BiocParallel)
})

root <- normalizePath(".", mustWork = TRUE)
work <- tempfile("corda2_exact_snow_package_")
pkg <- file.path(work, "Corda2ExactSnowCheck")
lib <- file.path(work, "library")
dir.create(file.path(pkg, "R"), recursive = TRUE)
dir.create(lib, recursive = TRUE)

files <- c(
  "layer2_corda_evidence.R",
  "layer2_corda_lp.R",
  "layer2_corda_output_contract.R",
  "layer2_corda_direction_contract.R",
  "layer2_corda2_algorithm.R",
  "layer2_corda2_algorithm_build.R",
  "layer2_corda2_options_contract.R"
)

writeLines(c(
  "Package: Corda2ExactSnowCheck",
  "Title: Exact CORDA2 Snow Namespace Check",
  "Version: 0.0.1",
  "Description: Installed namespace test for direct pinned Python CORDA2 semantics.",
  "License: MIT",
  "Encoding: UTF-8",
  "Imports: Matrix, methods, highs, BiocParallel",
  "Collate:",
  "    'support.R'",
  paste0("    '", files, "'"),
  "    'run.R'"
), file.path(pkg, "DESCRIPTION"))
writeLines("export(run_corda2_exact_snow_check)", file.path(pkg, "NAMESPACE"))
for (file in files) {
  stopifnot(file.copy(
    file.path(root, "R", file), file.path(pkg, "R", file),
    overwrite = TRUE
  ))
}

support <- c(
  "`%||%` <- function(x, y) if (is.null(x)) y else x",
  ".rc_as_dgCMatrix <- function(x) methods::as(x, 'dgCMatrix')",
  ".rc_bind_frames_fill <- function(values) {",
  "  values <- values[vapply(values, is.data.frame, logical(1))]",
  "  values <- values[vapply(values, nrow, integer(1)) > 0L]",
  "  if (!length(values)) return(data.frame())",
  "  columns <- unique(unlist(lapply(values, colnames), use.names = FALSE))",
  "  values <- lapply(values, function(value) {",
  "    missing <- setdiff(columns, colnames(value))",
  "    for (name in missing) value[[name]] <- NA",
  "    value[, columns, drop = FALSE]",
  "  })",
  "  answer <- do.call(rbind, values)",
  "  rownames(answer) <- NULL",
  "  answer",
  "}",
  ".rc_lp_status <- function(message = '', code = NA_integer_) {",
  "  text <- tolower(paste(message, collapse = ' '))",
  "  if (grepl('infeasible', text)) return('infeasible')",
  "  if (grepl('optimal', text)) return('optimal')",
  "  'error'",
  "}",
  "rc_validate_gem <- function(gem) {",
  "  S <- .rc_as_dgCMatrix(gem$S)",
  "  reactions <- colnames(S)",
  "  lb <- as.numeric(gem$lb[reactions]); ub <- as.numeric(gem$ub[reactions])",
  "  names(lb) <- names(ub) <- reactions",
  "  list(S = S, lb = lb, ub = ub, reactions = reactions)",
  "}",
  "rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub, solver = 'highs', time_limit = Inf) {",
  "  answer <- highs::highs_solve(",
  "    L = as.numeric(obj), lower = as.numeric(lb), upper = as.numeric(ub),",
  "    A = A, lhs = as.numeric(lhs), rhs = as.numeric(rhs), maximum = FALSE,",
  "    control = highs::highs_control(log_to_console = FALSE, output_flag = FALSE,",
  "      threads = 1L, solver = 'simplex', primal_feasibility_tolerance = 1e-7,",
  "      time_limit = as.numeric(time_limit)))",
  "  list(status = .rc_lp_status(answer$status_message, answer$status),",
  "       solution = as.numeric(answer$primal_solution),",
  "       objective = as.numeric(answer$objective_value),",
  "       solver_message = answer$status_message)",
  "}",
  "rc_parallel_config <- function(...) list(workers = 2L)",
  ".TEST_CONTEXT <- new.env(parent = emptyenv())",
  ".TEST_CONTEXT$bpparam <- FALSE",
  ".rc_layer2_task_bpparam <- function() .TEST_CONTEXT$bpparam",
  "rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {",
  "  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) return(lapply(X, FUN, ...))",
  "  BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)",
  "}"
)
writeLines(support, file.path(pkg, "R", "support.R"))

run <- c(
  "run_corda2_exact_snow_check <- function() {",
  "  S <- Matrix::Matrix(matrix(c(1, -1), nrow = 1), sparse = TRUE,",
  "                      dimnames = list('A', c('SRC', 'H')))",
  "  gem <- list(S = S, lb = c(SRC = 0, H = 0), ub = c(SRC = 1000, H = 1000))",
  "  classes <- list(hc = 'H', mc_module = character(), mc_evidence = character(),",
  "                  mc = character(), nc = character(), ot = 'SRC')",
  "  classes$confidence <- c(SRC = 'OT', H = 'HC')",
  "  classes$initial_confidence <- classes$confidence",
  "  options <- .rc_layer2_corda_options(list(model_completion = 'corda2'))",
  "  options$feasibility_tolerance <- .rc_corda2_solver_feasibility_tolerance('highs')",
  "  split <- .rc_corda_split_model(gem, tolerance = options$feasibility_tolerance)",
  "  .TEST_CONTEXT$bpparam <- BiocParallel::SnowParam(",
  "    workers = 2L, type = 'SOCK', progressbar = FALSE,",
  "    exportglobals = TRUE, exportvariables = TRUE)",
  "  BiocParallel::bpstart(.TEST_CONTEXT$bpparam)",
  "  on.exit({ try(BiocParallel::bpstop(.TEST_CONTEXT$bpparam), silent = TRUE);",
  "            .TEST_CONTEXT$bpparam <- FALSE }, add = TRUE)",
  "  result <- .rc_corda_build_three_stage(",
  "    split, classes, options, solver = 'highs', time_limit = Inf)",
  "  stopifnot(setequal(result$included, c('SRC', 'H')))",
  "  stopifnot(all(vapply(result$execution, function(x) {",
  "    identical(x$workers, 1L) && identical(x$target_parallelism, FALSE)",
  "  }, logical(1))))",
  "  stopifnot(identical(result$stage_update_policy, 'python_serial_mutation_order'))",
  "  TRUE",
  "}"
)
writeLines(run, file.path(pkg, "R", "run.R"))

status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--no-byte-compile", "-l", shQuote(lib), shQuote(pkg))
)
stopifnot(identical(status, 0L))
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(Corda2ExactSnowCheck))
stopifnot(isTRUE(run_corda2_exact_snow_check()))
cat("Direct exact Python CORDA2 Snow namespace check passed\n")
