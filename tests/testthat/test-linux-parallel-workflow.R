test_that("local FASTCORE completion is implemented as a parallel module loop", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "R", "workflow_utils.R") else character(),
    file.path("R", "workflow_utils.R"),
    file.path("..", "R", "workflow_utils.R"),
    file.path("..", "..", "R", "workflow_utils.R")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("workflow_utils.R is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")

  expect_match(text, ".rc_complete_stratum_meta_modules <- function", fixed = TRUE)
  expect_match(text, "rc_parallel_lapply(", fixed = TRUE)
  expect_match(text, "local_fastcore_by_meta_module", fixed = TRUE)
  expect_match(text, "parallel_backend", fixed = TRUE)
  expect_match(text, "parallel_workers", fixed = TRUE)
  expect_match(text, "backend = \"auto\"", fixed = TRUE)
})

test_that("canonical workflow injects layered workers and automatic backend", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "R", "regcompass.R") else character(),
    file.path("R", "regcompass.R"),
    file.path("..", "R", "regcompass.R"),
    file.path("..", "..", "R", "regcompass.R")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("regcompass.R is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")

  expect_match(text, "upstream_workers = 6L", fixed = TRUE)
  expect_match(text, "layer2_workers = 30L", fixed = TRUE)
  expect_match(text, ".rc_stage_worker_config(", fixed = TRUE)
  expect_match(text, ".rc_with_stage_workers(", fixed = TRUE)
  expect_match(text, "local_fastcore_args$workers <- upstream_config$workers", fixed = TRUE)
  expect_match(text, "local_fastcore_args$backend <- \"auto\"", fixed = TRUE)
  expect_match(text, "result$params$internal_threads_per_task <- 1L", fixed = TRUE)
  expect_match(text, "parallel_worker_lifecycle", fixed = TRUE)
  expect_false(grepl("parallel_backend =", text, fixed = TRUE))
})

test_that("stage worker lifecycle forces one-thread child environments", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "R", "stage_parallel_lifecycle.R") else character(),
    file.path("R", "stage_parallel_lifecycle.R"),
    file.path("..", "R", "stage_parallel_lifecycle.R"),
    file.path("..", "..", "R", "stage_parallel_lifecycle.R")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("stage_parallel_lifecycle.R is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")

  expect_match(text, "OMP_NUM_THREADS = \"1\"", fixed = TRUE)
  expect_match(text, "OPENBLAS_NUM_THREADS = \"1\"", fixed = TRUE)
  expect_match(text, "MKL_NUM_THREADS = \"1\"", fixed = TRUE)
  expect_match(text, "mc.cores = 1L", fixed = TRUE)
  expect_match(text, "BiocParallel::bpstart(param)", fixed = TRUE)
  expect_match(text, ".rc_release_bpparam(param)", fixed = TRUE)
  expect_match(text, "gc(verbose = FALSE, full = TRUE)", fixed = TRUE)
})

test_that("explicit Linux multicore backend creates a MulticoreParam", {
  skip_if(.Platform$OS.type == "windows")
  skip_if_not_installed("BiocParallel")
  param <- rc_default_bpparam(workers = 2L, backend = "multicore")
  expect_true(methods::is(param, "MulticoreParam"))
  expect_equal(BiocParallel::bpnworkers(param), 2L)
})
