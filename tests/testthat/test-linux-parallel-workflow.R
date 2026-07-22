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

test_that("one-shot upstream workers are propagated to meta-module completion", {
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

  expect_match(text, "local_fastcore_args$workers <- upstream_workers", fixed = TRUE)
  expect_match(text, "local_fastcore_args$backend <- parallel_backend", fixed = TRUE)
  expect_match(text, "result$params$upstream_workers <- upstream_workers", fixed = TRUE)
  expect_match(text, "result$params$layer2_workers <- layer2_workers", fixed = TRUE)
})

test_that("explicit Linux multicore backend creates a MulticoreParam", {
  skip_if(.Platform$OS.type == "windows")
  skip_if_not_installed("BiocParallel")
  param <- rc_default_bpparam(workers = 2L, backend = "multicore")
  expect_true(methods::is(param, "MulticoreParam"))
  expect_equal(BiocParallel::bpnworkers(param), 2L)
})
