test_that("Stage 3 does not run FASTCORE", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "R", "meta_module_construction.R") else character(),
    file.path("R", "meta_module_construction.R"),
    file.path("..", "R", "meta_module_construction.R"),
    file.path("..", "..", "R", "meta_module_construction.R")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("meta_module_construction.R is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")
  expect_match(text, ".rc_build_condition_meta_modules <- function", fixed = TRUE)
  expect_match(text, "none_at_meta_module_stage", fixed = TRUE)
  expect_false(grepl(".rc_complete_medium_union_gem(", text, fixed = TRUE))
  expect_false(grepl(".rc_fastcore_", text, fixed = TRUE))
})

test_that("canonical workflow injects stage-scoped workers", {
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
  expect_match(text, "layer2_args = layer2_args", fixed = TRUE)
  expect_match(text, "result$params$upstream_workers", fixed = TRUE)
  expect_match(text, "result$params$layer2_workers", fixed = TRUE)
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
