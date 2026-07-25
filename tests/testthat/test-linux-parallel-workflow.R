test_that("Stage 3 does not run a local FASTCORE loop", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "R", "workflow_z_union_gem.R") else character(),
    file.path("R", "workflow_z_union_gem.R"),
    file.path("..", "R", "workflow_z_union_gem.R"),
    file.path("..", "..", "R", "workflow_z_union_gem.R")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("workflow_z_union_gem.R is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")

  expect_match(text, ".rc_build_condition_meta_modules <- function", fixed = TRUE)
  expect_match(text, "none_at_meta_module_stage", fixed = TRUE)
  expect_false(grepl(".rc_complete_stratum_meta_modules(", text, fixed = TRUE))
  expect_false(grepl("local_fastcore_by_meta_module", text, fixed = TRUE))
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
  expect_match(text, "pando_infer_args$parallel <- FALSE", fixed = TRUE)
  expect_match(text, "result$params$pando_internal_parallel <- FALSE", fixed = TRUE)
  expect_match(text, "Local FASTCORE was removed", fixed = TRUE)
  expect_match(text, "layer2_args$model_params", fixed = TRUE)
  expect_match(text, "result$params$internal_threads_per_task <- 1L", fixed = TRUE)
  expect_match(text, "parallel_worker_lifecycle", fixed = TRUE)
  expect_false(grepl("local_fastcore_args$workers", text, fixed = TRUE))
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
