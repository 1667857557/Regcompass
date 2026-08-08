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
  expect_false(grepl(".rc_complete_celltype_medium_union_gem(", text, fixed = TRUE))
  expect_false(grepl(".rc_fastcore_", text, fixed = TRUE))
})

test_that("canonical workflow exposes one dynamic worker cap", {
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
  expect_match(text, "workers = 10L", fixed = TRUE)
  expect_match(text, ".rc_stage_worker_config(workers", fixed = TRUE)
  expect_match(text, "workers = worker_limit", fixed = TRUE)
  expect_match(text, "layer2_args = layer2_args", fixed = TRUE)
  expect_match(text, "result$params$workers", fixed = TRUE)
  expect_match(text, "result$params$parallel_backend", fixed = TRUE)
  expect_false(grepl("upstream_workers", text, fixed = TRUE))
  expect_false(grepl("layer2_workers", text, fixed = TRUE))
  expect_false(grepl("parallel_backend =", text, fixed = TRUE))
})

test_that("parallel runtime forces one-thread child environments and cleanup", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  stage_candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "R", "stage_parallel_lifecycle.R") else character(),
    file.path("R", "stage_parallel_lifecycle.R"),
    file.path("..", "R", "stage_parallel_lifecycle.R"),
    file.path("..", "..", "R", "stage_parallel_lifecycle.R")
  ))
  parallel_candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "R", "parallel.R") else character(),
    file.path("R", "parallel.R"),
    file.path("..", "R", "parallel.R"),
    file.path("..", "..", "R", "parallel.R")
  ))
  stage_candidates <- stage_candidates[file.exists(stage_candidates)]
  parallel_candidates <- parallel_candidates[file.exists(parallel_candidates)]
  if (!length(stage_candidates) || !length(parallel_candidates)) {
    skip("parallel runtime sources are unavailable.")
  }
  stage_text <- paste(readLines(stage_candidates[[1L]], warn = FALSE), collapse = "\n")
  parallel_text <- paste(readLines(parallel_candidates[[1L]], warn = FALSE), collapse = "\n")
  expect_match(stage_text, "OMP_NUM_THREADS = \"1\"", fixed = TRUE)
  expect_match(stage_text, "OPENBLAS_NUM_THREADS = \"1\"", fixed = TRUE)
  expect_match(stage_text, "MKL_NUM_THREADS = \"1\"", fixed = TRUE)
  expect_match(stage_text, "HIGHS_THREADS = \"1\"", fixed = TRUE)
  expect_match(stage_text, "mc.cores = 1L", fixed = TRUE)
  expect_match(parallel_text, "BiocParallel::bpstart(BPPARAM)", fixed = TRUE)
  expect_match(parallel_text, ".rc_release_bpparam(BPPARAM)", fixed = TRUE)
  expect_match(parallel_text, "gc(verbose = FALSE, full = TRUE)", fixed = TRUE)
  expect_match(parallel_text, "detected_cpu_capacity - 2L", fixed = TRUE)
})

test_that("explicit Linux multicore backend creates a MulticoreParam", {
  skip_if(.Platform$OS.type == "windows")
  skip_if_not_installed("BiocParallel")
  old_slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)
  old_nslots <- Sys.getenv("NSLOTS", unset = NA_character_)
  on.exit({
    if (is.na(old_slurm)) Sys.unsetenv("SLURM_CPUS_PER_TASK") else
      Sys.setenv(SLURM_CPUS_PER_TASK = old_slurm)
    if (is.na(old_nslots)) Sys.unsetenv("NSLOTS") else
      Sys.setenv(NSLOTS = old_nslots)
  }, add = TRUE)
  Sys.setenv(SLURM_CPUS_PER_TASK = "8")
  Sys.unsetenv("NSLOTS")
  param <- rc_default_bpparam(workers = 2L, backend = "multicore")
  expect_true(methods::is(param, "MulticoreParam"))
  expect_equal(BiocParallel::bpnworkers(param), 2L)
})
