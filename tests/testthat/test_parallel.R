test_that("rc_parallel_lapply supports forced sequential execution", {
  out <- rc_parallel_lapply(1:3, function(x) x + 1L, BPPARAM = FALSE)
  expect_equal(out, list(2L, 3L, 4L))
})

test_that("automatic backend selection is platform aware", {
  expect_identical(
    RegCompassR:::.rc_resolve_parallel_backend("auto", "windows"),
    "snow"
  )
  expect_identical(
    RegCompassR:::.rc_resolve_parallel_backend("auto", "unix"),
    "multicore"
  )
  expect_error(
    RegCompassR:::.rc_resolve_parallel_backend("multicore", "windows"),
    "not supported on Windows"
  )
})

test_that("parallel config reserves two detected CPUs", {
  old_slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)
  old_nslots <- Sys.getenv("NSLOTS", unset = NA_character_)
  on.exit({
    if (is.na(old_slurm)) Sys.unsetenv("SLURM_CPUS_PER_TASK") else
      Sys.setenv(SLURM_CPUS_PER_TASK = old_slurm)
    if (is.na(old_nslots)) Sys.unsetenv("NSLOTS") else
      Sys.setenv(NSLOTS = old_nslots)
  }, add = TRUE)
  Sys.setenv(SLURM_CPUS_PER_TASK = "32")
  Sys.unsetenv("NSLOTS")

  config <- rc_parallel_config(workers = 60L)
  expect_equal(config$detected_cpu_capacity, 32L)
  expect_equal(config$reserved_cpus, 2L)
  expect_equal(config$available_workers, 30L)
  expect_equal(config$worker_limit, 30L)
  expect_equal(config$requested_workers, 60L)
})

test_that("default requested worker cap is 10", {
  expect_identical(formals(rc_parallel_config)$workers, 10L)
  expect_identical(formals(rc_default_bpparam)$workers, 10L)
})

test_that("explicit serial backend remains available", {
  expect_null(rc_default_bpparam(workers = 4L, backend = "serial"))
})

test_that("task dispatch shrinks below the worker cap", {
  skip_if_not_installed("BiocParallel")
  old_slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)
  old_nslots <- Sys.getenv("NSLOTS", unset = NA_character_)
  on.exit({
    if (is.na(old_slurm)) Sys.unsetenv("SLURM_CPUS_PER_TASK") else
      Sys.setenv(SLURM_CPUS_PER_TASK = old_slurm)
    if (is.na(old_nslots)) Sys.unsetenv("NSLOTS") else
      Sys.setenv(NSLOTS = old_nslots)
  }, add = TRUE)
  Sys.setenv(SLURM_CPUS_PER_TASK = "64")
  Sys.unsetenv("NSLOTS")

  template <- rc_default_bpparam(workers = 10L, backend = "snow")
  tuned <- RegCompassR:::.rc_parallel_param_for_tasks(template, 3L)
  expect_equal(BiocParallel::bpnworkers(tuned), 3L)
  expect_equal(attr(tuned, "regcompass_worker_limit"), 10L)
})
