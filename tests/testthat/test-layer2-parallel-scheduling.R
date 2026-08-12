test_that("parallel dispatch tuning is isolated and load balanced", {
  skip_if_not_installed("BiocParallel")
  old_slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)
  old_nslots <- Sys.getenv("NSLOTS", unset = NA_character_)
  old_options <- options(RegCompassR.progress = FALSE)
  on.exit({
    options(old_options)
    if (is.na(old_slurm)) Sys.unsetenv("SLURM_CPUS_PER_TASK") else
      Sys.setenv(SLURM_CPUS_PER_TASK = old_slurm)
    if (is.na(old_nslots)) Sys.unsetenv("NSLOTS") else
      Sys.setenv(NSLOTS = old_nslots)
  }, add = TRUE)
  Sys.setenv(SLURM_CPUS_PER_TASK = "16")
  Sys.unsetenv("NSLOTS")

  template <- BiocParallel::SnowParam(
    workers = 8L, type = "SOCK", tasks = 3L, progressbar = FALSE
  )
  attr(template, "regcompass_worker_limit") <- 8L
  tuned <- RegCompassR:::.rc_parallel_param_for_tasks(template, 100L)

  expect_equal(BiocParallel::bpnworkers(template), 8L)
  expect_equal(BiocParallel::bptasks(template), 3L)
  expect_equal(BiocParallel::bpnworkers(tuned), 8L)
  expect_equal(BiocParallel::bptasks(tuned), 64L)
  expect_equal(attr(tuned, "regcompass_dispatch_tasks"), 64L)
  expect_equal(attr(tuned, "regcompass_dynamic_chunks_per_worker"), 8L)
})

test_that("CORDA2 structural task tuning cannot throttle later Layer2 scoring", {
  skip_if_not_installed("BiocParallel")
  old_slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)
  old_nslots <- Sys.getenv("NSLOTS", unset = NA_character_)
  old_options <- options(RegCompassR.progress = FALSE)
  on.exit({
    options(old_options)
    if (is.na(old_slurm)) Sys.unsetenv("SLURM_CPUS_PER_TASK") else
      Sys.setenv(SLURM_CPUS_PER_TASK = old_slurm)
    if (is.na(old_nslots)) Sys.unsetenv("NSLOTS") else
      Sys.setenv(NSLOTS = old_nslots)
  }, add = TRUE)
  Sys.setenv(SLURM_CPUS_PER_TASK = "16")
  Sys.unsetenv("NSLOTS")

  stage_template <- BiocParallel::SnowParam(
    workers = 8L, type = "SOCK", tasks = 0L, progressbar = FALSE
  )
  attr(stage_template, "regcompass_worker_limit") <- 8L

  structural <- RegCompassR:::.rc_corda_tune_task_bpparam(
    stage_template, 3L
  )
  expect_equal(BiocParallel::bptasks(stage_template), 0L)
  expect_equal(BiocParallel::bpnworkers(stage_template), 8L)
  expect_equal(BiocParallel::bptasks(structural), 3L)
  expect_equal(BiocParallel::bpnworkers(structural), 3L)

  scoring <- RegCompassR:::.rc_parallel_param_for_tasks(
    stage_template, 100L
  )
  expect_equal(BiocParallel::bpnworkers(scoring), 8L)
  expect_equal(BiocParallel::bptasks(scoring), 64L)
})

test_that("dispatch task count scales with expensive workload size", {
  expect_equal(
    RegCompassR:::.rc_parallel_dispatch_task_count(3L, 8L),
    3L
  )
  expect_equal(
    RegCompassR:::.rc_parallel_dispatch_task_count(100L, 8L),
    64L
  )
  expect_equal(
    RegCompassR:::.rc_parallel_dispatch_task_count(1000L, 20L),
    160L
  )
})
