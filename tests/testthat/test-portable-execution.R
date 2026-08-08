test_that("auto parallel backend follows the operating system", {
  expect_identical(.rc_resolve_parallel_backend("auto", "windows"), "snow")
  expect_identical(.rc_resolve_parallel_backend("auto", "unix"), "multicore")
  expect_identical(.rc_resolve_parallel_backend("snow", "unix"), "snow")
  expect_error(
    .rc_resolve_parallel_backend("multicore", "windows"),
    "not supported on Windows"
  )
})

test_that("parallel configuration records requested and protected execution", {
  one <- rc_parallel_config(workers = 1L, backend = "auto")
  expect_identical(one$actual_backend, "serial")
  expect_identical(one$workers, 1L)
  expect_identical(one$os_type, .Platform$OS.type)
  expect_true(one$detected_cpu_capacity >= 1L)
  expect_true(one$reserved_cpus %in% 0:2)
  expect_identical(one$available_workers,
                   max(1L, one$detected_cpu_capacity - 2L))

  expect_identical(.rc_resolve_parallel_backend("auto", "windows"), "snow")
  expect_identical(.rc_resolve_parallel_backend("auto", "unix"), "multicore")
})

test_that("canonical workflow exposes one adjustable worker cap", {
  args <- formals(rc_run_regcompass)
  expect_identical(args$workers, 10L)
  expect_false("upstream_workers" %in% names(args))
  expect_false("layer2_workers" %in% names(args))
  expect_false("parallel_backend" %in% names(args))
  expect_lt(match("species", names(args)), match("progress", names(args)))

  expect_identical(formals(rc_regcompass_step_grn)$workers, 10L)
  expect_identical(formals(rc_regcompass_step_layer1)$workers, 10L)
  expect_identical(formals(rc_regcompass_step_layer2)$workers, 10L)
  expect_identical(formals(rc_run_regcompass_one_shot)$workers, 10L)
  expect_false("BPPARAM" %in% names(formals(rc_regcompass_step_grn)))
  expect_false("BPPARAM" %in% names(formals(rc_regcompass_step_layer1)))
  expect_false("BPPARAM" %in% names(formals(rc_regcompass_step_layer2)))
  expect_false("parallel" %in% names(formals(rc_regcompass_step_grn)))
  expect_false("parallel" %in% names(formals(rc_regcompass_step_layer1)))
  expect_false("parallel" %in% names(formals(rc_regcompass_step_layer2)))
  expect_error(.rc_stage_worker_config(0L), "at least 1")
})

test_that("internal task thread settings are forced to one and restored", {
  variables <- names(.rc_internal_thread_env())
  before <- Sys.getenv(variables, unset = NA_character_)
  before_mc <- getOption("mc.cores")
  before_internal <- getOption("RegCompassR.internal_workers")

  state <- .rc_set_internal_single_thread()
  on.exit(.rc_restore_internal_threads(state), add = TRUE)

  expect_true(all(Sys.getenv(variables) %in% c("1", "FALSE")))
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "1")
  expect_identical(Sys.getenv("OPENBLAS_NUM_THREADS"), "1")
  expect_identical(Sys.getenv("MKL_NUM_THREADS"), "1")
  expect_identical(Sys.getenv("HIGHS_THREADS"), "1")
  expect_identical(getOption("mc.cores"), 1L)
  expect_identical(getOption("RegCompassR.internal_workers"), 1L)

  .rc_restore_internal_threads(state)
  state <- NULL
  expect_identical(Sys.getenv(variables, unset = NA_character_), before)
  expect_identical(getOption("mc.cores"), before_mc)
  expect_identical(
    getOption("RegCompassR.internal_workers"),
    before_internal
  )
})

test_that("dynamic task sizing does not start more workers than tasks", {
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
  tuned <- .rc_parallel_param_for_tasks(template, 3L)
  expect_equal(BiocParallel::bpnworkers(tuned), 3L)
  expect_equal(attr(tuned, "regcompass_worker_limit"), 10L)
  expect_false(BiocParallel::bpisup(tuned))
})

test_that("package-managed worker pool is stopped after dispatch", {
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

  template <- rc_default_bpparam(workers = 4L, backend = "snow")
  expect_false(BiocParallel::bpisup(template))
  result <- rc_parallel_lapply(1:2, function(x) x + 1L, BPPARAM = template)
  expect_equal(result, list(2L, 3L))
  expect_false(BiocParallel::bpisup(template))
})

test_that("bundled GEM manifest and files are complete", {
  manifest <- rc_bundled_gem_manifest()
  expect_setequal(manifest$species, c("human", "mouse"))
  expect_setequal(manifest$version, c("2.0.0", "1.8.0"))
  expect_true(all(manifest$size_bytes > 0))
  expect_true(all(nzchar(manifest$md5)))

  paths <- file.path(
    system.file("extdata", "gem", package = "RegCompassR"),
    manifest$file
  )
  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), manifest$md5)
})

test_that("bundled human and mouse GEMs persist to requested cache files", {
  human_path <- tempfile(fileext = ".rds")
  mouse_path <- tempfile(fileext = ".rds")
  human <- rc_prepare_gem(
    species = "human", version = "2.0.0", source = "bundled",
    save_rds = human_path
  )
  mouse <- rc_prepare_gem(
    species = "mouse", version = "1.8.0", source = "bundled",
    save_rds = mouse_path
  )
  expect_true(file.exists(human_path))
  expect_true(file.exists(mouse_path))
  expect_silent(rc_validate_species_gem(human, "human"))
  expect_silent(rc_validate_species_gem(mouse, "mouse"))
  expect_identical(human$model_info$distribution, "bundled_with_RegCompassR")
  expect_identical(mouse$model_info$distribution, "bundled_with_RegCompassR")
  expect_identical(readRDS(human_path)$model_info$species, "human")
  expect_identical(readRDS(mouse_path)$model_info$species, "mouse")
})

test_that("step monitor writes timing and can suppress progress", {
  outdir <- tempfile("regcompass-timing-")
  monitor <- .rc_step_monitor_start(
    "unit_test", outdir = outdir, progress = FALSE
  )
  value <- .rc_step_monitor_finish(list(ok = TRUE), monitor)
  expect_true(is.data.frame(value$timing))
  expect_identical(value$timing$stage, "unit_test")
  expect_true(value$timing$elapsed_seconds >= 0)
  expect_true(file.exists(file.path(outdir, "step_timing.tsv")))
})

test_that("known stages report success only after the final RDS is committed", {
  outdir <- tempfile("regcompass-stage-commit-")
  monitor <- .rc_step_monitor_start("grn", outdir, progress = FALSE)
  value <- .rc_step_monitor_finish(list(ok = TRUE), monitor)
  expect_false(file.exists(file.path(outdir, "step_timing.tsv")))
  saveRDS(value, file.path(outdir, "step_grn.rds"))
  .rc_step_monitor_fail(monitor)
  timing <- utils::read.delim(
    file.path(outdir, "step_timing.tsv"), stringsAsFactors = FALSE
  )
  expect_identical(timing$status, "success")

  failed_outdir <- tempfile("regcompass-stage-fail-")
  failed <- .rc_step_monitor_start("layer1", failed_outdir, progress = FALSE)
  .rc_step_monitor_finish(list(ok = TRUE), failed)
  .rc_step_monitor_fail(failed)
  failed_timing <- utils::read.delim(
    file.path(failed_outdir, "step_timing.tsv"), stringsAsFactors = FALSE
  )
  expect_identical(failed_timing$status, "error")
})

test_that("every public workflow stage exposes progress control", {
  functions <- list(
    rc_regcompass_step_grn,
    rc_regcompass_step_metacells,
    rc_regcompass_step_meta_modules,
    rc_regcompass_step_layer1,
    rc_regcompass_step_layer2,
    rc_regcompass_step_results,
    rc_regcompass_step_target_union,
    rc_run_regcompass,
    rc_run_regcompass_one_shot
  )
  expect_true(all(vapply(
    functions,
    function(fun) "progress" %in% names(formals(fun)),
    logical(1)
  )))
})
