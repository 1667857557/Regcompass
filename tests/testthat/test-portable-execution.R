test_that("auto parallel backend follows the operating system", {
  expect_identical(.rc_resolve_parallel_backend("auto", "windows"), "snow")
  expect_identical(.rc_resolve_parallel_backend("auto", "unix"), "multicore")
  expect_identical(.rc_resolve_parallel_backend("snow", "unix"), "snow")
  expect_error(
    .rc_resolve_parallel_backend("multicore", "windows"),
    "not supported on Windows"
  )
})

test_that("parallel configuration reserves two CPUs and caps every request", {
  capacity <- .rc_worker_capacity()
  expect_identical(capacity$reserved_cpus, 2L)
  expect_identical(
    capacity$worker_ceiling,
    max(1L, capacity$available_cpus - 2L)
  )
  expect_identical(
    .rc_normalize_worker_budget(capacity$available_cpus + 100L),
    capacity$worker_ceiling
  )

  one <- rc_parallel_config(workers = 1L, backend = "auto")
  expect_identical(one$actual_backend, "serial")
  expect_identical(one$workers, 1L)
  expect_identical(one$os_type, .Platform$OS.type)
  expect_identical(one$reserved_cpus, 2L)
  expect_identical(one$worker_ceiling, capacity$worker_ceiling)

  requested <- rc_parallel_config(workers = 100000L, backend = "auto")
  expect_lte(requested$worker_budget, capacity$worker_ceiling)
  expect_identical(requested$available_cpus, capacity$available_cpus)
})

test_that("canonical workflow exposes one global worker budget", {
  args <- formals(rc_run_regcompass)
  expect_true("workers" %in% names(args))
  expect_false(any(c(
    "upstream_workers", "layer2_workers", "parallel_backend", "BPPARAM",
    "parallel"
  ) %in% names(args)))
  expect_identical(
    eval(args$workers),
    getOption("RegCompassR.workers", 10L)
  )
  expect_lt(match("species", names(args)), match("progress", names(args)))

  parallel_steps <- list(
    grn = rc_regcompass_step_grn,
    layer1 = rc_regcompass_step_layer1,
    layer2 = rc_regcompass_step_layer2
  )
  for (fun in parallel_steps) {
    step_args <- names(formals(fun))
    expect_true("workers" %in% step_args)
    expect_false(any(c("BPPARAM", "parallel") %in% step_args))
  }
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

test_that("serial stage wrapper applies one-thread contract on errors", {
  before <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  expect_error(
    .rc_with_stage_workers(
      1L,
      function(param, config) {
        expect_identical(param, FALSE)
        expect_identical(config$actual_backend, "serial")
        expect_identical(Sys.getenv("OMP_NUM_THREADS"), "1")
        stop("expected-stage-error")
      }
    ),
    "expected-stage-error"
  )
  expect_identical(
    Sys.getenv("OMP_NUM_THREADS", unset = NA_character_),
    before
  )
})

test_that("package-managed worker pool is stopped after its stage", {
  skip_if_not_installed("BiocParallel")
  skip_if(.rc_worker_capacity()$worker_ceiling < 2L)
  param <- .rc_with_stage_workers(
    2L,
    function(param, config) {
      expect_false(identical(param, FALSE))
      expect_true(BiocParallel::bpisup(param))
      expect_identical(config$workers, 2L)
      param
    }
  )
  expect_false(BiocParallel::bpisup(param))
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
