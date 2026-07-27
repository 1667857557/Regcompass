test_that("auto parallel backend follows the operating system", {
  expect_identical(.rc_resolve_parallel_backend("auto", "windows"), "snow")
  expect_identical(.rc_resolve_parallel_backend("auto", "unix"), "multicore")
  expect_identical(.rc_resolve_parallel_backend("snow", "unix"), "snow")
  expect_error(
    .rc_resolve_parallel_backend("multicore", "windows"),
    "not supported on Windows"
  )
})

test_that("parallel configuration records requested and actual execution", {
  one <- rc_parallel_config(workers = 1L, backend = "auto")
  expect_identical(one$actual_backend, "serial")
  expect_identical(one$workers, 1L)
  expect_identical(one$os_type, .Platform$OS.type)

  windows <- .rc_resolve_parallel_backend("auto", "windows")
  linux <- .rc_resolve_parallel_backend("auto", "unix")
  expect_identical(windows, "snow")
  expect_identical(linux, "multicore")
})

test_that("canonical workflow exposes only two layered worker counts", {
  args <- formals(rc_run_regcompass)
  expect_identical(args$upstream_workers, 6L)
  expect_identical(args$layer2_workers, 30L)
  expect_false("parallel_backend" %in% names(args))
  expect_lt(match("species", names(args)), match("progress", names(args)))

  upstream <- .rc_stage_worker_config(1L, "upstream_workers")
  layer2 <- .rc_stage_worker_config(1L, "layer2_workers")
  expect_identical(upstream$actual_backend, "serial")
  expect_identical(layer2$actual_backend, "serial")
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

test_that("step monitor prints timing and does not persist it", {
  outdir <- tempfile("regcompass-timing-")
  monitor <- .rc_step_monitor_start(
    "unit_test", outdir = outdir, progress = FALSE
  )
  value <- expect_message(
    .rc_step_monitor_finish(list(ok = TRUE), monitor),
    "RegCompass timing: unit_test \\[success\\]"
  )
  expect_true(value$ok)
  expect_null(value$timing)
  expect_false(file.exists(file.path(outdir, "step_timing.tsv")))
})

test_that("known stages print final status only after artifact commit", {
  outdir <- tempfile("regcompass-stage-commit-")
  monitor <- .rc_step_monitor_start("grn", outdir, progress = FALSE)
  value <- .rc_step_monitor_finish(list(ok = TRUE), monitor)
  expect_false(file.exists(file.path(outdir, "step_timing.tsv")))
  saveRDS(value, file.path(outdir, "step_grn.rds"))
  expect_message(
    .rc_step_monitor_fail(monitor),
    "RegCompass timing: grn \\[success\\]"
  )
  expect_false(file.exists(file.path(outdir, "step_timing.tsv")))

  failed_outdir <- tempfile("regcompass-stage-fail-")
  failed <- .rc_step_monitor_start("layer1", failed_outdir, progress = FALSE)
  .rc_step_monitor_finish(list(ok = TRUE), failed)
  expect_message(
    .rc_step_monitor_fail(failed),
    "RegCompass timing: layer1 \\[error\\]"
  )
  expect_false(file.exists(file.path(failed_outdir, "step_timing.tsv")))
})

test_that("execution timing writer removes stale timing files", {
  outdir <- tempfile("regcompass-execution-timing-")
  dir.create(outdir, recursive = TRUE)
  stale <- file.path(outdir, "00_execution_timing.tsv")
  writeLines("stale", stale)
  timing <- data.frame(stage = "x", elapsed_seconds = 1)
  expect_invisible(.rc_write_execution_timing(timing, outdir))
  expect_false(file.exists(stale))
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
