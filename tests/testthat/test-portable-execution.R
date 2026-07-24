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

test_that("bundled human and mouse GEMs load without download", {
  human <- rc_prepare_gem(
    species = "human", version = "2.0.0", source = "bundled"
  )
  mouse <- rc_prepare_gem(
    species = "mouse", version = "1.8.0", source = "bundled"
  )
  expect_silent(rc_validate_species_gem(human, "human"))
  expect_silent(rc_validate_species_gem(mouse, "mouse"))
  expect_identical(human$model_info$distribution, "bundled_with_RegCompassR")
  expect_identical(mouse$model_info$distribution, "bundled_with_RegCompassR")
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
