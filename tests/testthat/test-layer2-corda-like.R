test_that("CORDA2 exposes only original adjustable parameters", {
  defaults <- RegCompassR:::.rc_layer2_corda_options(list())
  expect_identical(defaults$model_completion, "corda2")

  options <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda2",
    corda2_args = list(
      MCxNCthresh = 3,
      constraint = 20,
      constrainby = "perc",
      om = 1e5,
      ci = 0.02
    )
  ))
  expect_equal(options$MCxNCthresh, 3)
  expect_equal(options$constraint, 20)
  expect_identical(options$constrainby, "perc")
  expect_equal(options$om, 1e5)
  expect_equal(options$ci, 0.02)

  expect_error(
    RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = "corda2",
      corda2_args = list(n = 3L)
    )),
    "Allowed names"
  )
  expect_error(
    RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = "corda2",
      corda2_redundancies = 3L
    )),
    "Unsupported Python-CORDA"
  )
})

test_that("Layer 2 directly prepares and finalizes CORDA2", {
  implementation <- paste(
    deparse(body(rc_regcompass_step_layer2)), collapse = "\n"
  )
  expect_match(implementation, ".rc_prepare_corda_worker_pool", fixed = TRUE)
  expect_match(implementation, ".rc_layer2_prepare_completion", fixed = TRUE)
  expect_match(implementation, ".rc_layer2_finalize_completion", fixed = TRUE)
  expect_false(grepl("_base", implementation, fixed = TRUE))
  expect_false(grepl("before_", implementation, fixed = TRUE))
})

test_that("union cache directly dispatches original CORDA2 models to stage targets", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_build_celltype_medium_union_gem_cache_core)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    ".rc_build_celltype_medium_corda_cache",
    fixed = TRUE
  )
  helper <- paste(
    deparse(body(RegCompassR:::.rc_build_celltype_medium_corda_cache)),
    collapse = "\n"
  )
  expect_match(
    helper,
    "serial_cell_type_x_medium_models_stage_parallel_corda2_targets",
    fixed = TRUE
  )
  expect_match(
    helper,
    "serial_cell_type_x_medium_original_corda2",
    fixed = TRUE
  )
  expect_match(
    helper,
    "celltype_medium_original_matlab_corda2",
    fixed = TRUE
  )
})

test_that("CORDA2 is the default structural worker-template route", {
  expect_true(RegCompassR:::.rc_layer2_requested_corda2(list()))
  expect_true(RegCompassR:::.rc_layer2_requested_corda2(list(
    model_params = list()
  )))
  expect_false(RegCompassR:::.rc_layer2_requested_corda2(list(
    model_params = list(model_completion = "fastcore")
  )))
})

test_that("CORDA2 preparation leaves the template pool stopped", {
  skip_if_not_installed("BiocParallel")
  param <- BiocParallel::SnowParam(
    workers = 2L, type = "SOCK", progressbar = FALSE
  )
  state <- RegCompassR:::.rc_prepare_corda_worker_pool(
    layer2_args = list(), parallel = TRUE, BPPARAM = param
  )
  expect_identical(state$origin, "caller_supplied_stage_template")
  expect_false(isTRUE(state$started_here))
  expect_false(BiocParallel::bpisup(param))
  RegCompassR:::.rc_release_corda_worker_pool(state)
  expect_false(BiocParallel::bpisup(param))
})
