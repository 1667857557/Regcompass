test_that("Stage 1 min_cells defaults to 500 and preserves overrides", {
  defaulted <- RegCompassR:::.rc_resolve_stage1_min_cells_contract(list())
  expect_identical(defaulted$min_cells, 500L)
  expect_identical(defaulted$pando_args$min_cells, 500L)

  overridden <- RegCompassR:::.rc_resolve_stage1_min_cells_contract(
    list(min_cells = 300L)
  )
  expect_identical(overridden$min_cells, 300L)
  expect_identical(overridden$pando_args$min_cells, 300L)
})

test_that("Stage 1 min_cells requires a positive integer", {
  expect_error(
    RegCompassR:::.rc_resolve_stage1_min_cells_contract(list(min_cells = 0L)),
    "positive integer"
  )
  expect_error(
    RegCompassR:::.rc_resolve_stage1_min_cells_contract(list(min_cells = 10.5)),
    "positive integer"
  )
})

test_that("Stage 2 public route defaults min_cells_per_stratum to 500", {
  stage2 <- paste(deparse(body(rc_regcompass_step_metacells)), collapse = "\n")
  expect_match(
    stage2,
    "if (is.null(metacell_core_args$min_cells_per_stratum))",
    fixed = TRUE
  )
  expect_match(
    stage2,
    "metacell_core_args$min_cells_per_stratum <- 500L",
    fixed = TRUE
  )
})

test_that("Stage 1 provenance no longer declares min_cells fixed", {
  stage1 <- paste(deparse(body(rc_regcompass_step_grn)), collapse = "\n")
  expect_match(stage1, 'fixed = FALSE', fixed = TRUE)
  expect_match(stage1, 'configurable = TRUE', fixed = TRUE)
  expect_match(
    stage1,
    'threshold_source <- "configurable_pando_args_min_cells"',
    fixed = TRUE
  )
})
