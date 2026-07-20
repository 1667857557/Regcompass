test_that("stepwise workflow functions are exported", {
  expected <- c(
    "rc_regcompass_step_metacells",
    "rc_regcompass_step_meta_modules",
    "rc_regcompass_step_layer1",
    "rc_regcompass_step_layer2",
    "rc_regcompass_step_results"
  )
  expect_true(all(expected %in% getNamespaceExports("RegCompassR")))
  expect_true(all(vapply(expected, function(name) {
    is.function(getExportedValue("RegCompassR", name))
  }, logical(1))))
})

test_that("stepwise stages reject incompatible upstream objects", {
  expect_error(
    rc_regcompass_step_meta_modules(
      metacells = list(),
      gem = list(),
      outdir = tempfile(),
      pfm = list(),
      genome = list()
    ),
    "output of `rc_regcompass_step_metacells\\(\\)`"
  )
  expect_error(
    rc_regcompass_step_layer2(
      layer1 = list(),
      meta_modules = list(),
      gem = list(),
      medium_scenarios = data.frame(),
      outdir = tempfile()
    ),
    "output of `rc_regcompass_step_meta_modules\\(\\)`"
  )
  expect_error(
    rc_regcompass_step_results(
      metacells = list(),
      meta_modules = list(),
      layer1 = list(),
      layer2 = list(),
      gem = list(),
      outdir = tempfile()
    ),
    "require outputs from the metacell and meta-module stages"
  )
})

test_that("stepwise Layer 2 preserves configured metadata columns", {
  body_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  expect_match(body_text, "sample_col = params\\$sample_col")
  expect_match(body_text, "condition_col = params\\$condition_col")
  expect_match(body_text, "celltype_col = params\\$celltype_col")
})
