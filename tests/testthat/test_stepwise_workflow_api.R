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

test_that("stepwise computational stages expose optional parallelism", {
  stages <- list(
    meta_modules = rc_regcompass_step_meta_modules,
    layer1 = rc_regcompass_step_layer1,
    layer2 = rc_regcompass_step_layer2
  )
  for (stage in stages) {
    stage_formals <- formals(stage)
    expect_true(all(c("parallel", "BPPARAM") %in% names(stage_formals)))
    expect_identical(eval(stage_formals$parallel), TRUE)
    expect_null(eval(stage_formals$BPPARAM))
  }
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
    suppressWarnings(rc_regcompass_step_layer2(
      layer1 = list(),
      meta_modules = list(),
      gem = list(),
      medium_scenarios = data.frame(),
      outdir = tempfile()
    )),
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

test_that("stepwise Layer 2 delegates configured metadata columns", {
  body_text <- paste(
    deparse(body(.rc_regcompass_step_layer2_base)),
    collapse = "\n"
  )
  expect_match(body_text, "sample_col = params\\$sample_col")
  expect_match(body_text, "condition_col = params\\$condition_col")
  expect_match(body_text, "celltype_col = params\\$celltype_col")
})
