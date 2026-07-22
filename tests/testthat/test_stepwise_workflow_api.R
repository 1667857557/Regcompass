test_that("stepwise workflow functions are exported", {
  expected <- c(
    "rc_regcompass_step_grn", "rc_regcompass_step_metacells",
    "rc_regcompass_step_meta_modules", "rc_regcompass_step_layer1",
    "rc_regcompass_step_layer2", "rc_regcompass_step_results"
  )
  expect_true(all(expected %in% getNamespaceExports("RegCompassR")))
  expect_true(all(vapply(expected, function(name) {
    is.function(getExportedValue("RegCompassR", name))
  }, logical(1))))
})

test_that("single-cell GRN and computational stages expose optional parallelism", {
  stages <- list(
    grn = rc_regcompass_step_grn,
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

test_that("meta-module stage requires both GRN and metacell outputs", {
  expect_error(
    rc_regcompass_step_meta_modules(
      grn = list(), metacells = list(), gem = list(), outdir = tempfile()
    ),
    "output of `rc_regcompass_step_grn\\(\\)`"
  )
  fake_grn <- structure(list(), class = c("regcompass_grn_step", "list"))
  expect_error(
    rc_regcompass_step_meta_modules(
      grn = fake_grn, metacells = list(), gem = list(), outdir = tempfile()
    ),
    "output of `rc_regcompass_step_metacells\\(\\)`"
  )
})

test_that("meta-module stage no longer runs Pando", {
  f <- names(formals(rc_regcompass_step_meta_modules))
  expect_true(all(c("grn", "metacells", "gem", "outdir") %in% f))
  expect_false(any(c("pfm", "genome", "pando_args", "parallel", "BPPARAM") %in% f))
})

test_that("one-shot workflow refreshes the restartable Stage 6 result", {
  body_text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")

  stage6_save <- paste0(
    'file.path(outdir, "06_results", "regcompass_result.rds")'
  )
  root_save <- 'file.path(outdir, "regcompass_result.rds")'

  expect_match(body_text, stage6_save, fixed = TRUE)
  expect_match(body_text, root_save, fixed = TRUE)
  expect_lt(
    regexpr(stage6_save, body_text, fixed = TRUE)[[1L]],
    regexpr(root_save, body_text, fixed = TRUE)[[1L]]
  )
})
