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

test_that("computational stages expose optional parallelism", {
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

test_that("meta-module stage requires GRN and metacell outputs", {
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

test_that("one-shot workflow executes and saves Stage 6 output", {
  body_text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  expect_match(body_text, "rc_regcompass_step_results", fixed = TRUE)
  expect_match(
    body_text,
    'saveRDS(result, file.path(outdir, "regcompass_result.rds"))',
    fixed = TRUE
  )
})

test_that("GRN and metacell stages share automatic design resolution", {
  grn_text <- paste(c(
    deparse(body(rc_regcompass_step_grn)),
    deparse(body(RegCompassR:::.rc_original_step_grn_cell_set_contract))
  ), collapse = "\n")
  metacell_text <- paste(c(
    deparse(body(rc_regcompass_step_metacells)),
    deparse(body(RegCompassR:::.rc_original_step_metacells_cell_set_contract))
  ), collapse = "\n")
  resolver_text <- paste(
    deparse(body(.rc_resolve_condition_design)), collapse = "\n"
  )
  expect_match(grn_text, ".rc_resolve_condition_design", fixed = TRUE)
  expect_match(metacell_text, ".rc_resolve_condition_design", fixed = TRUE)
  expect_match(resolver_text, '"standard_pando"', fixed = TRUE)
  expect_match(resolver_text, '"condition_grn"', fixed = TRUE)
})
