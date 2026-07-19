test_that("workflow vignette follows the v1.7.0 public API", {
  candidates <- c(
    file.path("vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "..", "vignettes", "regcompass-workflow.Rmd")
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    skip("Source vignette is unavailable in the installed-package test context.")
  }
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")
  supported <- c(
    "rc_prepare_gem", "rc_prepare_human2_gem", "rc_prepare_mouse_gem",
    "rc_make_medium_scenarios", "rc_run_regcompass",
    "rc_run_regcompass_one_shot"
  )
  expect_true(all(vapply(
    supported,
    function(name) grepl(paste0(name, "\\("), text),
    logical(1)
  )))
  expect_match(text, "pooled across biological samples", fixed = TRUE)
  expect_match(text, "uses ATAC", fixed = TRUE)
  expect_match(text, 'fragment_files = FALSE')
  expect_match(text, 'species = "human"')
  expect_match(text, 'species = "mouse"')
  expect_match(text, 'medium_scenario = "normal_human_plasma"')
  expect_match(text, "medium_scenarios = medium")
  expect_match(text, 'inference_unit = "metacell"')
  expect_match(text, "pando_initiate_args = list")
  expect_match(text, "regions = SCREEN.ccRE.UCSC.hg38")
  expect_match(text, "microcompass\\$penalty")
  expect_match(text, "condition_contrast")
})
