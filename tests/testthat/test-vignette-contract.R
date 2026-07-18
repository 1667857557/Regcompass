test_that("workflow vignette follows the supported public API", {
  candidates <- c(
    file.path("vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "..", "vignettes", "regcompass-workflow.Rmd")
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    skip("Source vignette is unavailable in the installed-package test context.")
  }
  vignette_file <- candidates[[1L]]

  text <- paste(readLines(vignette_file, warn = FALSE), collapse = "\n")
  supported <- c(
    "rc_prepare_human2_gem",
    "rc_make_medium_scenarios",
    "rc_run_regcompass",
    "rc_run_regcompass_one_shot"
  )
  expect_true(all(vapply(
    supported,
    function(name) grepl(paste0(name, "\\("), text),
    logical(1)
  )))
  expect_match(text, 'inference_unit = "sample_celltype"')
  expect_match(text, "upstream_barrier\\$passed")
  expect_match(text, "microcompass\\$penalty")
})
