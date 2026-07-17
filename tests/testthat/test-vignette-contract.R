test_that("workflow vignette follows the supported public API", {
  root <- if (dir.exists("vignettes")) "." else file.path("..", "..")
  vignette_file <- file.path(root, "vignettes", "regcompass-workflow.Rmd")
  expect_true(file.exists(vignette_file))

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
