test_that("tutorials document the implemented TF-peak-gene penalty gate", {
  paths <- c(
    quick_start = testthat::test_path(
      "..", "..", "docs", "tutorial-01-quick-start.md"
    ),
    stepwise = testthat::test_path(
      "..", "..", "docs", "tutorial-02-stepwise-audit.md"
    ),
    mathematical = testthat::test_path(
      "..", "..", "docs", "mathematical-model.md"
    ),
    vignette = testthat::test_path(
      "..", "..", "vignettes", "regcompass-workflow.Rmd"
    )
  )

  documents <- lapply(paths, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  })

  for (document in documents) {
    expect_match(document, "padj < 0.05", fixed = TRUE)
    expect_match(document, "abs(corr) >= 0.05", fixed = TRUE)
    expect_match(document, "abs(estimate) >= 0.05", fixed = TRUE)
    expect_false(grepl("abs(corr) >= 0.1", document, fixed = TRUE))
    expect_false(grepl("abs(estimate) >= 0.01", document, fixed = TRUE))
  }

  expect_match(
    documents$quick_start,
    "tf_cor = 0.1` is a Pando candidate-discovery control",
    fixed = TRUE
  )
  expect_match(
    documents$stepwise,
    "tf_cor = 0.1` controls Pando candidate discovery",
    fixed = TRUE
  )
  expect_match(
    documents$mathematical,
    "\\rho_0=0.05,\\qquad \\beta_0=0.05",
    fixed = TRUE
  )
  expect_match(documents$stepwise, "corr_source", fixed = TRUE)
  expect_match(documents$stepwise, "penalty_effect", fixed = TRUE)
})
