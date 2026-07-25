test_that("stepwise tutorial defines portable BiocParallel backends", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")

  path <- file.path(root, "docs", "tutorial-02-stepwise-audit.md")
  expect_true(file.exists(path))
  text <- rc_read_doc(path)

  required <- c(
    "library(BiocParallel)",
    "upstream_workers <- 6L",
    "layer2_workers <- 30L",
    ".Platform$OS.type == \"windows\"",
    "SnowParam(",
    "type = \"SOCK\"",
    "MulticoreParam(",
    "BPPARAM = upstream_bp",
    "BPPARAM = layer2_bp"
  )

  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
})
