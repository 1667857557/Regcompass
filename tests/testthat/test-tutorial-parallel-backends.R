test_that("stepwise tutorial defines portable BiocParallel backends", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) file.exists(file.path(path, "DESCRIPTION")),
    logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")

  root <- normalizePath(roots[[1L]], mustWork = TRUE)
  path <- file.path(root, "docs", "tutorial-02-stepwise-audit.md")
  expect_true(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  required <- c(
    "library(BiocParallel)",
    ".Platform$OS.type == \"windows\"",
    "SnowParam(workers = 6L",
    "SnowParam(workers = 30L",
    "MulticoreParam(workers = 6L",
    "MulticoreParam(workers = 30L",
    "BPPARAM = upstream_bp",
    "BPPARAM = layer2_bp"
  )
  missing <- required[!vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )]
  expect_length(missing, 0L, info = paste("Missing tutorial entries:", paste(missing, collapse = ", ")))
})
