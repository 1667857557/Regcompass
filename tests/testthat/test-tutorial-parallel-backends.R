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
