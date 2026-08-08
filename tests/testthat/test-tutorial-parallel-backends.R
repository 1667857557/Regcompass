test_that("stepwise tutorial documents one portable worker budget", {
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
    "workers <- 10L",
    "options(RegCompassR.workers = 10L)",
    "max(1, available CPUs - 2)",
    "workers = workers",
    "Windows uses a SOCK cluster",
    "Linux/macOS use multicore workers"
  )
  missing <- required[!vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )]
  expect_length(missing, 0L)
  expect_false(grepl("upstream_bp", text, fixed = TRUE))
  expect_false(grepl("layer2_bp", text, fixed = TRUE))
  expect_false(grepl("BPPARAM =", text, fixed = TRUE))
})
