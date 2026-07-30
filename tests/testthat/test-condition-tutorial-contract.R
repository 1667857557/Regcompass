test_that("condition statistics remain documented for multi-condition runs", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) {
      file.exists(file.path(path, "README.md")) &&
        file.exists(file.path(path, "docs", "condition-reaction-statistics.md"))
    },
    logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")
  root <- normalizePath(roots[[1L]], mustWork = TRUE)
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "run-modes-and-stepwise-workflow.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "condition-reaction-statistics.md"),
    file.path(root, "vignettes", "regcompass-workflow.Rmd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  expect_match(text, "rc_test_condition_reactions", fixed = TRUE)
  expect_match(text, "rc_plot_condition_reaction", fixed = TRUE)
  expect_match(text, "Kruskal-Wallis", fixed = TRUE)
  expect_match(text, "Wilcoxon", fixed = TRUE)
  expect_match(text, "metacell", fixed = TRUE)
  expect_match(text, "standard_pando", fixed = TRUE)
  expect_match(text, "condition_grn", fixed = TRUE)
})

test_that("single-condition documentation does not promise artificial tests", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "README.md") else character(),
    "README.md", file.path("..", "README.md"),
    file.path("..", "..", "README.md")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("README is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")
  expect_match(text, "condition_col = NULL", fixed = TRUE)
  expect_match(text, "condition_contrast` is empty", fixed = TRUE)
  expect_match(text, "reaction_ranking", fixed = TRUE)
  expect_false(grepl("manufactured condition coefficient", text, fixed = TRUE))
})
