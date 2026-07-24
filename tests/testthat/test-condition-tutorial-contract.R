test_that("condition statistics and plots are linked from tutorial entry points", {
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
  expect_match(text, "condition-reaction-statistics.md", fixed = TRUE)
  expect_match(text, "rc_test_condition_reactions(", fixed = TRUE)
  expect_match(text, "rc_plot_condition_reaction(", fixed = TRUE)
  expect_match(text, "Kruskal-Wallis", fixed = TRUE)
  expect_match(text, "Wilcoxon", fixed = TRUE)
  expect_match(text, "p_adjust_scope", fixed = TRUE)
  expect_match(text, "annotation_p = \"p_adj\"", fixed = TRUE)
  expect_match(text, "one point per metacell", fixed = TRUE)
  expect_match(text, "significance brackets", fixed = TRUE)
  expect_match(text, "metacell-level", fixed = TRUE)
})

test_that("README demonstrates multi-condition testing and plotting", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "README.md") else character(),
    "README.md", file.path("..", "README.md"), file.path("..", "..", "README.md")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("README is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")

  expect_match(text, "Compare the same reaction across conditions", fixed = TRUE)
  expect_match(text, "control_24hr", fixed = TRUE)
  expect_match(text, "JQ1_24hr", fixed = TRUE)
  expect_match(text, "MS177_24hr", fixed = TRUE)
  expect_match(text, "condition_stats$omnibus", fixed = TRUE)
  expect_match(text, "condition_stats$pairwise", fixed = TRUE)
  expect_match(text, "target_direction = \"reverse\"", fixed = TRUE)
  expect_match(text, "annotation_p = \"p_adj\"", fixed = TRUE)
})