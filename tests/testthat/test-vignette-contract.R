test_that("workflow vignette documents the GRN-first API and required inputs", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "vignettes", "regcompass-workflow.Rmd") else character(),
    file.path("vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "..", "vignettes", "regcompass-workflow.Rmd")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("Source vignette is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")
  expect_match(text, "RegCompassR 1.8.2", fixed = TRUE)
  expect_match(text, "rc_regcompass_step_grn(", fixed = TRUE)
  expect_match(text, "rc_regcompass_step_metacells(", fixed = TRUE)
  expect_match(text, "rc_regcompass_step_meta_modules(", fixed = TRUE)
  expect_match(text, "peak_cor = 0.01", fixed = TRUE)
  expect_match(text, "gamma = 75", fixed = TRUE)
  expect_match(text, "RNA is normalized once", fixed = TRUE)
  expect_match(text, "ATAC TF-IDF", fixed = TRUE)
  expect_match(text, "condition as its only varying stratum", fixed = TRUE)
  expect_match(text, "dominant-cell-type ties are rejected", fixed = TRUE)
  expect_match(text, "Pando_regcompass.tar.gz", fixed = TRUE)
  expect_match(text, "rc_prepare_gem", fixed = TRUE)
  expect_match(text, "rc_make_medium_scenarios", fixed = TRUE)
  expect_match(text, "ChromatinAssay", fixed = TRUE)
  expect_match(text, "do not pass the `motif2tf`", fixed = TRUE)
  expect_match(text, "table(step5$feasible)", fixed = TRUE)
  expect_false(grepl("sample_balance = TRUE", text, fixed = TRUE))
  expect_false(grepl("min_metacells", text, fixed = TRUE))
})

test_that("README and full tutorial document local Pando and solver behavior", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) dir.exists(file.path(path, "man")) &&
      dir.exists(file.path(path, "docs")),
    logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")
  root <- normalizePath(roots[[1L]], mustWork = TRUE)
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "run-modes-and-stepwise-workflow.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "vignettes", "regcompass-workflow.Rmd"),
    file.path(root, "man", "rc_run_regcompass.Rd"),
    file.path(root, "man", "rc_run_regcompass_one_shot.Rd"),
    file.path(root, "man", "rc_regcompass_stepwise.Rd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  expect_match(text, "peak_cor = 0.01", fixed = TRUE)
  expect_match(text, "gamma = 75", fixed = TRUE)
  expect_match(text, "condition-only", fixed = TRUE)
  expect_match(text, "dominant", fixed = TRUE)
  expect_match(text, "Pando_regcompass.tar.gz", fixed = TRUE)
  expect_match(text, "GitHub remote metadata are not required", fixed = TRUE)
  expect_match(text, "solver installation", fixed = TRUE)
  expect_match(text, "core reactions", fixed = TRUE)
  expect_match(text, "master-Rhea", fixed = TRUE)
  expect_false(grepl("v170_sample_balance", text, fixed = TRUE))
  expect_false(grepl("sample_balance_seed", text, fixed = TRUE))
})
