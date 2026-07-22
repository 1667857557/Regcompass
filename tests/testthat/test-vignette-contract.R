test_that("workflow vignette documents the GRN-first API and tutorial levels", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) {
      file.path(workspace, "vignettes", "regcompass-workflow.Rmd")
    } else {
      character()
    },
    file.path("vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "..", "vignettes", "regcompass-workflow.Rmd")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) skip("Source vignette is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")

  expect_match(text, "RegCompassR 1.8.1", fixed = TRUE)
  expect_match(text, "Tutorial levels", fixed = TRUE)
  expect_match(text, "Level 1", fixed = TRUE)
  expect_match(text, "Level 2", fixed = TRUE)
  expect_match(text, "Level 3", fixed = TRUE)
  expect_match(text, "tutorial-01-quick-start.md", fixed = TRUE)
  expect_match(text, "tutorial-02-stepwise-audit.md", fixed = TRUE)
  expect_match(text, "tutorial-03-advanced-restart.md", fixed = TRUE)
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
  expect_match(text, "MulticoreParam", fixed = TRUE)
  expect_match(text, "local_fastcore_by_meta_module", fixed = TRUE)
  expect_match(text, "parallel = FALSE", fixed = TRUE)
  expect_match(text, "table(step5$feasible)", fixed = TRUE)
  expect_false(grepl("sample_balance = TRUE", text, fixed = TRUE))
  expect_false(grepl("min_metacells", text, fixed = TRUE))
})

test_that("three tutorial levels exist and have distinct Linux parallel scopes", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) dir.exists(file.path(path, "docs")),
    logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")
  root <- normalizePath(roots[[1L]], mustWork = TRUE)

  tutorial_paths <- file.path(
    root,
    "docs",
    c(
      "tutorial-01-quick-start.md",
      "tutorial-02-stepwise-audit.md",
      "tutorial-03-advanced-restart.md"
    )
  )
  expect_true(all(file.exists(tutorial_paths)))

  level1 <- paste(readLines(tutorial_paths[[1L]], warn = FALSE), collapse = "\n")
  level2 <- paste(readLines(tutorial_paths[[2L]], warn = FALSE), collapse = "\n")
  level3 <- paste(readLines(tutorial_paths[[3L]], warn = FALSE), collapse = "\n")

  expect_match(level1, "Tutorial Level 1", fixed = TRUE)
  expect_match(level1, "rc_run_regcompass_one_shot(", fixed = TRUE)
  expect_match(level1, "upstream_workers = upstream_workers", fixed = TRUE)
  expect_match(level1, "layer2_workers = layer2_workers", fixed = TRUE)
  expect_match(level1, "parallel_backend = \"multicore\"", fixed = TRUE)
  expect_match(level1, "automatically passed to SuperCell2", fixed = TRUE)
  expect_match(level1, "Confirm that the run completed", fixed = TRUE)

  expect_match(level2, "Tutorial Level 2", fixed = TRUE)
  expect_match(level2, "Stage map", fixed = TRUE)
  expect_match(level2, "BiocParallel::MulticoreParam", fixed = TRUE)
  expect_match(level2, "backend = \"multicore\"", fixed = TRUE)
  expect_match(level2, "local_fastcore_by_meta_module", fixed = TRUE)
  expect_match(level2, "shared-model × metacell", fixed = TRUE)
  expect_match(level2, "Gate before Stage 2 or 3", fixed = TRUE)
  expect_match(level2, "GRN/metacell group coverage", fixed = TRUE)
  expect_match(level2, "no separate label parameter is required", fixed = TRUE)

  expect_match(level3, "Tutorial Level 3", fixed = TRUE)
  expect_match(level3, "Linux process and thread controls", fixed = TRUE)
  expect_match(level3, "OMP_NUM_THREADS=1", fixed = TRUE)
  expect_match(level3, "REGCOMPASS_WORKERS=16", fixed = TRUE)
  expect_match(level3, "Parallel units by workflow stage", fixed = TRUE)
  expect_match(level3, "Minimal rerun matrix", fixed = TRUE)
  expect_match(level3, "Serial troubleshooting", fixed = TRUE)
  expect_match(level3, "Failure classification", fixed = TRUE)
  expect_match(level3, "Distinguish medium infeasibility from target blockage", fixed = TRUE)
  expect_match(level3, "automatically as", fixed = TRUE)

  combined <- paste(level1, level2, level3, collapse = "\n")
  expect_match(combined, "Pando_regcompass.tar.gz", fixed = TRUE)
  expect_match(combined, "peak_cor = 0.01", fixed = TRUE)
  expect_match(combined, "gamma = 75", fixed = TRUE)
  expect_match(combined, "pando_infer_args", fixed = TRUE)
  expect_match(combined, "parallel = FALSE", fixed = TRUE)
  expect_false(grepl("sample_balance = TRUE", combined, fixed = TRUE))
  expect_false(grepl("sample_balance_seed", combined, fixed = TRUE))
})

test_that("README and tutorial index expose Linux multicore controls", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) file.exists(file.path(path, "README.md")) &&
      file.exists(file.path(path, "docs", "run-modes-and-stepwise-workflow.md")),
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

  expect_match(text, "Choose a tutorial level", fixed = TRUE)
  expect_match(text, "Level 1", fixed = TRUE)
  expect_match(text, "Level 2", fixed = TRUE)
  expect_match(text, "Level 3", fixed = TRUE)
  expect_match(text, "tutorial-01-quick-start.md", fixed = TRUE)
  expect_match(text, "tutorial-02-stepwise-audit.md", fixed = TRUE)
  expect_match(text, "tutorial-03-advanced-restart.md", fixed = TRUE)
  expect_match(text, "upstream_workers", fixed = TRUE)
  expect_match(text, "layer2_workers", fixed = TRUE)
  expect_match(text, "MulticoreParam", fixed = TRUE)
  expect_match(text, "parallel_backend = \"multicore\"", fixed = TRUE)
  expect_match(text, "local FASTCORE", fixed = TRUE)
  expect_match(text, "shared-model", fixed = TRUE)
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
