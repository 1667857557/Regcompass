test_that("workflow vignette documents the 1.8.2 staged API", {
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
  required <- c(
    "RegCompassR 1.8.2",
    "rc_prepare_gem",
    "rc_make_medium_scenarios",
    "scenario = \"physiologic\"",
    "Pando_regcompass.tar.gz",
    "ChromatinAssay",
    "peak_cor = 0.01",
    "gamma = 75",
    "rc_regcompass_step_grn(",
    "rc_regcompass_step_metacells(",
    "rc_regcompass_step_meta_modules(",
    "rc_regcompass_step_layer1(",
    "rc_regcompass_step_layer2(",
    "rc_regcompass_step_target_union(",
    "rc_regcompass_step_results(",
    "regcompass_layer1_step",
    "regcompass_layer2_step",
    "structural_model_reused_exactly",
    "result$version, \"1.8.2\""
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  forbidden <- c(
    "metacell_label_col",
    "label_col =",
    "sample_balance = TRUE",
    "_v170",
    "RegCompassR.inference_unit"
  )
  expect_false(any(vapply(forbidden, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("tutorials cover quick start, strict stage audit, and restart", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots, function(path) dir.exists(file.path(path, "docs")), logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")
  root <- normalizePath(roots[[1L]], mustWork = TRUE)
  paths <- file.path(root, "docs", c(
    "tutorial-01-quick-start.md",
    "tutorial-02-stepwise-audit.md",
    "tutorial-03-advanced-restart.md",
    "target-union-scoring.md"
  ))
  expect_true(all(file.exists(paths)))
  text <- lapply(paths, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  })
  expect_match(text[[1L]], "Tutorial Level 1", fixed = TRUE)
  expect_match(text[[1L]], "rc_run_regcompass_one_shot(", fixed = TRUE)
  expect_match(text[[2L]], "Tutorial Level 2", fixed = TRUE)
  expect_match(text[[2L]], "regcompass_grn_step", fixed = TRUE)
  expect_match(text[[2L]], "regcompass_layer1_step", fixed = TRUE)
  expect_match(text[[2L]], "regcompass_layer2_step", fixed = TRUE)
  expect_match(text[[2L]], "rc_regcompass_step_target_union(", fixed = TRUE)
  expect_match(text[[2L]], "structural_model_reused_exactly", fixed = TRUE)
  expect_match(text[[3L]], "Tutorial Level 3", fixed = TRUE)
  expect_match(text[[3L]], "Earliest stage to rerun", fixed = TRUE)
  expect_match(text[[3L]], "Serial troubleshooting", fixed = TRUE)
  expect_match(text[[4L]], "exact cached model file", fixed = TRUE)
  expect_match(text[[4L]], "source_model_md5", fixed = TRUE)
  combined <- paste(unlist(text), collapse = "\n")
  expect_match(combined, "peak_cor = 0.01", fixed = TRUE)
  expect_match(combined, "gamma = 75", fixed = TRUE)
  expect_match(combined, "OMP_NUM_THREADS=1", fixed = TRUE)
  expect_false(grepl("_v170", combined, fixed = TRUE))
})

test_that("README and API index expose current public workflow only", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) file.exists(file.path(path, "README.md")) &&
      file.exists(file.path(path, "docs", "functions.md")),
    logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")
  root <- normalizePath(roots[[1L]], mustWork = TRUE)
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "man", "rc_regcompass_stepwise.Rd"),
    file.path(root, "man", "rc_regcompass_step_target_union.Rd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  required <- c(
    "RegCompassR 1.8.2",
    "rc_run_regcompass_one_shot",
    "rc_regcompass_step_target_union",
    "GEM fingerprint",
    "ordered metacell IDs",
    "exact cached global union GEM",
    "master-Rhea",
    "medium-presets.md"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  forbidden <- c(
    "v170_microcompass_contract",
    "internal_apply",
    "metacell_label_col",
    "sample_balance_seed"
  )
  expect_false(any(vapply(forbidden, grepl, logical(1), x = text, fixed = TRUE)))
})
