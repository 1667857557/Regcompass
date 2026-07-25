rc_doc_root <- function() {
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
  if (!length(roots)) return(NULL)
  normalizePath(roots[[1L]], mustWork = TRUE)
}

rc_read_doc <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

rc_current_user_docs <- function(root) {
  c(
    file.path(root, "README.md"),
    list.files(
      file.path(root, "docs"),
      pattern = "\\.md$",
      recursive = TRUE,
      full.names = TRUE
    ),
    file.path(root, "vignettes", "regcompass-workflow.Rmd"),
    file.path(root, "man", c(
      "rc_regcompass_stepwise.Rd",
      "rc_regcompass_step_target_union.Rd",
      "rc_run_regcompass.Rd"
    ))
  )
}

test_that("workflow vignette documents global-only FASTCORE", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  path <- file.path(root, "vignettes", "regcompass-workflow.Rmd")
  expect_true(file.exists(path))
  text <- rc_read_doc(path)

  required <- c(
    "RegCompassR 1.8.4",
    "rc_run_regcompass_one_shot(",
    "rc_regcompass_step_meta_modules(",
    "rc_regcompass_step_layer2(",
    "rc_regcompass_step_target_union(",
    "merged_modules$merged_core_reactions",
    "merged_modules$merged_reaction_membership",
    "medium-specific union GEM",
    "single global FASTCORE completion",
    "global_fastcore_support",
    "file_checksum",
    "structural_model_reused_exactly",
    "fastcore_rerun",
    "model_rebuild"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
})

test_that("all five tutorials use the current stage contract", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- file.path(root, "docs", sprintf(
    "tutorial-%02d-%s.md",
    1:5,
    c(
      "quick-start",
      "stepwise-audit",
      "advanced-restart",
      "targeted-reaction-remapping",
      "condition-differential-analysis"
    )
  ))
  expect_true(all(file.exists(paths)))
  text <- lapply(paths, rc_read_doc)
  combined <- paste(unlist(text), collapse = "\n")

  expect_match(text[[1L]], "rc_run_regcompass_one_shot(", fixed = TRUE)
  expect_match(text[[2L]], "rc_regcompass_step_results(", fixed = TRUE)
  expect_match(text[[2L]], "merged_modules", fixed = TRUE)
  expect_match(text[[3L]], "global_fastcore_support", fixed = TRUE)
  expect_match(text[[4L]], "rc_regcompass_step_target_union(", fixed = TRUE)
  expect_match(text[[4L]], "available_in_all_cached_union_gems", fixed = TRUE)
  expect_match(text[[4L]], "file_checksum", fixed = TRUE)
  expect_match(text[[4L]], "fastcore_rerun", fixed = TRUE)
  expect_match(text[[4L]], "model_rebuild", fixed = TRUE)
  expect_match(text[[5L]], "rc_test_condition_reactions(", fixed = TRUE)
  expect_match(combined, "medium-specific union GEM", fixed = TRUE)
  expect_match(combined, "global FASTCORE", fixed = TRUE)
  expect_match(combined, "peak_cor = 0.01", fixed = TRUE)
  expect_match(combined, "gamma = 30", fixed = TRUE)
})

test_that("user documentation contains no retired architecture names", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- rc_current_user_docs(root)
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  forbidden <- c(
    "local_fastcore",
    "local_fastcore_args",
    "global_modules",
    "global_core_reactions",
    "global_reaction_membership",
    "previous_union_membership",
    "GLOBAL_UNION",
    "workflow_z_union_gem",
    "rc_build_meta_module_gem",
    "reaction_meta$fastcore_support",
    "available_in_all_cached_union_models"
  )
  expect_false(any(vapply(
    forbidden, grepl, logical(1), x = text, fixed = TRUE
  )))
})

test_that("README API index and Rd files expose final model reuse", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "workflow.md"),
    file.path(root, "docs", "stage-interface-contracts.md"),
    file.path(root, "docs", "target-union-scoring.md"),
    file.path(root, "man", "rc_regcompass_stepwise.Rd"),
    file.path(root, "man", "rc_regcompass_step_target_union.Rd"),
    file.path(root, "man", "rc_run_regcompass.Rd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  required <- c(
    "RegCompassR 1.8.4",
    "merged_core_reactions",
    "merged_reaction_membership",
    "medium-specific union GEM",
    "global FASTCORE",
    "layer2_args$model_params",
    "file_checksum",
    "structural_model_reused_exactly",
    "fastcore_rerun",
    "model_rebuild"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
})
