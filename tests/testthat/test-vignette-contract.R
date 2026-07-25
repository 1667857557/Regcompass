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
    "merged_modules$merged_core_reactions",
    "merged_modules$merged_reaction_membership",
    "medium-specific union GEM",
    "single global FASTCORE completion",
    "layer2_args = list(",
    "model_params = list(",
    "fastcore_epsilon = 1e-4"
  )
  expect_true(all(vapply(
    required,
    grepl,
    logical(1),
    x = text,
    fixed = TRUE
  )))

  forbidden <- c(
    "layer1_args = list(local_fastcore",
    "local_fastcore_args = list(",
    "$global_modules",
    "$global_core_reactions",
    "$global_reaction_membership"
  )
  expect_false(any(vapply(
    forbidden,
    grepl,
    logical(1),
    x = text,
    fixed = TRUE
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
  expect_match(text[[3L]], "Rerun Stage 5 onward", fixed = TRUE)
  expect_match(text[[4L]], "rc_regcompass_step_target_union(", fixed = TRUE)
  expect_match(text[[5L]], "rc_test_condition_reactions(", fixed = TRUE)
  expect_match(combined, "medium-specific union GEM", fixed = TRUE)
  expect_match(combined, "global FASTCORE", fixed = TRUE)
  expect_match(combined, "peak_cor = 0.01", fixed = TRUE)
  expect_match(combined, "gamma = 30", fixed = TRUE)

  forbidden <- c(
    "local_fastcore = TRUE",
    "local_fastcore_args = list(",
    "$global_modules",
    "$global_core_reactions",
    "$global_reaction_membership",
    "global union meta-module GEM"
  )
  expect_false(any(vapply(
    forbidden,
    grepl,
    logical(1),
    x = combined,
    fixed = TRUE
  )))
})

test_that("README API index and Rd files expose current terminology", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "workflow.md"),
    file.path(root, "docs", "stage-interface-contracts.md"),
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
    "layer2_args$model_params"
  )
  expect_true(all(vapply(
    required,
    grepl,
    logical(1),
    x = text,
    fixed = TRUE
  )))

  forbidden <- c(
    "Worker count for GRN inference, local FASTCORE",
    "Global union meta-module GEM",
    "local FASTCORE completion."
  )
  expect_false(any(vapply(
    forbidden,
    grepl,
    logical(1),
    x = text,
    fixed = TRUE
  )))
})
