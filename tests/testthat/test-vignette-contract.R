rc_doc_root <- function() {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots, function(path) file.exists(file.path(path, "DESCRIPTION")),
    logical(1)
  )]
  if (!length(roots)) return(NULL)
  normalizePath(roots[[1L]], mustWork = TRUE)
}

rc_read_doc <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

test_that("canonical documentation describes RegCompassR 1.8.9", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md"),
    file.path(root, "docs", "tutorial-05-condition-differential-analysis.md")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")
  required <- c(
    "RegCompassR 1.8.9", "multitask_shared_backbone",
    "condition-stratified", "with replacement", "selection_frequency",
    "sign_stability", "complete-GPR", "shared medium-specific union GEM",
    "candidate_screen_threshold = 0", "max_edges_per_target = Inf",
    "gpr_and_method = \"min\"", "compact", "table_manifest",
    "stage_provenance", "biological-replicate"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("tutorials form one continuous staged workflow", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- file.path(root, "docs", sprintf(
    "tutorial-%02d-%s.md",
    1:5,
    c(
      "quick-start", "stepwise-audit", "advanced-restart",
      "targeted-reaction-remapping", "condition-differential-analysis"
    )
  ))
  expect_true(all(file.exists(paths)))
  text <- lapply(paths, rc_read_doc)
  expect_match(text[[1L]], "Next:", fixed = TRUE)
  expect_match(text[[2L]], "Previous:", fixed = TRUE)
  expect_match(text[[2L]], "Next:", fixed = TRUE)
  expect_match(text[[3L]], "Previous:", fixed = TRUE)
  expect_match(text[[3L]], "Next:", fixed = TRUE)
  expect_match(text[[4L]], "Previous:", fixed = TRUE)
  expect_match(text[[4L]], "Next:", fixed = TRUE)
  expect_match(text[[5L]], "Previous:", fixed = TRUE)
  expect_match(text[[5L]], "bootstrap-active TF–peak–target edge", fixed = TRUE)
  expect_match(text[[5L]], "complete GPR branch", fixed = TRUE)
  expect_match(text[[5L]], "direction-specific LP support", fixed = TRUE)
})

test_that("quick-start and stepwise examples use current contracts", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- file.path(root, "docs", c(
    "tutorial-01-quick-start.md", "tutorial-02-stepwise-audit.md"
  ))
  for (value in lapply(paths, rc_read_doc)) {
    expect_match(value, 'grn_mode = "multitask_shared_backbone"', fixed = TRUE)
    expect_match(value, "pando_design_args = list(", fixed = TRUE)
    expect_match(value, "multitask_args = list(", fixed = TRUE)
    expect_match(value, "n_bootstrap = 100", fixed = TRUE)
    expect_match(value, "min_bootstrap_success_fraction = 0.8", fixed = TRUE)
    expect_match(value, "candidate_screen_threshold = 0", fixed = TRUE)
    expect_match(value, "max_edges_per_target = Inf", fixed = TRUE)
    expect_false(grepl("(?m)^\\s*sample_col\\s*=", value, perl = TRUE))
    expect_false(grepl("(?m)^\\s*pool_col\\s*=", value, perl = TRUE))
  }
})

test_that("documentation explains current SuperCell2 and compact results", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- paste(
    rc_read_doc(file.path(root, "README.md")),
    rc_read_doc(file.path(root, "docs", "tutorial-02-stepwise-audit.md")),
    collapse = "\n"
  )
  required <- c(
    "strata_cols", "label = celltype_col", "no artificial condition-pool",
    "reaction_catalog", "reaction_evidence", "detailed stage checkpoints"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
})
