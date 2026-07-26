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
      "rc_run_regcompass.Rd",
      "rc_run_regcompass_one_shot.Rd"
    ))
  )
}

test_that("canonical documentation describes bootstrap-stable 1.8.8 GRNs", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md"),
    file.path(root, "docs", "multitask-shared-grn.md"),
    file.path(root, "man", "rc_run_regcompass.Rd"),
    file.path(root, "man", "rc_regcompass_stepwise.Rd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  required <- c(
    "RegCompassR 1.8.8",
    "multitask_shared_backbone",
    "legacy_condition_pando",
    "prepare_grn_design",
    "tf_peak_gene_candidates",
    "tf_peak_gene_global",
    "tf_peak_gene_condition_all",
    "tf_peak_gene_significant",
    "condition_target_genes",
    "selection_frequency",
    "sign_stability",
    "n_bootstrap",
    "with replacement",
    "re-centred",
    "complete-GPR",
    "medium-specific union GEM",
    "same stoichiometric",
    "max_edges_per_target = Inf",
    "alpha = 0.5",
    "gpr_and_method = \"min\"",
    "completion_time_limit"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
})

test_that("quick-start and stepwise tutorials use the current Stage 1 contract", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- file.path(root, "docs", c(
    "tutorial-01-quick-start.md",
    "tutorial-02-stepwise-audit.md"
  ))
  text <- lapply(paths, rc_read_doc)

  for (value in text) {
    expect_match(value, 'grn_mode = "multitask_shared_backbone"', fixed = TRUE)
    expect_match(value, "pando_design_args = list(", fixed = TRUE)
    expect_match(value, "multitask_args = list(", fixed = TRUE)
    expect_match(value, "n_bootstrap = 100", fixed = TRUE)
    expect_match(value, "candidate_screen_threshold = 0", fixed = TRUE)
    expect_match(value, "max_edges_per_target = Inf", fixed = TRUE)
    expect_match(value, "seed = 12345L", fixed = TRUE)
    expect_false(grepl("sample_col", value, fixed = TRUE))
    expect_false(grepl("n_stability", value, fixed = TRUE))
    expect_false(grepl("stability_fraction", value, fixed = TRUE))
  }
  expect_match(text[[1L]], "rc_run_regcompass_one_shot(", fixed = TRUE)
  expect_match(text[[2L]], "rc_regcompass_step_results(", fixed = TRUE)
})

test_that("current user examples contain no retired or sample-column assignments", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- rc_current_user_docs(root)
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  forbidden_names <- c(
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
    "available_in_all_cached_union_models",
    "top_k_neighbors",
    "min_shared_tfs",
    "min_tf_jaccard",
    "max_targets_per_tf",
    "rc_project_metabolic_grn",
    "metabolic_gene_edges",
    "n_stability",
    "stability_fraction",
    "sample-blocked",
    "sample-aware"
  )
  expect_false(any(vapply(
    forbidden_names, grepl, logical(1), x = text, fixed = TRUE
  )))

  forbidden_assignments <- c(
    "(?m)^\\s*sample_col\\s*=",
    "(?m)^\\s*expansion_mode\\s*=",
    "(?m)^\\s*max_iterations\\s*=",
    "(?m)^\\s*tau\\s*=",
    "(?m)^\\s*and_method\\s*=\\s*[\"']boltzmann[\"']"
  )
  expect_false(any(vapply(
    forbidden_assignments, grepl, logical(1), x = text, perl = TRUE
  )))
})

test_that("generated help exposes bootstrap controls and no sample formal", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- file.path(root, "man", c(
    "rc_regcompass_stepwise.Rd",
    "rc_run_regcompass.Rd",
    "rc_run_regcompass_one_shot.Rd"
  ))
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  required <- c(
    "grn_mode",
    "multitask_args",
    "multitask_shared_backbone",
    "legacy_condition_pando",
    "n_bootstrap",
    "with replacement",
    "condition deviations",
    "complete GPR",
    "union GEM"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
  expect_false(grepl("sample_col", text, fixed = TRUE))
})
