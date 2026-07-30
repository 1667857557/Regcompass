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
    file.path(root, "docs", c(
      "workflow.md",
      "functions.md",
      "mathematical-model.md",
      "condition-comparable-grn.md",
      "condition-comparability-safeguards.md",
      "tutorial-01-quick-start.md",
      "tutorial-02-stepwise-audit.md",
      "tutorial-03-advanced-restart.md",
      "tutorial-04-targeted-reaction-remapping.md",
      "tutorial-05-condition-differential-analysis.md"
    )),
    file.path(root, "vignettes", "regcompass-workflow.Rmd"),
    file.path(root, "man", c(
      "rc_regcompass_stepwise.Rd",
      "rc_run_regcompass.Rd",
      "rc_run_regcompass_one_shot.Rd"
    ))
  )
}

test_that("workflow vignette exposes the canonical executable workflow", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  path <- file.path(root, "vignettes", "regcompass-workflow.Rmd")
  expect_true(file.exists(path))
  text <- rc_read_doc(path)

  required <- c(
    "rc_prepare_gem(",
    "rc_make_medium_scenarios(",
    "rc_run_regcompass_one_shot(",
    'candidate_screen = "motif_domain"',
    "condition_mix = 0.5",
    "gamma = 30L",
    'comparison_support = "auto"',
    "regulatory_alpha = 1",
    'gpr_and_method = "min"',
    "rc_regcompass_step_grn(",
    "rc_regcompass_step_metacells(",
    "rc_regcompass_step_meta_modules(",
    "rc_regcompass_step_layer1(",
    "grn = step1",
    "rc_regcompass_step_layer2(",
    "rc_regcompass_step_results(",
    "rc_test_condition_reactions(",
    "Mathematical model"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
  expect_false(grepl("reference_condition", text, fixed = TRUE))
})

test_that("quick-start and stepwise tutorials retain runnable essentials", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- file.path(root, "docs", c(
    "tutorial-01-quick-start.md",
    "tutorial-02-stepwise-audit.md"
  ))
  text <- lapply(paths, rc_read_doc)

  shared <- c(
    'candidate_screen = "motif_domain"',
    "min_model_rsq = 0.1",
    'rna_reduction = "pca"',
    'atac_reduction = "lsi"',
    "gamma = 30L",
    "seed = 12345L",
    "regulatory_alpha = 1",
    'gpr_and_method = "min"',
    "[functions.md](functions.md)"
  )
  for (one in text) {
    expect_true(all(vapply(shared, grepl, logical(1), x = one, fixed = TRUE)))
    expect_false(grepl("reference_condition", one, fixed = TRUE))
  }

  expect_match(text[[1L]], "rc_run_regcompass_one_shot(", fixed = TRUE)
  expect_match(text[[1L]], "mm10_regulatory_regions.rds", fixed = TRUE)
  expect_match(text[[2L]], "rc_regcompass_step_grn(", fixed = TRUE)
  expect_match(text[[2L]], "BPPARAM = upstream_bp", fixed = TRUE)
  expect_match(text[[2L]], "grn = step1", fixed = TRUE)
  expect_match(text[[2L]], "rc_regcompass_step_results(", fixed = TRUE)
})

test_that("current documentation rejects retired or unsupported routes", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- rc_current_user_docs(root)
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  forbidden <- c(
    "reference_condition",
    "reference-condition coefficient",
    "TF–peak–GEM-gene",
    "local_fastcore",
    "local_fastcore_args",
    "global_modules",
    "previous_union_membership",
    "workflow_z_union_gem",
    "rc_build_meta_module_gem",
    "rc_project_metabolic_grn",
    "metabolic_gene_edges",
    "condition_grn_fit_v2.rds",
    "Stage 4 uses only comparable condition effects",
    "No metacell-wise robust rescaling",
    "regulatory_alpha = 0.5"
  )
  expect_false(any(vapply(
    forbidden, grepl, logical(1), x = text, fixed = TRUE
  )))

  unsupported_assignments <- c(
    "(?m)^\\s*candidate_screen\\s*=\\s*[\"']pooled_within_condition[\"']",
    "(?m)^\\s*aggregate_rna_col\\s*=",
    "(?m)^\\s*aggregate_peaks_col\\s*=",
    "(?m)^\\s*projection_component\\s*=\\s*[\"'](shared|deviation)[\"']",
    "(?m)^\\s*tau\\s*=",
    "(?m)^\\s*depth_balance\\s*=\\s*TRUE"
  )
  expect_false(any(vapply(
    unsupported_assignments, grepl, logical(1), x = text, perl = TRUE
  )))
})

test_that("README API and help agree on enforced arguments", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "condition-comparable-grn.md"),
    file.path(root, "man", "rc_regcompass_stepwise.Rd"),
    file.path(root, "man", "rc_run_regcompass.Rd"),
    file.path(root, "man", "rc_run_regcompass_one_shot.Rd")
  )
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  required <- c(
    "ConditionGRNFit v5",
    "condition_grn_fits",
    "tf_peak_gene_condition",
    "motif_domain",
    "outer-heldout",
    "absolute condition",
    "structural zero",
    "RNA-only",
    "gamma = 30L",
    "regulatory_alpha = 1",
    "condition_mix",
    "condition_weight",
    'projection_component = "condition"',
    "pairwise_common",
    "global_common",
    "gpr_and_method",
    "GRanges"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
})
