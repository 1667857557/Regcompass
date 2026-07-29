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
      "functions.md",
      "condition-comparable-grn.md",
      "condition-comparability-safeguards.md",
      "tutorial-01-quick-start.md",
      "tutorial-02-stepwise-audit.md"
    )),
    file.path(root, "vignettes", "regcompass-workflow.Rmd"),
    file.path(root, "man", c(
      "rc_regcompass_stepwise.Rd",
      "rc_run_regcompass.Rd",
      "rc_run_regcompass_one_shot.Rd"
    ))
  )
}

test_that("workflow vignette documents the current Pando bridge", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  path <- file.path(root, "vignettes", "regcompass-workflow.Rmd")
  expect_true(file.exists(path))
  text <- rc_read_doc(path)

  required <- c(
    "RegCompassR 2.1.0",
    "rc_run_regcompass_one_shot(",
    'candidate_screen = "motif_domain"',
    "condition_mix = 0.5",
    "cell_type",
    'reference_condition = "Control"',
    "comparison_mask",
    "comparable_to_reference",
    "pando_initiate_args",
    "pando_motif_args",
    "pando_infer_args",
    "Stage 2 owns",
    "rc_regcompass_step_grn(",
    "BPPARAM = upstream_bp",
    "step1$params$pando_parallel",
    "grn = step1",
    'comparison_support = "auto"',
    "condition × broad cell type",
    "outer-heldout",
    "gpr_and_method = \"min\"",
    "global FASTCORE",
    "Mouse runs must supply"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
})

test_that("the first two tutorials expose executable Pando routing", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- file.path(root, "docs", c(
    "tutorial-01-quick-start.md",
    "tutorial-02-stepwise-audit.md"
  ))
  expect_true(all(file.exists(paths)))
  text <- lapply(paths, rc_read_doc)

  for (one in text) {
    expect_match(one, 'candidate_screen = "motif_domain"', fixed = TRUE)
    expect_match(one, "cell_type", fixed = TRUE)
    expect_match(one, 'reference_condition = "Control"', fixed = TRUE)
    expect_match(one, "comparison_mask", fixed = TRUE)
    expect_match(one, "comparable_to_reference", fixed = TRUE)
    expect_match(one, "mm10_regulatory_regions.rds", fixed = TRUE)
    expect_match(one, "min_model_rsq = 0.1", fixed = TRUE)
    expect_match(one, 'rna_reduction = "pca"', fixed = TRUE)
    expect_match(one, 'atac_reduction = "lsi"', fixed = TRUE)
    expect_match(one, "seed = 12345L", fixed = TRUE)
    expect_match(one, 'gpr_and_method = "min"', fixed = TRUE)
  }

  expect_match(text[[1L]], "rc_run_regcompass_one_shot(", fixed = TRUE)
  expect_match(text[[1L]], "upstream_workers = 6L", fixed = TRUE)
  expect_match(text[[2L]], "rc_regcompass_step_grn(", fixed = TRUE)
  expect_match(text[[2L]], "BPPARAM = upstream_bp", fixed = TRUE)
  expect_match(text[[2L]], "Pando native map backend", fixed = TRUE)
  expect_match(
    text[[2L]], "condition × broad cell type", fixed = TRUE
  )
  expect_match(text[[2L]], "rc_regcompass_step_results(", fixed = TRUE)
})

test_that("canonical examples do not reintroduce unsafe Pando assignments", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- rc_current_user_docs(root)
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  forbidden_names <- c(
    "local_fastcore",
    "local_fastcore_args",
    "global_modules",
    "previous_union_membership",
    "workflow_z_union_gem",
    "rc_build_meta_module_gem",
    "rc_project_metabolic_grn",
    "metabolic_gene_edges"
  )
  expect_false(any(vapply(
    forbidden_names, grepl, logical(1), x = text, fixed = TRUE
  )))

  forbidden_assignments <- c(
    "(?m)^\\s*candidate_screen\\s*=\\s*[\"']pooled_within_condition[\"']",
    "(?m)^\\s*aggregate_rna_col\\s*=",
    "(?m)^\\s*aggregate_peaks_col\\s*=",
    "(?m)^\\s*expansion_mode\\s*=",
    "(?m)^\\s*max_iterations\\s*=",
    "(?m)^\\s*tau\\s*="
  )
  expect_false(any(vapply(
    forbidden_assignments, grepl, logical(1), x = text, perl = TRUE
  )))
  expect_false(grepl(
    "mouse: phastConsElements20Mammals.UCSC.hg38 only",
    text,
    fixed = TRUE
  ))
})

test_that("README API index and help agree on the Pando contract", {
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
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")

  required <- c(
    "RegCompassR 2.1.0",
    "condition_grn_fits",
    "tf_peak_gene_condition_effect",
    "comparison_mask",
    "comparable_to_reference",
    "motif_domain",
    "ConditionGRNFit v5",
    "condition × broad cell type",
    "outer-heldout",
    "reference_condition",
    "condition_mix",
    "condition_weight",
    "scale = TRUE",
    "pando_initiate_args",
    "pando_motif_args",
    "pando_infer_args",
    "BPPARAM = TRUE",
    "Pando native",
    "mouse",
    "GRanges",
    "complete-GPR",
    "global FASTCORE",
    "gpr_and_method"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )))
})
