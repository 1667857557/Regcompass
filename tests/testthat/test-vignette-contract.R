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

rc_expect_doc_terms <- function(text, required) {
  missing <- required[!vapply(
    required, grepl, logical(1), x = text, fixed = TRUE
  )]
  expect_length(missing, 0L)
}

test_that("workflow vignette exposes canonical condition-full route and graph scope", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  path <- file.path(root, "vignettes", "regcompass-workflow.Rmd")
  expect_true(file.exists(path))
  text <- rc_read_doc(path)
  required <- c(
    "rc_prepare_gem(",
    "rc_make_medium_scenarios(",
    "rc_run_regcompass(",
    'condition_col = "Group"',
    "condition_grn",
    "standard_pando",
    "Pando::infer_grn()",
    "condition-full outer-heldout projection",
    "cell.graph.group",
    "cell.split.condition",
    "gamma = 30L",
    "regulatory_alpha = 1",
    'gpr_and_method = "min"',
    "penalty_condition_full_oof",
    "penalty_common_oof",
    "penalty_condition_unique_increment",
    "Tutorial 3"
  )
  rc_expect_doc_terms(text, required)
})

test_that("workflow vignette exposes supported media and user composition", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- rc_read_doc(file.path(root, "vignettes", "regcompass-workflow.Rmd"))
  required <- c(
    'scenario = "normal_human_plasma"',
    "mouse_plasma",
    "high_glucose",
    "low_glucose",
    "high_lactate",
    "low_lactate",
    "low_glutamine",
    "custom_medium",
    "custom_metabolites",
    "background_reference_doi",
    "challenge_reference_doi",
    "published_background_plus_named_nutrient_override"
  )
  rc_expect_doc_terms(text, required)
})

test_that("vignette rejects removed runtime and guardrail architecture", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- rc_read_doc(file.path(root, "vignettes", "regcompass-workflow.Rmd"))
  forbidden <- c(
    "ConditionGRNFit v5",
    "reference_condition",
    "regulatory_alpha = 0.5",
    "supercell_stratum_col",
    "cell.annotation",
    "zzz00_absolute_pando_contract",
    "zzz04_canonical_pando_fit_schema",
    "penalty_depth_matched_rna",
    "penalty_common_depth_interval_rna",
    "penalty_alpha_sensitivity",
    "reaction_zero_support_sensitivity",
    "reaction_link_saturation_sensitivity",
    'scenario = "physiologic"',
    'scenario = "minimal"',
    'scenario = "compass_model_bounds"',
    'scenario = "permissive_all_exchange"'
  )
  present <- forbidden[vapply(
    forbidden, grepl, logical(1), x = text, fixed = TRUE
  )]
  expect_length(present, 0L)
})

test_that("primary documentation agrees on enforced arguments and semantics", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "stage-interface-contracts.md"),
    file.path(root, "docs", "condition-comparable-grn.md"),
    file.path(root, "docs", "metacell-graph-contract.md"),
    file.path(root, "vignettes", "regcompass-workflow.Rmd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")
  required <- c(
    "standard_pando",
    "condition_grn",
    "motif_domain",
    "outer-heldout",
    "condition_full_oof",
    "projectable structural zero",
    "SCimplify_by_graph_group_from_embedding",
    "cell.graph.group",
    "cell.split.condition",
    "one_independent_graph_per_cell_type",
    "all_conditions_joint_within_cell_type_graph",
    "gamma = 30L",
    "regulatory_alpha = 1",
    "gpr_and_method"
  )
  rc_expect_doc_terms(text, required)
})