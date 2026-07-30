documentation_root <- function() {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  candidates <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (!length(candidates)) return(NULL)
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

read_documentation <- function(paths) {
  paths <- paths[file.exists(paths)]
  paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
}

normalize_documentation <- function(x) {
  trimws(gsub("[[:space:]]+", " ", x))
}

test_that("all exported APIs have Rd aliases", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  namespace <- trimws(readLines(file.path(root, "NAMESPACE"), warn = FALSE))
  exports <- grep("^export\\(", namespace, value = TRUE)
  exports <- sub("^export\\((.*)\\)$", "\\1", exports)
  rd_files <- list.files(
    file.path(root, "man"), pattern = "\\.Rd$", full.names = TRUE
  )
  rd <- read_documentation(rd_files)
  aliases <- regmatches(
    rd, gregexpr("\\\\alias\\{[^}]+\\}", rd, perl = TRUE)
  )[[1L]]
  aliases <- sub("^\\\\alias\\{", "", aliases)
  aliases <- sub("\\}$", "", aliases)
  expect_setequal(exports, aliases)
})

test_that("primary documentation describes routing, graph scope and condition-full OOF", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  docs <- normalize_documentation(read_documentation(c(
    file.path(root, "README.md"),
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md"),
    file.path(root, "docs", "stage-interface-contracts.md"),
    file.path(root, "docs", "metacell-graph-contract.md")
  )))
  required <- c(
    "standard_pando",
    "condition_grn",
    "Pando::infer_grn()",
    "pando_condition_grn_fit",
    "condition_full_oof",
    "projectable structural zero",
    "SCimplify_by_graph_group_from_embedding",
    "cell.graph.group",
    "cell.split.condition",
    "one_independent_graph_per_cell_type",
    "all_conditions_joint_within_cell_type_graph",
    "temporary_combined_stratum = FALSE",
    "No condition coefficients"
  )
  for (term in required) {
    expect_match(docs, term, fixed = TRUE, info = term)
  }
})

test_that("primary documentation rejects obsolete runtime and guardrail schemas", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  docs <- read_documentation(c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "stage-interface-contracts.md"),
    file.path(root, "DESCRIPTION")
  ))
  retired <- c(
    "ConditionGRNFit v5",
    "zzz00_absolute_pando_contract.R",
    "zzz01_fixed_gamma_metacells.R",
    "zzz02_layer1_policy.R",
    "zzz03_compass_gpr_penalty.R",
    "zzz04_canonical_pando_fit_schema.R",
    "penalty_depth_matched_rna",
    "penalty_common_depth_interval_rna",
    "penalty_alpha_sensitivity",
    "reaction_zero_support_sensitivity",
    "reaction_link_saturation_sensitivity"
  )
  expect_false(any(vapply(
    retired, grepl, logical(1), x = docs, fixed = TRUE
  )))
})

test_that("five canonical tutorials include targeted remapping", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  expected <- c(
    "tutorial-01-quick-start.md",
    "tutorial-02-stepwise-audit.md",
    "tutorial-03-mathematical-model.md",
    "tutorial-04-targeted-reaction-remapping.md",
    "tutorial-05-condition-differential-analysis.md"
  )
  tutorials <- file.path(root, "docs", expected)
  expect_true(all(file.exists(tutorials)))
  expect_true(all(vapply(
    tutorials,
    function(path) grepl(
      "[functions.md](functions.md)",
      paste(readLines(path, warn = FALSE), collapse = "\n"),
      fixed = TRUE
    ),
    logical(1)
  )))
  removed <- c(
    "tutorial-03-advanced-restart.md",
    "tutorial-04-condition-differential-analysis.md"
  )
  expect_false(any(file.exists(file.path(root, "docs", removed))))
})

test_that("targeted remapping retains exact model and condition-full contracts", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- read_documentation(c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "tutorial-04-targeted-reaction-remapping.md"),
    file.path(root, "docs", "target-union-scoring.md"),
    file.path(root, "docs", "workflow.md")
  ))
  required <- c(
    "rc_regcompass_step_target_union",
    "reaction_expression_condition_full_oof",
    "exact cached Stage 5 union GEM",
    "does not rerun FASTCORE",
    "KEGG",
    "Reactome",
    "master-Rhea"
  )
  for (term in required) {
    expect_match(text, term, fixed = TRUE, info = term)
  }
})

test_that("mathematical details remain centralized in Tutorial 3", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  mathematical <- file.path(
    root, "docs", "tutorial-03-mathematical-model.md"
  )
  expect_true(file.exists(mathematical))
  text <- paste(readLines(mathematical, warn = FALSE), collapse = "\n")
  required <- c(
    "Sparse-group multitask objective",
    "Estimability and projectable structural zeros",
    "Primary condition-full OOF projection",
    "Reliability and calibration",
    "GPR aggregation and reaction penalty",
    "Shared metabolic model"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  other_tutorials <- list.files(
    file.path(root, "docs"),
    pattern = "^tutorial-0[1245].*\\.md$",
    full.names = TRUE
  )
  other_text <- read_documentation(other_tutorials)
  expect_false(grepl("\\\\sum", other_text, fixed = TRUE))
})
