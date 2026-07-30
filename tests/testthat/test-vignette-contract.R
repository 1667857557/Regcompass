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

test_that("workflow vignette exposes both Pando modes", {
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
    "condition_col = NULL",
    '"condition_grn"',
    '"standard_pando"',
    "Pando::infer_grn()",
    "pando_condition_grn_fit",
    "cell.annotation",
    "cell.split.condition",
    "gamma = 30L",
    "regulatory_alpha = 1",
    'gpr_and_method = "min"',
    "rc_regcompass_step_grn(",
    "rc_regcompass_step_metacells(",
    "rc_regcompass_step_meta_modules(",
    "rc_regcompass_step_layer1(",
    "rc_regcompass_step_layer2(",
    "rc_regcompass_step_results(",
    "Mathematical model"
  )
  rc_expect_doc_terms(text, required)
})

test_that("vignette rejects removed runtime architecture", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- rc_read_doc(file.path(root, "vignettes", "regcompass-workflow.Rmd"))
  forbidden <- c(
    "ConditionGRNFit v5",
    "reference_condition",
    "regulatory_alpha = 0.5",
    "supercell_stratum_col",
    "zzz00_absolute_pando_contract",
    "zzz04_canonical_pando_fit_schema"
  )
  present <- forbidden[vapply(
    forbidden, grepl, logical(1), x = text, fixed = TRUE
  )]
  expect_length(present, 0L)
})

test_that("primary documentation agrees on enforced arguments", {
  root <- rc_doc_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "stage-interface-contracts.md"),
    file.path(root, "docs", "condition-comparable-grn.md"),
    file.path(root, "vignettes", "regcompass-workflow.Rmd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, rc_read_doc)), collapse = "\n")
  required <- c(
    "standard_pando",
    "condition_grn",
    "condition_coefficients_calculated",
    "motif_domain",
    "outer-heldout",
    "cell.annotation",
    "cell.split.condition",
    "gamma = 30L",
    "regulatory_alpha = 1",
    "pairwise_common",
    "global_common",
    "gpr_and_method"
  )
  rc_expect_doc_terms(text, required)
})
