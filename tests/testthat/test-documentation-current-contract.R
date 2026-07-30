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

test_that("primary documentation describes automatic Pando routing", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  docs <- read_documentation(c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "run-modes-and-stepwise-workflow.md"),
    file.path(root, "docs", "stage-interface-contracts.md")
  ))
  required <- c(
    "standard_pando",
    "condition_grn",
    "Pando::infer_grn()",
    "pando_condition_grn_fit",
    "cell.annotation",
    "cell.split.condition",
    "temporary_combined_stratum = FALSE",
    "No condition coefficients"
  )
  expect_true(all(vapply(
    required, grepl, logical(1), x = docs, fixed = TRUE
  )))
})

test_that("primary documentation rejects obsolete runtime implementation", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  docs <- read_documentation(c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "run-modes-and-stepwise-workflow.md"),
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
    "reference_condition"
  )
  expect_false(any(vapply(
    retired, grepl, logical(1), x = docs, fixed = TRUE
  )))
})

test_that("each tutorial links the current API index", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  tutorials <- list.files(
    file.path(root, "docs"), pattern = "^tutorial-0[1-5].*\\.md$",
    full.names = TRUE
  )
  expect_length(tutorials, 5L)
  expect_true(all(vapply(
    tutorials,
    function(path) grepl(
      "[functions.md](functions.md)",
      paste(readLines(path, warn = FALSE), collapse = "\n"),
      fixed = TRUE
    ),
    logical(1)
  )))
})

test_that("mathematical details remain centralized", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  mathematical <- file.path(root, "docs", "mathematical-model.md")
  expect_true(file.exists(mathematical))
  text <- paste(readLines(mathematical, warn = FALSE), collapse = "\n")
  required <- c(
    "Sparse-group multitask objective",
    "Outer-heldout regulatory projection",
    "Reliability and calibration",
    "GPR aggregation and reaction penalty",
    "Shared metabolic model"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
})
