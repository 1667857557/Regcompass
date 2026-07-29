documentation_root <- function() {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  candidates <- candidates[
    file.exists(file.path(candidates, "DESCRIPTION"))
  ]
  if (!length(candidates)) return(NULL)
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

read_documentation <- function(paths) {
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
    rd,
    gregexpr("\\\\alias\\{[^}]+\\}", rd, perl = TRUE)
  )[[1L]]
  aliases <- sub("^\\\\alias\\{", "", aliases)
  aliases <- sub("\\}$", "", aliases)
  expect_setequal(exports, aliases)
})

test_that("tutorials and Rd use the current condition-GRN vocabulary", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  tutorials <- file.path(root, "docs", sprintf(
    "tutorial-%02d-%s.md",
    1:5,
    c(
      "quick-start", "stepwise-audit", "advanced-restart",
      "targeted-reaction-remapping", "condition-differential-analysis"
    )
  ))
  vignette <- file.path(root, "vignettes", "regcompass-workflow.Rmd")
  rd_files <- list.files(
    file.path(root, "man"), pattern = "\\.Rd$", full.names = TRUE
  )
  expect_true(all(file.exists(c(tutorials, vignette, rd_files))))

  user_workflows <- read_documentation(c(
    file.path(root, "README.md"), tutorials, vignette
  ))
  all_docs <- read_documentation(c(
    file.path(root, "README.md"),
    list.files(file.path(root, "docs"), pattern = "\\.md$",
               full.names = TRUE),
    vignette, rd_files
  ))

  current <- c(
    "ConditionGRNFit v5",
    "condition × broad cell type", "outer-heldout",
    "condition_grn_fits",
    "condition_fit_status", "tf_peak_gene_condition",
    "tf_peak_gene_condition_effect", "common-support",
    "TF RNA × peak ATAC", "single global FASTCORE completion",
    "penalty / (omega * vmax)"
  )
  expect_true(all(vapply(
    current, grepl, logical(1), x = all_docs, fixed = TRUE
  )))

  retired <- c(
    "RegCompassR 1.8.3", "RegCompassR 1.9.3",
    "Pando 1.2.1", "ConditionGRNFit v4", "Pando 1.4.0",
    "shared_design_independent",
    "condition_multitask_grn.md",
    "condition × cell type Pando evidence",
    "condition-by-cell-type Pando GRNs", "local FASTCORE",
    "locally completed meta-modules", "significantly supported GEM target",
    "number of significant edges", "minimum adjusted P value"
  )
  expect_false(any(vapply(
    retired, grepl, logical(1), x = all_docs, fixed = TRUE
  )))
  expect_false(grepl("padj_threshold", user_workflows, fixed = TRUE))
  expect_false(grepl("require_padj", user_workflows, fixed = TRUE))
  expect_false(file.exists(file.path(root, "docs", "recent-pr-audit.md")))
  expect_false(file.exists(
    file.path(root, "docs", "condition_multitask_grn.md")
  ))
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

test_that("current condition fit status excludes retired aliases", {
  implementation <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  expect_match(implementation, "condition_fit_status = status", fixed = TRUE)
  expect_false(grepl("sample_status = status", implementation, fixed = TRUE))
})

test_that("condition GRN architecture documentation matches current functions", {
  root <- documentation_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  path <- file.path(root, "docs", "condition-comparable-grn.md")
  expect_true(file.exists(path))
  documentation <- paste(readLines(path, warn = FALSE), collapse = "\n")

  current <- c(
    "coupled sparse-group multitask objective",
    "not independent elastic-net fits",
    "absolute condition coefficients, not `Delta beta`",
    "condition-pooled OOF reliability",
    "broad-cell-type pooled robust projection scale",
    "not 1-homogeneous",
    "condition_grn_fit_v5.rds"
  )
  expect_true(all(vapply(
    current, grepl, logical(1), x = documentation, fixed = TRUE
  )))

  retired <- c(
    "solves an ordinary elastic-net problem for each condition",
    "Condition coefficient columns are separable at fixed lambda",
    "Stage 4 uses only comparable condition effects",
    "No metacell-wise robust rescaling",
    "condition_grn_fit_v2.rds"
  )
  expect_false(any(vapply(
    retired, grepl, logical(1), x = documentation, fixed = TRUE
  )))
})
