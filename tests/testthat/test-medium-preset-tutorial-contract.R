medium_tutorial_root <- function() {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  candidates <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (!length(candidates)) return(NULL)
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

read_medium_tutorials <- function(paths) {
  paths <- paths[file.exists(paths)]
  paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
}

test_that("README begins with a runnable publication-bound minimal workflow", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  readme <- paste(
    readLines(file.path(root, "README.md"), warn = FALSE),
    collapse = "\n"
  )
  minimal_position <- regexpr("## Minimal workflow", readme, fixed = TRUE)[[1L]]
  details_position <- regexpr("For condition-aware analysis", readme,
                              fixed = TRUE)[[1L]]
  expect_gt(minimal_position, 0L)
  expect_gt(details_position, minimal_position)
  required <- c(
    "rc_prepare_gem(",
    "rc_make_medium_scenarios(",
    'scenario = "cantor2017_hplm"',
    "rc_run_regcompass_one_shot("
  )
  for (term in required) {
    expect_match(readme, term, fixed = TRUE, info = term)
  }
})

test_that("canonical tutorials retain the publication-bound medium contract", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")

  tutorials <- read_medium_tutorials(c(
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md")
  ))

  required <- c(
    "rc_make_medium_scenarios(",
    "cantor2017_hplm",
    "10.1016/j.cell.2017.03.023",
    "scenario = NULL",
    "reference_label",
    "reference_doi",
    "medium-presets.md",
    "uptake_scale",
    "Ambiguous salt-to-free-ion conversions"
  )

  for (term in required) {
    expect_match(tutorials, term, fixed = TRUE, info = term)
  }
})

test_that("published medium reference documents retained and removed scopes", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")

  path <- file.path(root, "docs", "medium-presets.md")
  expect_true(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  required <- c(
    "Published extracellular medium scenarios",
    "cantor2017_hplm",
    "published_paper_bound_presets_only",
    "manufacturer-only formulations",
    "Partial implementations are not retained",
    "10.1016/j.cell.2017.03.023"
  )
  for (term in required) {
    expect_match(text, term, fixed = TRUE, info = term)
  }
})