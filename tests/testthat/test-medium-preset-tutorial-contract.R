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

test_that("README begins with a runnable minimal workflow", {
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
    'scenario = "normal_human_plasma"',
    "rc_run_regcompass_one_shot("
  )
  for (term in required) {
    expect_match(readme, term, fixed = TRUE, info = term)
  }
})

test_that("canonical tutorials retain all biological medium scenarios", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")

  tutorials <- read_medium_tutorials(c(
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md")
  ))

  scenarios <- c(
    "normal_human_plasma", "mouse_plasma", "high_glucose", "low_glucose",
    "high_lactate", "low_lactate", "low_glutamine", "custom"
  )
  required <- c(
    "rc_make_medium_scenarios(",
    scenarios,
    "scenario = NULL",
    "custom_medium",
    "custom_metabolites",
    "background_reference_doi",
    "challenge_reference_doi",
    "medium-presets.md"
  )
  for (term in required) {
    expect_match(tutorials, term, fixed = TRUE, info = term)
  }
})

test_that("medium reference documents evidence and composite construction", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")

  path <- file.path(root, "docs", "medium-presets.md")
  expect_true(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  required <- c(
    "Medium scenarios and published evidence",
    "published_background_plus_named_nutrient_override",
    "published_plasma_or_culture_background_with_explicit_overrides",
    "normal_human_plasma",
    "mouse_plasma",
    "high_glucose",
    "low_lactate",
    "low_glutamine",
    "User-defined medium composition",
    "Concentration is not uptake flux"
  )
  for (term in required) {
    expect_match(text, term, fixed = TRUE, info = term)
  }
})

test_that("technical boundary modes are not documented as biological scenarios", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  tutorials <- read_medium_tutorials(c(
    file.path(root, "README.md"),
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md")
  ))
  for (term in c("minimal", "compass_model_bounds", "permissive_all_exchange")) {
    expect_false(grepl(paste0('scenario = "', term, '"'), tutorials,
                       fixed = TRUE), info = term)
  }
})