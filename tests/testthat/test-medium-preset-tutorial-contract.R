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

test_that("canonical tutorials retain the complete medium preset interface", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")

  tutorials <- read_medium_tutorials(c(
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md")
  ))

  scenarios <- c(
    "physiologic",
    "normal_human_plasma",
    "mouse_plasma",
    "rpmi1640",
    "dmem_high_glucose",
    "high_glucose",
    "low_glucose",
    "high_lactate",
    "low_lactate",
    "low_glutamine",
    "minimal",
    "compass_model_bounds",
    "permissive_all_exchange",
    "custom"
  )

  required <- c(
    "rc_make_medium_scenarios(",
    scenarios,
    "medium-presets.md",
    "preset_diagnostics",
    "identical exchange bounds"
  )

  for (term in required) {
    expect_match(tutorials, term, fixed = TRUE, info = term)
  }
})

test_that("medium preset reference retains species and interpretation policy", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")

  path <- file.path(root, "docs", "medium-presets.md")
  expect_true(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  required <- c(
    "Predefined extracellular medium scenarios",
    "Accepted GEM species",
    "Concentration is not uptake flux",
    "Culture formulations are species-neutral, not physiological",
    "Use a custom medium for the actual experiment"
  )

  for (term in required) {
    expect_match(text, term, fixed = TRUE, info = term)
  }
})