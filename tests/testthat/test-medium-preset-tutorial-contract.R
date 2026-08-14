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

test_that("quick-start documents the runnable workflow and predefined medium", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- read_medium_tutorials(c(
    file.path(root, "README.md"),
    file.path(root, "docs", "tutorial-01-quick-start.md")
  ))
  required <- c(
    "rc_prepare_gem(", "rc_make_medium_scenarios(",
    'scenario = "normal_human_plasma"', "rc_run_regcompass_one_shot(",
    "target_rsq_threshold", "min_cells_per_stratum"
  )
  for (term in required) expect_match(text, term, fixed = TRUE, info = term)
})

test_that("tutorials list supported biological medium identifiers", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- read_medium_tutorials(c(
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md"),
    file.path(root, "docs", "functions.md")
  ))
  for (term in c(
    "normal_human_plasma", "mouse_plasma", "high_glucose", "low_glucose",
    "high_lactate", "low_lactate", "low_glutamine", "custom",
    "medium-presets.md"
  )) expect_match(text, term, fixed = TRUE, info = term)
})

test_that("detailed medium evidence remains in the dedicated medium document", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  medium <- read_medium_tutorials(file.path(root, "docs", "medium-presets.md"))
  for (term in c(
    "10.1016/j.cell.2017.03.023", "10.1016/j.cmet.2021.02.005",
    "10.1126/sciadv.aau7314", "10.1038/s41586-025-09898-9",
    "concentration", "flux"
  )) expect_match(medium, term, fixed = TRUE, info = term)
})

test_that("tutorials do not expose technical boundary modes as biological media", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- read_medium_tutorials(c(
    file.path(root, "README.md"),
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md")
  ))
  for (term in c("minimal", "compass_model_bounds", "permissive_all_exchange")) {
    expect_false(grepl(paste0('scenario = "', term, '"'), text, fixed = TRUE),
                 info = term)
  }
})

test_that("medium implementation no longer aliases the public function", {
  root <- medium_tutorial_root()
  if (is.null(root)) skip("Source tree is unavailable.")
  contract <- paste(
    readLines(file.path(root, "R", "published_medium_contract.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_false(grepl(
    ".rc_make_medium_scenarios_unrestricted <- rc_make_medium_scenarios",
    contract, fixed = TRUE
  ))
})
