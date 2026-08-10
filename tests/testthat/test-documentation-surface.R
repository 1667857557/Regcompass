documentation_contract_root <- function() {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  candidates <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (!length(candidates)) return(NULL)
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

read_documentation_contract <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("tutorials keep mathematics in the mathematical specification", {
  root <- documentation_contract_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  tutorial_paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md"),
    file.path(root, "vignettes", "regcompass-workflow.Rmd")
  )
  tutorials <- paste(vapply(
    tutorial_paths, read_documentation_contract, character(1)
  ), collapse = "\n")
  forbidden <- c(
    "P^*_{", "\\widetilde P", "\\sum_jp_", "beta\\times mean(TF)"
  )
  for (term in forbidden) {
    expect_false(grepl(term, tutorials, fixed = TRUE), info = term)
  }
  expect_false(grepl("completion_time_limit = 3000", tutorials, fixed = TRUE))

  math <- read_documentation_contract(
    file.path(root, "docs", "mathematical-model.md")
  )
  required <- c(
    "beta\\times mean(TF)\\times mean(ATAC)",
    "immutable structural requirement",
    "no second parent/final closure LP pass",
    "P^*_{r,d,u,m}",
    "X^{RNA}_{g,u}"
  )
  for (term in required) {
    expect_match(math, term, fixed = TRUE, info = term)
  }
})

test_that("retired GEM helpers are absent from public documentation", {
  root <- documentation_contract_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  retired_rd <- c(
    "rc_prepare_human2_gem.Rd", "rc_prepare_mouse_gem.Rd",
    "rc_download_species_gem.Rd", "rc_bundled_gem_manifest.Rd"
  )
  expect_false(any(file.exists(file.path(root, "man", retired_rd))))

  namespace <- read_documentation_contract(file.path(root, "NAMESPACE"))
  retired_exports <- c(
    "export(rc_prepare_human2_gem)", "export(rc_prepare_mouse_gem)",
    "export(rc_download_species_gem)", "export(rc_bundled_gem_manifest)"
  )
  for (term in retired_exports) {
    expect_false(grepl(term, namespace, fixed = TRUE), info = term)
  }
})
