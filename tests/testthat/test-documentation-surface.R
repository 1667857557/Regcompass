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

test_that("interface documentation keeps mathematics in one specification", {
  root <- documentation_contract_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  interface_paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "tutorial-01-quick-start.md"),
    file.path(root, "docs", "tutorial-02-stepwise-audit.md"),
    file.path(root, "docs", "tutorial-04-post-analysis.md"),
    file.path(root, "docs", "layer2-corda.md"),
    file.path(root, "docs", "layer2-model-builders.md"),
    file.path(root, "vignettes", "regcompass-workflow.Rmd")
  )
  interface_text <- paste(vapply(
    interface_paths, read_documentation_contract, character(1)
  ), collapse = "\n")
  forbidden <- c(
    "P^*_{", "\\widetilde P", "\\sum_jp_", "beta\\times mean(TF)",
    "v^{max}", "X^{RNA}_{"
  )
  for (term in forbidden) {
    expect_false(grepl(term, interface_text, fixed = TRUE), info = term)
  }
  expect_false(grepl("completion_time_limit = 3000", interface_text, fixed = TRUE))

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

test_that("function index follows the supported public API", {
  root <- documentation_contract_root()
  if (is.null(root)) skip("Source documentation is unavailable.")
  text <- read_documentation_contract(file.path(root, "docs", "functions.md"))
  expected <- c(
    "rc_prepare_gem", "rc_parallel_config", "rc_make_medium_scenarios",
    "rc_run_regcompass", "rc_run_regcompass_one_shot",
    "rc_regcompass_step_grn", "rc_regcompass_step_metacells",
    "rc_regcompass_step_meta_modules", "rc_regcompass_step_layer1",
    "rc_regcompass_step_layer2", "rc_regcompass_step_target_union",
    "rc_regcompass_step_results", "rc_test_condition_reactions",
    "rc_plot_condition_reaction", "rc_build_reaction_annotations",
    "rc_attach_reaction_annotations", "rc_select_gene_reactions",
    "rc_plot_condition_gene_reactions", "plot_top_celltype_reaction_rank"
  )
  for (name in expected) {
    expect_match(text, paste0("`", name, "()"), fixed = TRUE, info = name)
  }
  retired <- c(
    "`rc_prepare_human2_gem()`", "`rc_prepare_mouse_gem()`",
    "`rc_download_species_gem()`", "`rc_bundled_gem_manifest()`"
  )
  for (name in retired) expect_false(grepl(name, text, fixed = TRUE), info = name)
})

test_that("retired GEM helpers remain undocumented and cannot be re-exported by roxygen", {
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

  sources <- paste(vapply(
    c(file.path(root, "R", "humangem.R"),
      file.path(root, "R", "bundled_gems.R"),
      file.path(root, "R", "gem_compat.R")),
    read_documentation_contract, character(1)
  ), collapse = "\n")
  for (name in c(
    "rc_prepare_human2_gem", "rc_prepare_mouse_gem",
    "rc_download_species_gem", "rc_bundled_gem_manifest"
  )) {
    expect_false(grepl(
      paste0("#' @export\n", name, " <- function"),
      sources, fixed = TRUE
    ), info = name)
  }
})
