test_that("workflow vignette follows the v1.7.0 public API", {
  candidates <- c(
    file.path("vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "..", "vignettes", "regcompass-workflow.Rmd")
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    skip("Source vignette is unavailable in the installed-package test context.")
  }
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")
  supported <- c(
    "rc_prepare_gem", "rc_prepare_human2_gem", "rc_prepare_mouse_gem",
    "rc_make_medium_scenarios", "rc_run_regcompass",
    "rc_run_regcompass_one_shot"
  )
  expect_true(all(vapply(
    supported,
    function(name) grepl(paste0(name, "\\("), text),
    logical(1)
  )))
  expect_match(text, "pooled across biological samples", fixed = TRUE)
  expect_match(text, "uses ATAC", fixed = TRUE)
  expect_match(text, 'fragment_files = FALSE')
  expect_match(text, 'species = "human"')
  expect_match(text, 'species = "mouse"')
  expect_match(text, 'medium_scenario = "normal_human_plasma"')
  expect_match(text, "medium_scenarios = medium")
  expect_match(text, 'inference_unit = "metacell"')
  expect_match(text, "pando_initiate_args = list")
  expect_match(text, "regions = SCREEN.ccRE.UCSC.hg38")
  expect_match(text, "microcompass\\$penalty")
  expect_match(text, "condition_contrast")
  expect_match(text, "No metabolite-neighbour", fixed = TRUE)
  expect_match(
    text,
    'meta_module_expansion = "core_subsystem_plus_kegg_reactome_master_rhea_only"',
    fixed = TRUE
  )
  expect_match(text, 'feasibility_completion = "local_fastcore_only"', fixed = TRUE)
})

test_that("tutorial and man pages exclude retired one-hop interfaces", {
  candidate_roots <- c(".", "..", file.path("..", ".."))
  root <- candidate_roots[vapply(candidate_roots, function(path) {
    dir.exists(file.path(path, "man")) &&
      dir.exists(file.path(path, "vignettes")) &&
      dir.exists(file.path(path, "docs"))
  }, logical(1))]
  if (!length(root)) {
    skip("Source documentation is unavailable in the installed-package context.")
  }
  root <- normalizePath(root[[1L]], mustWork = TRUE)
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "vignettes", "regcompass-workflow.Rmd"),
    file.path(root, "docs", "workflow.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "docs", "architecture-corrections.md"),
    file.path(root, "docs", "v1.7.0-condition-pooled-architecture.md"),
    file.path(root, "man", "rc_run_regcompass.Rd"),
    file.path(root, "man", "rc_run_regcompass_one_shot.Rd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  retired_tokens <- c(
    ".rc_meta_module_one_hop",
    "include_one_hop",
    "one_hop_max_metabolite_degree",
    "one_hop_metabolite_neighbor",
    "n_one_hop_added",
    "meta_module_one_hop = TRUE"
  )
  expect_false(any(vapply(
    retired_tokens,
    function(token) grepl(token, text, fixed = TRUE),
    logical(1)
  )))
  expect_match(
    text,
    "core_subsystem_plus_kegg_reactome_master_rhea_only",
    fixed = TRUE
  )
  expect_match(text, "local_fastcore_only", fixed = TRUE)
})
