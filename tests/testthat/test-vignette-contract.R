test_that("workflow vignette follows the v1.7.0 public API", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) {
      file.path(workspace, "vignettes", "regcompass-workflow.Rmd")
    } else {
      character()
    },
    file.path("vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "vignettes", "regcompass-workflow.Rmd"),
    file.path("..", "..", "vignettes", "regcompass-workflow.Rmd")
  ))
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
  expect_match(text, "balances biological-sample cell contributions", fixed = TRUE)
  expect_match(text, "uses coefficient-weighted peak accessibility only", fixed = TRUE)
  expect_match(text, 'fragment_files = FALSE')
  expect_match(text, 'species = "human"')
  expect_match(text, 'species = "mouse"')
  expect_match(text, 'medium_scenario = "normal_human_plasma"')
  expect_match(text, "medium_scenarios = medium")
  expect_match(text, 'model_mode = "meta_module_gem"')
  expect_match(text, 'model_mode = "full_gem"')
  expect_match(text, "pando_initiate_args = list")
  expect_match(text, "regions = SCREEN.ccRE.UCSC.hg38")
  expect_match(text, "gamma = 20", fixed = TRUE)
  expect_match(text, "peak_cor = 0", fixed = TRUE)
  expect_match(text, "sample_balance = TRUE", fixed = TRUE)
  expect_match(text, "microcompass\\$penalty")
  expect_match(text, "reaction_ranking")
  expect_match(text, "single_condition_reaction_ranking", fixed = TRUE)
  expect_match(text, "condition_contrast")
  expect_match(text, "No metabolite-neighbour", fixed = TRUE)
  expect_false(grepl("strict_biological_defaults", text, fixed = TRUE))
  expect_false(grepl("inference_unit =", text, fixed = TRUE))
  expect_match(
    text,
    'meta_module_expansion = "core_subsystem_plus_kegg_reactome_master_rhea_only"',
    fixed = TRUE
  )
  expect_match(
    text,
    paste0(
      'feasibility_completion = ',
      '"local_unconstrained_fastcore_then_global_union_medium_specific_fastcore"'
    ),
    fixed = TRUE
  )
})

test_that("tutorial and man pages exclude retired interface usage", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidate_roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
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
    file.path(root, "man", "rc_run_regcompass_one_shot.Rd"),
    file.path(root, "man", "rc_regcompass_stepwise.Rd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  retired_usage_tokens <- c(
    ".rc_meta_module_one_hop(",
    "include_one_hop =",
    "one_hop_max_metabolite_degree =",
    "meta_module_one_hop = TRUE",
    "strict_biological_defaults =",
    "inference_unit = \"metacell\""
  )
  expect_false(any(vapply(
    retired_usage_tokens,
    function(token) grepl(token, text, fixed = TRUE),
    logical(1)
  )))
  expect_match(
    text,
    "core_subsystem_plus_kegg_reactome_master_rhea_only",
    fixed = TRUE
  )
  expect_match(
    text,
    "local_unconstrained_fastcore_then_global_union_medium_specific_fastcore",
    fixed = TRUE
  )
  expect_match(text, "penalty / (omega * vmax)", fixed = TRUE)
  expect_match(text, "gamma = 20", fixed = TRUE)
  expect_match(text, "peak_cor = 0", fixed = TRUE)
})
