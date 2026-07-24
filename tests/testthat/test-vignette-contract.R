test_that("workflow vignette documents the canonical one-shot API", {
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
  if (!length(candidates)) skip("Source vignette is unavailable.")
  text <- paste(readLines(candidates[[1L]], warn = FALSE), collapse = "\n")

  required <- c(
    "RegCompassR 1.8.3",
    "rc_prepare_gem",
    "rc_make_medium_scenarios",
    "scenario = \"physiologic\"",
    "Pando_regcompass.tar.gz",
    "ChromatinAssay",
    "peak_cor = 0.01",
    "gamma = 30",
    "rc_run_regcompass_one_shot(",
    "upstream_workers = 6L",
    "layer2_workers = 30L",
    "parallel = FALSE",
    "internal_threads_per_task",
    "stage_scoped_create_start_stop_release_full_gc",
    "regcompass_grn_step",
    "regcompass_layer1_step",
    "regcompass_layer2_step",
    "rc_regcompass_step_target_union(",
    "shared_kegg_reaction",
    "shared_reactome_reaction",
    "shared_master_rhea_reaction",
    "direct_kegg_reactome_master_rhea_noncore_only",
    "structural_model_reused_exactly",
    "result$version, \"1.8.3\""
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))

  forbidden <- c(
    "parallel_backend =",
    "metacell_label_col",
    "label_col =",
    "sample_balance = TRUE",
    "expansion_mode =",
    "subsystem_table =",
    "max_iterations =",
    "_v170",
    "RegCompassR.inference_unit"
  )
  expect_false(any(vapply(forbidden, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("tutorials preserve one-shot, true stepwise, targeted, and differential modes", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) dir.exists(file.path(path, "docs")),
    logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")
  root <- normalizePath(roots[[1L]], mustWork = TRUE)
  paths <- file.path(root, "docs", c(
    "tutorial-01-quick-start.md",
    "tutorial-02-stepwise-audit.md",
    "tutorial-03-advanced-restart.md",
    "tutorial-04-targeted-reaction-remapping.md",
    "tutorial-05-condition-differential-analysis.md",
    "target-union-scoring.md",
    "portable-execution.md"
  ))
  expect_true(all(file.exists(paths)))
  text <- lapply(paths, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  })

  expect_match(text[[1L]], "Tutorial Level 1", fixed = TRUE)
  expect_match(text[[1L]], "rc_run_regcompass_one_shot(", fixed = TRUE)
  expect_match(text[[1L]], "gamma = 30", fixed = TRUE)
  expect_match(text[[1L]], "upstream_workers = 6L", fixed = TRUE)
  expect_match(text[[1L]], "layer2_workers = 30L", fixed = TRUE)

  expect_match(text[[2L]], "Tutorial Level 2", fixed = TRUE)
  expect_false(grepl("rc_run_regcompass_one_shot(", text[[2L]], fixed = TRUE))
  for (fun in c(
    "rc_regcompass_step_grn(",
    "rc_regcompass_step_metacells(",
    "rc_regcompass_step_meta_modules(",
    "rc_regcompass_step_layer1(",
    "rc_regcompass_step_layer2(",
    "rc_regcompass_step_results("
  )) {
    expect_match(text[[2L]], fun, fixed = TRUE)
  }
  expect_match(text[[2L]], "gamma = 30", fixed = TRUE)
  expect_match(text[[2L]], "BPPARAM = upstream_bp", fixed = TRUE)
  expect_match(text[[2L]], "BPPARAM = layer2_bp", fixed = TRUE)

  expect_match(text[[3L]], "Tutorial Level 3", fixed = TRUE)
  expect_match(text[[3L]], "Earliest stage to rerun", fixed = TRUE)
  expect_match(text[[3L]], "Serial troubleshooting", fixed = TRUE)

  expect_match(text[[4L]], "Tutorial Level 4", fixed = TRUE)
  expect_match(text[[4L]], "rc_regcompass_step_target_union(", fixed = TRUE)
  expect_match(text[[4L]], "core_genes", fixed = TRUE)
  expect_match(text[[4L]], "core_reaction_ids", fixed = TRUE)
  expect_match(text[[4L]], "shared_kegg_reaction", fixed = TRUE)
  expect_match(
    text[[4L]],
    "direct_kegg_reactome_master_rhea_noncore_only",
    fixed = TRUE
  )
  expect_match(text[[4L]], "structural_model_reused_exactly", fixed = TRUE)

  expect_match(text[[5L]], "Tutorial Level 5", fixed = TRUE)
  expect_match(text[[5L]], "rc_test_condition_reactions(", fixed = TRUE)
  expect_match(text[[5L]], "rc_plot_condition_reaction(", fixed = TRUE)
  expect_match(text[[5L]], "Kruskal-Wallis", fixed = TRUE)
  expect_match(text[[5L]], "Wilcoxon", fixed = TRUE)
  expect_match(text[[5L]], "metacell_within_dataset", fixed = TRUE)

  expect_match(text[[6L]], "No same-subsystem expansion", fixed = TRUE)
  expect_match(text[[6L]], "anchor_core_reaction_id", fixed = TRUE)
  expect_match(text[[6L]], "source_model_md5", fixed = TRUE)
  expect_match(
    text[[7L]],
    "One outer worker equals one single-thread task",
    fixed = TRUE
  )
  expect_match(text[[7L]], "gc(full = TRUE)", fixed = TRUE)

  combined <- paste(unlist(text), collapse = "\n")
  expect_match(combined, "peak_cor = 0.01", fixed = TRUE)
  expect_match(combined, "gamma = 30", fixed = TRUE)
  expect_match(combined, "OMP_NUM_THREADS=1", fixed = TRUE)
  expect_match(combined, "internal_threads_per_task", fixed = TRUE)
  expect_false(grepl("gamma = 75", combined, fixed = TRUE))
  expect_false(grepl("gamma = 100", combined, fixed = TRUE))
  expect_false(grepl("expansion_mode =", combined, fixed = TRUE))
  expect_false(grepl("subsystem_table =", combined, fixed = TRUE))
  expect_false(grepl("max_iterations =", combined, fixed = TRUE))
  expect_false(grepl("_v170", combined, fixed = TRUE))
})

test_that("README and API index expose current public workflow", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) {
      file.exists(file.path(path, "README.md")) &&
        file.exists(file.path(path, "docs", "functions.md"))
    },
    logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")
  root <- normalizePath(roots[[1L]], mustWork = TRUE)
  paths <- c(
    file.path(root, "README.md"),
    file.path(root, "docs", "functions.md"),
    file.path(root, "man", "rc_regcompass_stepwise.Rd"),
    file.path(root, "man", "rc_regcompass_step_target_union.Rd"),
    file.path(root, "man", "rc_run_regcompass.Rd")
  )
  expect_true(all(file.exists(paths)))
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")

  required <- c(
    "RegCompassR 1.8.3",
    "rc_run_regcompass_one_shot",
    "rc_regcompass_step_target_union",
    "upstream_workers",
    "layer2_workers",
    "GEM fingerprint",
    "ordered metacell IDs",
    "KEGG",
    "Reactome",
    "master-Rhea",
    "Same-subsystem",
    "medium-presets.md"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))

  forbidden <- c(
    "parallel_backend = c(",
    "v170_microcompass_contract",
    "internal_apply",
    "metacell_label_col",
    "sample_balance_seed",
    "expansion_mode =",
    "subsystem_table =",
    "max_iterations ="
  )
  expect_false(any(vapply(forbidden, grepl, logical(1), x = text, fixed = TRUE)))
})
