test_that("public API exposes the GRN-first restartable workflow", {
  expect_setequal(
    getNamespaceExports("RegCompassR"),
    c(
      "rc_prepare_gem", "rc_prepare_human2_gem", "rc_prepare_mouse_gem",
      "rc_bundled_gem_manifest", "rc_download_species_gem",
      "rc_parallel_config", "rc_make_medium_scenarios", "rc_run_regcompass",
      "rc_run_regcompass_one_shot", "rc_regcompass_step_grn",
      "rc_regcompass_step_metacells", "rc_regcompass_step_meta_modules",
      "rc_regcompass_step_layer1", "rc_regcompass_step_layer2",
      "rc_regcompass_step_target_union", "rc_regcompass_step_results",
      "rc_test_condition_reactions", "rc_plot_condition_reaction",
      "rc_build_reaction_annotations", "rc_attach_reaction_annotations",
      "rc_select_gene_reactions", "rc_plot_condition_gene_reactions"
    )
  )
})

test_that("canonical source architecture has no retired compatibility layers", {
  description <- utils::packageDescription("RegCompassR")
  collate <- description$Collate %||% ""
  retired <- c(
    "v170_sample_balance.R", "v170_aliases.R",
    "v170_stepwise_parallel.R", "v170_tfidf.R",
    "v170_pando_reuse.R", "v170_rsq_metadata.R",
    "v170_microcompass_contract.R", "internal_apply.R",
    "pando_rsq_reliability.R", "workflow_stage_", "zzz"
  )
  expect_false(any(vapply(retired, grepl, logical(1), x = collate, fixed = TRUE)))
  required <- c(
    "stage_contracts.R", "shared_tfidf.R", "grn_inference.R",
    "regulatory_modifier.R", "reaction_annotations.R", "reaction_evidence.R",
    "reaction_gene_plots.R", "execution_monitor.R", "bundled_gems.R",
    "parallel.R"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = collate, fixed = TRUE)))

  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) file.path(workspace, "R") else character(),
    "R", file.path("..", "R"), file.path("..", "..", "R")
  ))
  candidates <- candidates[dir.exists(candidates)]
  if (!length(candidates)) skip("Source R files are unavailable.")
  source_dir <- normalizePath(candidates[[1L]], mustWork = TRUE)
  source_retired <- retired[grepl("[.]R$", retired)]
  expect_false(any(file.exists(file.path(source_dir, source_retired))))
  source_text <- paste(
    unlist(lapply(list.files(source_dir, full.names = TRUE), readLines,
                  warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  expect_false(grepl("_v170", source_text, fixed = TRUE))
})

test_that("canonical order is GRN then metacells then meta-modules", {
  run_body <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  positions <- vapply(
    c("rc_regcompass_step_grn", "rc_regcompass_step_metacells", "rc_regcompass_step_meta_modules"),
    function(x) regexpr(x, run_body, fixed = TRUE)[[1L]], integer(1)
  )
  expect_true(all(positions > 0L))
  expect_true(positions[[1L]] < positions[[2L]])
  expect_true(positions[[2L]] < positions[[3L]])
})

test_that("GRN and metacell defaults match the canonical design", {
  grn_body <- paste(deparse(body(.rc_run_condition_single_cell_grns)), collapse = "\n")
  metacell_body <- paste(deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n")
  expect_match(grn_body, "peak_cor = 0.01", fixed = TRUE)
  expect_match(grn_body, "condition_col, celltype_col", fixed = TRUE)
  expect_match(metacell_body, "gamma <- 75L", fixed = TRUE)
  expect_match(metacell_body, 'pooling_scope <- "condition_only"', fixed = TRUE)
  expect_match(metacell_body, "metacell_grouping = condition_col", fixed = TRUE)
  expect_match(metacell_body, "Sample balancing is not part", fixed = TRUE)
  expect_null(eval(formals(rc_regcompass_step_metacells)$sample_col))
})

test_that("Seurat stack retains required versions", {
  description <- utils::packageDescription("RegCompassR")
  imports <- description$Imports %||% ""
  expect_match(imports, "SeuratObject (>= 4.1.4)", fixed = TRUE)
  expect_match(imports, "Seurat (>= 4.4.0)", fixed = TRUE)
  expect_match(imports, "Signac (>= 1.11.0)", fixed = TRUE)
})
