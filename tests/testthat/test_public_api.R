test_that("public API exposes the restartable workflow and current metacell builder", {
  expected <- c(
    "rc_prepare_gem", "rc_prepare_human2_gem", "rc_prepare_mouse_gem",
    "rc_bundled_gem_manifest", "rc_download_species_gem",
    "rc_parallel_config", "rc_make_medium_scenarios",
    "rc_make_supercell2_metacells", "rc_run_regcompass",
    "rc_run_regcompass_one_shot", "rc_regcompass_step_grn",
    "rc_regcompass_step_metacells", "rc_regcompass_step_meta_modules",
    "rc_regcompass_step_layer1", "rc_regcompass_step_layer2",
    "rc_regcompass_step_target_union", "rc_regcompass_step_results",
    "rc_test_condition_reactions", "rc_report_condition_directions",
    "rc_plot_condition_reaction", "rc_build_reaction_annotations",
    "rc_attach_reaction_annotations", "rc_select_gene_reactions",
    "rc_plot_condition_gene_reactions"
  )
  expect_setequal(getNamespaceExports("RegCompassR"), expected)
})

test_that("canonical source architecture loads current contracts", {
  description <- utils::packageDescription("RegCompassR")
  collate <- description$Collate %||% ""
  required <- c(
    "condition_metacell_cache.R", "supercell2_current_contract.R",
    "multitask_grn_bootstrap_contract.R", "multitask_grn_cv_contract.R",
    "meta_module_core_contract.R", "result_compaction.R",
    "reaction_evidence.R", "reaction_annotations.R"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = collate, fixed = TRUE)))

  removed <- c(
    "condition_pooling.R", "condition_pooling_no_sample.R",
    "metacell_no_sample_contract.R"
  )
  expect_false(any(vapply(removed, grepl, logical(1), x = collate, fixed = TRUE)))
})

test_that("canonical order is GRN then metacells then meta-modules", {
  run_body <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  positions <- vapply(
    c(
      "rc_regcompass_step_grn",
      "rc_regcompass_step_metacells",
      "rc_regcompass_step_meta_modules"
    ),
    function(x) regexpr(x, run_body, fixed = TRUE)[[1L]], integer(1)
  )
  expect_true(all(positions > 0L))
  expect_true(positions[[1L]] < positions[[2L]])
  expect_true(positions[[2L]] < positions[[3L]])
})

test_that("canonical formals separate stage settings", {
  run_formals <- names(formals(rc_run_regcompass))
  stage3_formals <- names(formals(rc_regcompass_step_meta_modules))
  expect_true("meta_module_args" %in% run_formals)
  expect_true("layer1_args" %in% run_formals)
  expect_true("meta_module_args" %in% stage3_formals)
  expect_false("layer1_args" %in% stage3_formals)
  expect_false("sample_col" %in% run_formals)
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_grn)))
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_metacells)))
  expect_false("sample_col" %in% names(formals(rc_make_supercell2_metacells)))
  expect_false("pool_col" %in% names(formals(rc_make_supercell2_metacells)))
})

test_that("current SuperCell2 defaults are explicit", {
  defaults <- formals(rc_make_supercell2_metacells)
  expect_identical(eval(defaults$strata_cols), "condition")
  expect_null(eval(defaults$label_col))
  expect_identical(eval(defaults$rna_reduction), "pca")
  expect_identical(eval(defaults$atac_reduction), "lsi")
  expect_identical(eval(defaults$rna_dims), 1:30)
  expect_identical(eval(defaults$atac_dims), 2:30)
  expect_identical(defaults$gamma, 30)
  expect_identical(defaults$seed, 12345L)
  expect_identical(defaults$min_cells_per_stratum, 100L)
  expect_identical(defaults$min_metacell_size, 20L)
  expect_identical(defaults$min_metacells_per_stratum, 2L)
})

test_that("Pando defaults and structural design remain explicit", {
  expect_null(eval(formals(rc_run_regcompass)$pfm))
  expect_null(eval(formals(rc_run_regcompass_one_shot)$pfm))
  expect_null(eval(formals(rc_regcompass_step_grn)$pfm))
  region_helper <- paste(deparse(body(.rc_default_pando_regions)), collapse = "\n")
  expect_match(region_helper, "phastConsElements20Mammals.UCSC.hg38", fixed = TRUE)
  expect_match(region_helper, "SCREEN.ccRE.UCSC.hg38", fixed = TRUE)
})

test_that("dependency pins reference merged Pando and SuperCell2 contracts", {
  description <- utils::packageDescription("RegCompassR")
  remotes <- description$Remotes %||% ""
  expect_match(
    remotes,
    "1667857557/SuperCell_Seurat_V4@c8b94949cd8a5ff7403f9f186c516f8efbac9b6f",
    fixed = TRUE
  )
  expect_match(
    remotes,
    "1667857557/Pando_regcompass@6f42c8143bec6610b001e714a51627337f6d9ba9",
    fixed = TRUE
  )
})
