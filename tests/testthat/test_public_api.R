test_that("public API exposes the restartable workflow", {
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
      "rc_report_condition_directions",
      "rc_build_reaction_annotations", "rc_attach_reaction_annotations",
      "rc_select_gene_reactions", "rc_plot_condition_gene_reactions"
    )
  )
})

test_that("canonical source architecture has no runtime override layers", {
  description <- utils::packageDescription("RegCompassR")
  collate <- description$Collate %||% ""
  retired <- c(
    "zzz00_absolute_pando_contract.R",
    "zzz01_fixed_gamma_metacells.R",
    "zzz02_layer1_policy.R",
    "zzz03_compass_gpr_penalty.R",
    "zzz04_canonical_pando_fit_schema.R",
    "workflow_stage_", "v170_"
  )
  expect_false(any(vapply(retired, grepl, logical(1), x = collate, fixed = TRUE)))
  required <- c(
    "condition_grn_contract.R", "standard_pando.R",
    "condition_pooling.R", "metacell_object_merge.R",
    "step_layer1_oof.R", "penalty.R"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = collate, fixed = TRUE)))
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

test_that("complete workflow exposes automatic condition routing", {
  run_formals <- formals(rc_run_regcompass)
  expect_identical(eval(run_formals$condition_col), "condition")
  expect_true(all(c(
    "pando_args", "metacell_args", "meta_module_args",
    "layer1_args", "layer2_args"
  ) %in% names(run_formals)))

  step_formals <- formals(rc_regcompass_step_grn)
  expect_identical(eval(step_formals$condition_col), "condition")
  expect_true(all(c("parallel", "BPPARAM") %in% names(step_formals)))
})

test_that("condition and standard Pando implementations are separate", {
  condition_text <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  standard_text <- paste(
    deparse(body(.rc_fit_standard_pando_by_cell_type)), collapse = "\n"
  )
  expect_match(condition_text, "Pando::infer_condition_grn", fixed = TRUE)
  expect_match(condition_text, "pando_condition_grn_fit", fixed = TRUE)
  expect_match(standard_text, "Pando::infer_grn", fixed = TRUE)
  expect_match(
    standard_text, "condition_coefficients_calculated = FALSE", fixed = TRUE
  )
  expect_false(grepl("Pando::infer_condition_grn", standard_text, fixed = TRUE))
})

test_that("native SuperCell API uses cell-type graph groups", {
  native <- paste(
    deparse(body(.rc_native_supercell_membership)), collapse = "\n"
  )
  pooling <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)), collapse = "\n"
  )
  expect_match(
    native, "SCimplify_by_graph_group_from_embedding", fixed = TRUE
  )
  expect_match(native, "cell.graph.group", fixed = TRUE)
  expect_match(native, "cell.split.condition", fixed = TRUE)
  expect_match(native, ".rc_scale_embedding_block_by_group", fixed = TRUE)
  expect_false(grepl("cell.annotation", native, fixed = TRUE))
  expect_false(grepl("supercell_stratum_col", pooling, fixed = TRUE))
  expect_false(grepl("stratum_col", pooling, fixed = TRUE))
})

test_that("Layer 1 retains canonical controls", {
  formals_layer1 <- formals(rc_regcompass_step_layer1)
  expect_identical(eval(formals_layer1$gpr_and_method), c("min", "median", "mean"))
  expect_identical(formals_layer1$regulatory_alpha, 1)
  expect_false("tau" %in% names(formals_layer1))
})

test_that("Pando defaults use bundled human inputs and guard mouse regions", {
  expect_null(eval(formals(rc_run_regcompass)$pfm))
  expect_null(eval(formals(rc_regcompass_step_grn)$pfm))
  expect_null(eval(formals(.rc_fit_condition_grns_by_cell_type)$pfm))
  expect_null(eval(formals(.rc_fit_standard_pando_by_cell_type)$pfm))
  motif_helper <- paste(deparse(body(.rc_default_pando_motifs)), collapse = "\n")
  region_guard <- paste(deparse(body(.rc_default_pando_regions)), collapse = "\n")
  expect_match(motif_helper, 'list = "motifs"', fixed = TRUE)
  expect_match(region_guard, "phastConsElements20Mammals.UCSC.hg38", fixed = TRUE)
  expect_match(region_guard, "SCREEN.ccRE.UCSC.hg38", fixed = TRUE)
  expect_match(region_guard, 'identical(species, "mouse")', fixed = TRUE)
})

test_that("Seurat stack retains required versions", {
  description <- utils::packageDescription("RegCompassR")
  imports <- description$Imports %||% ""
  expect_match(imports, "SeuratObject (>= 4.1.4)", fixed = TRUE)
  expect_match(imports, "Seurat (>= 4.4.0)", fixed = TRUE)
  expect_match(imports, "Signac (>= 1.11.0)", fixed = TRUE)
})
