test_that("public API exposes the restartable workflow", {
  expected <- c(
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
  expect_length(setdiff(expected, getNamespaceExports("RegCompassR")), 0L)
})

test_that("canonical source architecture has one direct routing layer", {
  description <- utils::packageDescription("RegCompassR")
  collate <- description$Collate %||% ""
  expect_false(grepl("zzz", collate, fixed = TRUE))
  required <- c(
    "stage1_input_contract.R", "condition_grn_contract.R",
    "standard_pando.R", "mixed_pando.R", "condition_pooling.R",
    "shared_tfidf.R", "step_grn_common_dictionary.R",
    "stepwise_workflow.R", "step_layer2.R", "step_results.R"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = collate, fixed = TRUE)))
  expect_false(grepl("condition_full_contract.R", collate, fixed = TRUE))
  expect_lt(
    regexpr("stage1_input_contract.R", collate, fixed = TRUE)[[1L]],
    regexpr("step_grn_common_dictionary.R", collate, fixed = TRUE)[[1L]]
  )
  expect_lt(
    regexpr("standard_pando.R", collate, fixed = TRUE)[[1L]],
    regexpr("mixed_pando.R", collate, fixed = TRUE)[[1L]]
  )
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

test_that("complete workflow exposes parallel cell-type routing", {
  run_formals <- formals(rc_run_regcompass)
  expect_identical(eval(run_formals$condition_col), "condition")
  expect_true(all(c(
    "pando_args", "metacell_args", "meta_module_args",
    "layer1_args", "layer2_args"
  ) %in% names(run_formals)))
  routing_text <- paste(
    deparse(body(.rc_fit_pando_by_celltype_route)), collapse = "\n"
  )
  worker_text <- paste(
    deparse(body(.rc_run_pando_celltype_job)), collapse = "\n"
  )
  expect_match(routing_text, "rc_parallel_lapply", fixed = TRUE)
  expect_match(worker_text, ".rc_fit_condition_grns_by_cell_type", fixed = TRUE)
  expect_match(worker_text, ".rc_fit_standard_pando_by_cell_type", fixed = TRUE)
  expect_match(worker_text, "outer_parallel", fixed = TRUE)
})

test_that("condition and standard Pando implementations are separate", {
  condition_text <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  standard_text <- paste(
    deparse(body(.rc_fit_standard_pando_by_cell_type)), collapse = "\n"
  )
  expect_match(condition_text, "Pando::infer_condition_grn", fixed = TRUE)
  expect_match(condition_text, ".rc_extract_condition_grn_contract", fixed = TRUE)
  expect_match(standard_text, "Pando::infer_grn", fixed = TRUE)
  expect_match(
    standard_text, "condition_coefficients_calculated = FALSE", fixed = TRUE
  )
  expect_false(grepl("Pando::infer_condition_grn", standard_text, fixed = TRUE))
})

test_that("parallel condition projection preserves cell-type feature spaces", {
  projection_text <- paste(
    deparse(body(.rc_condition_pando_projection)), collapse = "\n"
  )
  selector_text <- paste(
    deparse(body(.rc_condition_pando_object_for_fit)), collapse = "\n"
  )
  expect_match(
    projection_text, ".rc_condition_pando_object_for_fit", fixed = TRUE
  )
  expect_match(selector_text, "pando_grn_data_by_cell_type", fixed = TRUE)
  expect_match(selector_text, "fit$condition_cell_ids", fixed = TRUE)
})

test_that("canonical SuperCell API uses cell-type grouped WNN", {
  builder <- paste(
    deparse(body(.rc_build_grouped_wnn_membership)), collapse = "\n"
  )
  pooling <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)), collapse = "\n"
  )
  expect_match(builder, "SCimplify_by_graph_group", fixed = TRUE)
  expect_match(builder, "cell.graph.group", fixed = TRUE)
  expect_match(builder, "cell.split.condition", fixed = TRUE)
  expect_match(
    builder,
    "assay = c\\(rna_assay,[[:space:]]*atac_assay\\)"
  )
  expect_false(grepl("cell.annotation", builder, fixed = TRUE))
  expect_false(grepl("supercell_stratum_col", pooling, fixed = TRUE))
  expect_false(grepl("stratum_col", pooling, fixed = TRUE))
  expect_identical(.rc_condition_metacell_defaults()$gamma, 30L)
})

test_that("Layer 1 exposes only current controls", {
  formals_layer1 <- formals(rc_regcompass_step_layer1)
  expect_identical(eval(formals_layer1$gpr_and_method), c("min", "median", "mean"))
  expect_true("gene_half_saturation" %in% names(formals_layer1))
  expect_false(any(c(
    "tau", "projection_component", "comparison_support", "regulatory_alpha"
  ) %in% names(formals_layer1)))
})

test_that("companion repositories are unpinned", {
  description <- utils::packageDescription("RegCompassR")
  remotes <- description$Remotes %||% ""
  expect_match(remotes, "1667857557/Pando_regcompass", fixed = TRUE)
  expect_match(remotes, "1667857557/SuperCell_Seurat_V4", fixed = TRUE)
  expect_false(grepl("@", remotes, fixed = TRUE))
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