test_that("public API exposes the condition-comparable restartable workflow", {
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
  expect_true(is.function(rc_build_reaction_annotations))
})

test_that("canonical source architecture has no retired compatibility layers", {
  description <- utils::packageDescription("RegCompassR")
  collate <- description$Collate %||% ""
  retired_files <- c(
    "v170_sample_balance.R", "v170_aliases.R",
    "v170_stepwise_parallel.R", "v170_tfidf.R",
    "v170_pando_reuse.R", "v170_rsq_metadata.R",
    "v170_microcompass_contract.R", "internal_apply.R",
    "pando_rsq_reliability.R", "workflow_stage_", "zzz"
  )
  expect_false(any(vapply(
    retired_files, grepl, logical(1), x = collate, fixed = TRUE
  )))
  required <- c(
    "stage_contracts.R", "shared_tfidf.R", "grn_inference.R",
    "condition_grn_contract.R", "reaction_evidence.R", "reaction_annotations.R",
    "reaction_annotation_api.R", "reaction_gene_plots.R",
    "execution_monitor.R", "bundled_gems.R", "parallel.R"
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
  source_retired <- retired_files[grepl("[.]R$", retired_files)]
  expect_false(any(file.exists(file.path(source_dir, source_retired))))
  source_text <- paste(
    unlist(lapply(list.files(source_dir, full.names = TRUE), readLines,
                  warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  expect_false(grepl("_v170", source_text, fixed = TRUE))
  retired_functions <- c(
    ".rc_build_metabolic_projection_graph",
    ".rc_mm_empty_edges",
    ".rc_mm_components",
    ".rc_signed_relation",
    "rc_project_metabolic_grn",
    ".rc_remap_projection_metadata",
    "rc_run_pando_meta_modules",
    "rc_boltzmann_minavg",
    ".rc_meta_module_one_hop"
  )
  expect_false(any(vapply(
    retired_functions, grepl, logical(1), x = source_text, fixed = TRUE
  )))
})

test_that("canonical order is Pando then metacells then meta-modules", {
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

test_that("canonical formals separate Stage 3 and Layer 1 settings", {
  run_formals <- names(formals(rc_run_regcompass))
  stage3_formals <- names(formals(rc_regcompass_step_meta_modules))
  expect_true("meta_module_args" %in% run_formals)
  expect_true("layer1_args" %in% run_formals)
  expect_true("meta_module_args" %in% stage3_formals)
  expect_false("layer1_args" %in% stage3_formals)
  expect_false("expansion_mode" %in% names(formals(rc_expand_meta_module_reactions)))
  expect_false("max_iterations" %in% names(formals(rc_expand_meta_module_reactions)))
  expect_false("tau" %in% names(formals(rc_regcompass_step_layer1)))
  expect_identical(
    eval(formals(rc_regcompass_step_layer1)$gpr_and_method),
    c("min", "median", "mean")
  )
})

test_that("Pando defaults use bundled motifs and species-specific regions", {
  expect_null(eval(formals(rc_run_regcompass)$pfm))
  expect_null(eval(formals(rc_run_regcompass_one_shot)$pfm))
  expect_null(eval(formals(rc_regcompass_step_grn)$pfm))
  expect_null(eval(formals(.rc_run_condition_single_cell_grns)$pfm))

  grn_body <- paste(
    deparse(body(.rc_run_condition_single_cell_grns)), collapse = "\n"
  )
  motif_helper <- paste(deparse(body(.rc_default_pando_motifs)), collapse = "\n")
  region_helper <- paste(deparse(body(.rc_default_pando_regions)), collapse = "\n")
  expect_match(grn_body, ".rc_default_pando_motifs", fixed = TRUE)
  expect_match(grn_body, ".rc_default_pando_regions(species)", fixed = TRUE)
  expect_match(motif_helper, 'list = "motifs"', fixed = TRUE)
  expect_match(region_helper, "phastConsElements20Mammals.UCSC.hg38", fixed = TRUE)
  expect_match(region_helper, "SCREEN.ccRE.UCSC.hg38", fixed = TRUE)
  expect_match(region_helper, 'identical(species, "mouse")', fixed = TRUE)
  expect_match(region_helper, "return(phast_cons)", fixed = TRUE)
  expect_match(region_helper, "BiocGenerics::union", fixed = TRUE)
  expect_identical(
    eval(formals(.rc_default_pando_regions)$species),
    c("human", "mouse")
  )

  grn_formals <- formals(.rc_run_condition_single_cell_grns)
  expect_identical(grn_formals$min_cells, 20L)
  expect_identical(grn_formals$padj_threshold, 0.05)
  expect_identical(grn_formals$min_abs_estimate, 0)
  expect_identical(grn_formals$min_model_rsq, 0.1)
  expect_false(isTRUE(grn_formals$require_padj))
  infer_defaults <- eval(grn_formals$pando_infer_args)
  expect_identical(infer_defaults$method, "shared_design_independent")
  expect_identical(infer_defaults$condition_mix, 1)
  expect_identical(infer_defaults$condition_weight, "equal")
  expect_true(infer_defaults$scale)
})

test_that("metacell defaults expose reductions dimensions seed and thresholds", {
  defaults <- formals(rc_make_supercell2_metacells)
  expect_identical(eval(defaults$rna_reduction), "pca")
  expect_identical(eval(defaults$atac_reduction), "lsi")
  expect_identical(eval(defaults$rna_dims), 1:30)
  expect_identical(eval(defaults$atac_dims), 2:30)
  expect_identical(defaults$gamma, 30)
  expect_identical(defaults$seed, 12345L)
  expect_identical(defaults$min_cells_per_stratum, 100)
  expect_identical(defaults$min_metacell_size, 20)
  expect_identical(defaults$min_metacells_per_stratum, 2L)

  metacell_body <- paste(
    deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n"
  )
  builder_body <- paste(
    deparse(body(.rc_build_supercell2_strata)), collapse = "\n"
  )
  expect_match(metacell_body, "gamma <- 30L", fixed = TRUE)
  expect_match(metacell_body, 'pooling_scope <- "condition_only"', fixed = TRUE)
  expect_match(metacell_body, "metacell_grouping = condition_col", fixed = TRUE)
  expect_match(metacell_body, "Sample balancing is not part", fixed = TRUE)
  expect_match(
    builder_body,
    paste0(
      "seed_i <- as.integer\\(seed\\) \\+ ",
      "match\\(key, names\\(groups\\)\\) -\\s+1L"
    )
  )
  expect_null(eval(formals(rc_regcompass_step_metacells)$sample_col))
})

test_that("Seurat stack retains required versions", {
  description <- utils::packageDescription("RegCompassR")
  imports <- description$Imports %||% ""
  expect_match(imports, "SeuratObject (>= 4.1.4)", fixed = TRUE)
  expect_match(imports, "Seurat (>= 4.4.0)", fixed = TRUE)
  expect_match(imports, "Signac (>= 1.11.0)", fixed = TRUE)
})
