test_that("public workflow is GRN first", {
  text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  stages <- c(
    "rc_regcompass_step_grn",
    "rc_regcompass_step_metacells",
    "rc_regcompass_step_meta_modules"
  )
  positions <- vapply(
    stages,
    function(x) regexpr(x, text, fixed = TRUE)[[1L]],
    integer(1)
  )
  expect_true(all(positions > 0L))
  expect_true(
    positions[[1L]] < positions[[2L]] &&
      positions[[2L]] < positions[[3L]]
  )
  expect_false("inference_unit" %in% names(formals(rc_run_regcompass)))
  expect_false("fragment_files" %in% names(formals(rc_run_regcompass)))
  expect_false("fragment_files" %in% names(formals(rc_regcompass_step_grn)))
  expect_false("fragment_files" %in% names(formals(rc_regcompass_step_metacells)))
  expect_true("meta_module_args" %in% names(formals(rc_run_regcompass)))
})

test_that("condition-aware Pando shares one cell-type fit across conditions", {
  implementation <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  expect_match(implementation, "Pando::infer_condition_grn", fixed = TRUE)
  expect_match(implementation, "condition_grn_fits", fixed = TRUE)
  expect_match(implementation, "cell_type = cell_type", fixed = TRUE)
  expect_match(implementation, "pando_condition_grn_fit", fixed = TRUE)
})

test_that("standard fallback uses original Pando without condition coefficients", {
  implementation <- paste(
    deparse(body(.rc_fit_standard_pando_by_cell_type)), collapse = "\n"
  )
  expect_match(implementation, "Pando::infer_grn", fixed = TRUE)
  expect_match(
    implementation, "condition_coefficients_calculated = FALSE", fixed = TRUE
  )
  expect_false(grepl("Pando::infer_condition_grn", implementation, fixed = TRUE))
})

test_that("Stage 3 uses active targets rather than target projection", {
  text <- paste(
    deparse(body(.rc_build_condition_meta_modules)), collapse = "\n"
  )
  expect_match(text, ".rc_summarize_supported_metabolic_genes", fixed = TRUE)
  expect_match(text, "rc_map_meta_module_core_reactions", fixed = TRUE)
  expect_false(grepl("rc_project_metabolic_grn", text, fixed = TRUE))
  expect_false(grepl("top_k_neighbors", text, fixed = TRUE))
})

test_that("meta-modules merge conditions within cell type only", {
  condition_modules <- list(
    condition_fit_status = data.frame(
      group_id = c("A|T", "B|T", "A|B"),
      condition = c("A", "B", "A"),
      cell_type = c("T", "T", "B"),
      status = "ok",
      stringsAsFactors = FALSE
    ),
    tf_peak_gene_condition_all = data.frame(),
    tf_peak_gene_condition = data.frame(),
    supported_metabolic_genes = data.frame(),
    core_gene_reaction = data.frame(
      group_id = c("A|T", "B|T", "A|B"),
      condition = c("A", "B", "A"),
      cell_type = c("T", "T", "B"),
      module_id = c(
        "A|T::SUPPORTED_METABOLIC_GENES",
        "B|T::SUPPORTED_METABOLIC_GENES",
        "A|B::SUPPORTED_METABOLIC_GENES"
      ),
      reaction_id = c("RT1", "RT2", "RB1"),
      is_core = TRUE,
      stringsAsFactors = FALSE
    ),
    reaction_membership = data.frame(
      group_id = c("A|T", "A|T", "B|T", "A|B"),
      condition = c("A", "A", "B", "A"),
      cell_type = c("T", "T", "T", "B"),
      module_id = c(
        "A|T::SUPPORTED_METABOLIC_GENES",
        "A|T::SUPPORTED_METABOLIC_GENES",
        "B|T::SUPPORTED_METABOLIC_GENES",
        "A|B::SUPPORTED_METABOLIC_GENES"
      ),
      reaction_id = c("RT1", "RT3", "RT2", "RB1"),
      stringsAsFactors = FALSE
    ),
    meta_module_summary = data.frame()
  )
  out <- .rc_merge_meta_modules_by_cell_type(
    condition_modules,
    celltype_col = "cell_type",
    condition_col = "condition"
  )
  expect_setequal(names(out$cell_type_catalogues), c("B", "T"))
  expect_setequal(
    out$cell_type_catalogues$T$merged_reaction_membership$reaction_id,
    c("RT1", "RT2", "RT3")
  )
  expect_identical(
    out$cell_type_catalogues$B$merged_reaction_membership$reaction_id,
    "RB1"
  )
  expect_false("RB1" %in%
    out$cell_type_catalogues$T$merged_reaction_membership$reaction_id)
  expect_false(any(c("RT1", "RT2", "RT3") %in%
    out$cell_type_catalogues$B$merged_reaction_membership$reaction_id))
  expect_identical(out$merge_scope, "cell_type")
  expect_false(out$cross_celltype_merge)
  expect_false(out$is_gem)
  expect_false(out$fastcore_applied)
})

test_that("metacells use one grouped WNN graph per cell type", {
  builder <- paste(
    deparse(body(.rc_build_grouped_wnn_membership)), collapse = "\n"
  )
  wrapper <- paste(
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
  expect_identical(.rc_condition_metacell_defaults()$gamma, 30L)
  expect_false("overwrite" %in% names(.rc_condition_metacell_defaults()))
  expect_match(
    wrapper, "celltype_grouped_joint_condition_WNN", fixed = TRUE
  )
  expect_false(grepl("supercell_stratum_col", wrapper, fixed = TRUE))
  expect_false(grepl("stratum_col", wrapper, fixed = TRUE))
  expect_false(grepl("metacell_object.rds", wrapper, fixed = TRUE))
  expect_false(grepl("rna_counts.rds", wrapper, fixed = TRUE))
})

test_that("grouped construction does not use posthoc cell-type assignment", {
  text <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)), collapse = "\n"
  )
  expect_false(grepl(
    ".rc_assign_metacell_dominant_celltype", text, fixed = TRUE
  ))
  expect_match(text, "Grouped WNN produced impure metacells", fixed = TRUE)
})

test_that("canonical metacells use in-memory assays without fragment pooling", {
  expect_false(
    "fragment_files" %in% names(formals(.rc_make_condition_celltype_metacells))
  )
  text <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)), collapse = "\n"
  )
  expect_false(grepl("fragment_manifest", text, fixed = TRUE))
})
