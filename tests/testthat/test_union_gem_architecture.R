test_that("meta-modules merge conditions within cell type only", {
  biological <- data.frame(
    group_id = c("C1|T", "C2|T", "C1|B"),
    condition = c("C1", "C2", "C1"),
    cell_type = c("T", "T", "B"),
    module_id = c("T1", "T2", "B1"),
    reaction_id = c("RT1", "RT2", "RB1"),
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  condition_modules <- list(
    condition_fit_status = biological[, c("group_id", "condition", "cell_type")],
    tf_peak_gene_condition_all = data.frame(),
    tf_peak_gene_condition = data.frame(),
    supported_metabolic_genes = data.frame(),
    core_gene_reaction = biological,
    reaction_membership = biological,
    meta_module_summary = data.frame()
  )
  merged <- .rc_merge_meta_modules_by_cell_type(
    condition_modules, "cell_type", "condition"
  )
  expect_setequal(names(merged$cell_type_catalogues), c("T", "B"))
  expect_setequal(
    merged$cell_type_catalogues$T$merged_core_reactions$reaction_id,
    c("RT1", "RT2")
  )
  expect_identical(
    merged$cell_type_catalogues$B$merged_core_reactions$reaction_id,
    "RB1"
  )
  expect_true(all(
    merged$merged_core_reactions$cell_type ==
      c("B", "T", "T")[match(
        merged$merged_core_reactions$reaction_id,
        c("RB1", "RT1", "RT2")
      )]
  ))
  expect_identical(merged$merge_scope, "cell_type")
  expect_false(merged$cross_celltype_merge)
  expect_false(merged$is_gem)
  expect_false(merged$fastcore_applied)
})

test_that("Stage 3 contains no FASTCORE execution path", {
  construction <- paste(
    deparse(body(.rc_build_condition_meta_modules)), collapse = "\n"
  )
  stage <- paste(deparse(body(rc_regcompass_step_meta_modules)), collapse = "\n")
  expect_false(grepl(".rc_complete_celltype_medium_union_gem", construction,
                     fixed = TRUE))
  expect_false(grepl(".rc_fastcore_", construction, fixed = TRUE))
  expect_false(grepl(".rc_complete_celltype_medium_union_gem", stage,
                     fixed = TRUE))
  expect_false(grepl(".rc_fastcore_", stage, fixed = TRUE))
  expect_true(grepl("none_at_meta_module_stage", construction, fixed = TRUE))
  expect_true(grepl("merge_creates_gem = FALSE", stage, fixed = TRUE))
})

test_that("union GEM and FASTCORE scopes are cell type by medium", {
  cache_body <- paste(
    deparse(body(.rc_build_celltype_medium_union_gem_cache)), collapse = "\n"
  )
  completion_body <- paste(
    deparse(body(.rc_complete_celltype_medium_union_gem)), collapse = "\n"
  )
  engine_body <- paste(
    deparse(body(.rc_run_celltype_microcompass_engine)), collapse = "\n"
  )
  expect_true(grepl("for (cell_type in scoped$cell_types)", cache_body,
                    fixed = TRUE))
  expect_true(grepl("celltype_specific_fastcore", cache_body, fixed = TRUE))
  expect_true(grepl("cell_type = cell_type", completion_body, fixed = TRUE))
  expect_true(grepl("celltype_fastcore_support", completion_body, fixed = TRUE))
  expect_true(grepl("unit_celltype == cell_type", engine_body, fixed = TRUE))
  expect_true(grepl("shared_across_cell_types = FALSE", engine_body,
                    fixed = TRUE))
})

test_that("microCOMPASS row IDs retain cell-type structural scope", {
  parsed <- rc_parse_microcompass_row_id(
    paste0(
      "celltype=Tumor%20cell::reaction=R1::direction=forward::medium=base"
    )
  )
  expect_identical(parsed$cell_type, "Tumor cell")
  expect_identical(parsed$reaction_id, "R1")
})
