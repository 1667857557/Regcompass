test_that("Stage 2 inherits omitted assay names from Stage 1", {
  counts <- matrix(
    c(1, 0, 2, 3), nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2"))
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object[["peaks"]] <- SeuratObject::CreateAssayObject(counts = counts)

  contract <- list(
    retained_cells = colnames(object),
    retained_cell_types = "T",
    condition_pando_cell_types = "T",
    standard_pando_cell_types = character(),
    diagnostics = data.frame(),
    analysis_mode = "condition_grn",
    condition_levels = c("baseline", "treated")
  )
  grn <- list(
    params = list(
      requested_condition_col = "timepoint",
      condition_col = "timepoint",
      celltype_col = "Cell type",
      rna_assay = "RNA",
      atac_assay = "peaks",
      analysis_mode = "condition_grn",
      fallback_reason = NA_character_
    )
  )

  testthat::local_mocked_bindings(
    .rc_step_monitor_start = function(...) list(),
    .rc_step_monitor_fail = function(...) invisible(NULL),
    .rc_stage_worker_config = function(...) list(worker_limit = 1L),
    .rc_fragment_input_enabled = function(...) FALSE,
    .rc_resolve_fragment_aggregation_args = function(...) list(),
    .rc_validate_stage1_cell_set = function(...) contract,
    .rc_subset_to_stage1_cell_set = function(...) {
      stop("passed Stage 1 assay contract", call. = FALSE)
    },
    .package = "RegCompassR"
  )

  expect_error(
    rc_regcompass_step_metacells(
      object = object,
      grn = grn,
      outdir = tempfile("stage2-assay-inherit-"),
      condition_col = "timepoint",
      celltype_col = "Cell type",
      workers = 1L,
      progress = FALSE
    ),
    "passed Stage 1 assay contract",
    fixed = TRUE
  )
})

test_that("Stage 2 still rejects an explicitly different assay", {
  counts <- matrix(
    c(1, 0, 2, 3), nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2"))
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object[["peaks"]] <- SeuratObject::CreateAssayObject(counts = counts)

  contract <- list(
    retained_cells = colnames(object),
    retained_cell_types = "T",
    condition_pando_cell_types = "T",
    standard_pando_cell_types = character(),
    diagnostics = data.frame(),
    analysis_mode = "condition_grn",
    condition_levels = c("baseline", "treated")
  )
  grn <- list(
    params = list(
      requested_condition_col = "timepoint",
      condition_col = "timepoint",
      celltype_col = "Cell type",
      rna_assay = "RNA",
      atac_assay = "peaks",
      analysis_mode = "condition_grn",
      fallback_reason = NA_character_
    )
  )

  testthat::local_mocked_bindings(
    .rc_step_monitor_start = function(...) list(),
    .rc_step_monitor_fail = function(...) invisible(NULL),
    .rc_stage_worker_config = function(...) list(worker_limit = 1L),
    .rc_fragment_input_enabled = function(...) FALSE,
    .rc_resolve_fragment_aggregation_args = function(...) list(),
    .rc_validate_stage1_cell_set = function(...) contract,
    .package = "RegCompassR"
  )

  expect_error(
    rc_regcompass_step_metacells(
      object = object,
      grn = grn,
      outdir = tempfile("stage2-assay-mismatch-"),
      condition_col = "timepoint",
      celltype_col = "Cell type",
      atac_assay = "ATAC",
      workers = 1L,
      progress = FALSE
    ),
    "Stage 2 assay names differ from Stage 1",
    fixed = TRUE
  )
})
