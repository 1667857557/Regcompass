test_that("parallel condition jobs cannot share paired cells", {
  results <- list(
    list(paired_cell_metadata = data.frame(
      cell_id = c("cell1", "cell2"), stringsAsFactors = FALSE
    )),
    list(paired_cell_metadata = data.frame(
      cell_id = c("cell2", "cell3"), stringsAsFactors = FALSE
    ))
  )
  expect_error(
    .rc_validate_pando_result_cell_partition(results),
    "more than one Pando route"
  )
})

test_that("canonical condition merge keeps one Pando object per cell type", {
  fake_object <- structure(list(), class = "GRNData")
  make_result <- function(cell_type, cell_id) {
    list(
      pando_grn_data = fake_object,
      condition_grn_fits = list(structure(
        list(cell_type = cell_type),
        class = c("ConditionGRNFit", "list")
      )),
      paired_cell_metadata = data.frame(
        cell_id = cell_id, stringsAsFactors = FALSE
      ),
      target_metabolic_genes = paste0("GENE_", cell_type),
      pando_execution_summary = list(
        fit_engine = "fixed_dictionary_glm",
        targets_total = 1L,
        targets_failed = 0L
      )
    )
  }

  answer <- .rc_merge_condition_job_results(list(
    make_result("epithelial_like", "cell1"),
    make_result("stem_like", "cell2")
  ))

  expect_setequal(
    names(answer$pando_grn_data_by_cell_type),
    c("epithelial_like", "stem_like")
  )
  expect_null(answer$pando_grn_data)
  expect_true(answer$pando_object_scope$preserves_cell_type_peak_space)
  expect_false(answer$pando_object_scope$combined_grndata)
})

test_that("canonical condition merge rejects duplicated cell-type results", {
  fake_object <- structure(list(), class = "GRNData")
  make_result <- function(cell_id) {
    list(
      pando_grn_data = fake_object,
      condition_grn_fits = list(structure(
        list(cell_type = "epithelial_like"),
        class = c("ConditionGRNFit", "list")
      )),
      paired_cell_metadata = data.frame(
        cell_id = cell_id, stringsAsFactors = FALSE
      ),
      target_metabolic_genes = "GENE",
      pando_execution_summary = list(
        fit_engine = "fixed_dictionary_glm",
        targets_total = 1L,
        targets_failed = 0L
      )
    )
  }
  expect_error(
    .rc_merge_condition_job_results(list(
      make_result("cell1"), make_result("cell2")
    )),
    "duplicated"
  )
})
