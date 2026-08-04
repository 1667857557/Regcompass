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

test_that("parallel output maps fits to the same unique cell types", {
  fake_object <- structure(list(), class = "GRNData")
  answer <- list(
    condition_grn_fits = list(
      epithelial_like = structure(
        list(cell_type = "epithelial_like"),
        class = c("ConditionGRNFit", "list")
      )
    ),
    pando_grn_data_by_cell_type = list(
      epithelial_like = fake_object
    ),
    paired_cell_metadata = data.frame(
      cell_id = "cell1", stringsAsFactors = FALSE
    )
  )
  expect_true(.rc_validate_condition_job_merge_output(answer))

  names(answer$pando_grn_data_by_cell_type) <- "wrong_type"
  expect_error(
    .rc_validate_condition_job_merge_output(answer),
    "same unique cell-type set"
  )
})
