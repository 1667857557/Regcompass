test_that("condition layer overrides are never silently ignored", {
  normalized <- .rc_normalize_condition_pando_layer_args(
    list(
      rna_layer = "data",
      peak_layer = "data",
      peak_value_type = "normalized",
      tf_cor = 0.1
    ),
    condition_types = "epithelial_like"
  )
  expect_identical(
    normalized$supplied,
    c("rna_layer", "peak_layer", "peak_value_type")
  )
  expect_identical(names(normalized$args), "tf_cor")

  expect_error(
    .rc_normalize_condition_pando_layer_args(
      list(rna_layer = "counts"),
      condition_types = "epithelial_like"
    ),
    "Unsupported override"
  )
})

test_that("standard-only routing leaves condition layer arguments untouched", {
  args <- list(rna_layer = "counts")
  normalized <- .rc_normalize_condition_pando_layer_args(
    args, condition_types = character()
  )
  expect_identical(normalized$args, args)
  expect_length(normalized$supplied, 0L)
})

test_that("condition fit cells agree with stored Pando metadata", {
  metadata <- data.frame(
    dataset = c("control", "treated"),
    cell_type = c("epithelial_like", "epithelial_like"),
    row.names = c("cell1", "cell2"),
    stringsAsFactors = FALSE
  )
  fit <- structure(list(
    condition_col = "dataset",
    cell_type_col = "cell_type",
    cell_type = "epithelial_like",
    condition_levels = c("control", "treated"),
    condition_cell_ids = list(control = "cell1", treated = "cell2")
  ), class = c("ConditionGRNFit", "list"))

  expect_true(.rc_validate_pando_fit_metadata_frame(
    metadata, list(fit), "dataset", "cell_type"
  ))

  wrong_column <- fit
  wrong_column$condition_col <- "condition"
  expect_error(
    .rc_validate_pando_fit_metadata_frame(
      metadata, list(wrong_column), "dataset", "cell_type"
    ),
    "do not match the RegCompass request"
  )

  wrong_assignment <- metadata
  wrong_assignment["cell2", "dataset"] <- "control"
  expect_error(
    .rc_validate_pando_fit_metadata_frame(
      wrong_assignment, list(fit), "dataset", "cell_type"
    ),
    "disagree with stored object metadata"
  )
})
