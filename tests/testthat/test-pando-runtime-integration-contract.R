test_that("condition layer overrides are validated by canonical routing", {
  skip_if_not_installed("Pando")
  routed <- .rc_route_pando_infer_args(
    list(
      rna_layer = "data",
      peak_layer = "data",
      peak_value_type = "normalized",
      tf_cor = 0.1
    ),
    condition_types = "epithelial_like",
    standard_types = character()
  )
  expect_identical(routed$condition$rna_layer, "data")
  expect_identical(routed$condition$peak_layer, "data")
  expect_identical(routed$condition$peak_value_type, "normalized")
  expect_identical(routed$condition$tf_cor, 0.1)

  expect_error(
    .rc_route_pando_infer_args(
      list(rna_layer = "counts"),
      condition_types = "epithelial_like",
      standard_types = character()
    ),
    "Unsupported override"
  )
})

test_that("standard-only routing disables condition-only layer arguments", {
  routed <- .rc_route_pando_infer_args(
    list(rna_layer = "counts", tf_cor = 0.1),
    condition_types = character(),
    standard_types = "epithelial_like"
  )
  expect_false("rna_layer" %in% names(routed$standard))
  expect_identical(routed$standard$tf_cor, 0.1)
  expect_true(any(
    routed$diagnostics$route == "standard_pando" &
      routed$diagnostics$argument == "rna_layer" &
      routed$diagnostics$action == "disabled_condition_only"
  ))
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
