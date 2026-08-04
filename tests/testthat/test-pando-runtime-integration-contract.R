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
