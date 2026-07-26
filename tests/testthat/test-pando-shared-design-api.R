test_that("Pando exposes the shared GRN design contract", {
  skip_if_not_installed("Pando")
  expect_true(utils::packageVersion("Pando") >= "1.1.2")
  exports <- getNamespaceExports("Pando")
  expect_true("prepare_grn_design" %in% exports)
  expect_true("validate_grn_design" %in% exports)
  expect_true(exists(
    "prepare_grn_design.GRNData",
    envir = asNamespace("Pando"),
    inherits = FALSE
  ))
})
