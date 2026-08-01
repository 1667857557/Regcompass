test_that("RegCompass source requires the native Pando refit contract", {
  description <- read.dcf("DESCRIPTION")
  expect_identical(unname(description[1L, "Version"]), "2.2.5")
  expect_match(
    unname(description[1L, "Suggests"]),
    "Pando \\(>= 1\\.5\\.2\\)"
  )
  expect_identical(
    unname(description[1L, "Remotes"]),
    paste(
      "1667857557/SuperCell_Seurat_V4@agent/canonical-celltype-wnn-metacells,",
      "1667857557/Pando_regcompass"
    )
  )

  source_text <- paste(
    readLines("R/condition_grn_contract.R", warn = FALSE),
    collapse = "\n"
  )
  expect_match(source_text, 'package_version\\("1\\.5\\.2"\\)')
  expect_match(source_text, 'Config/Pando/NativeSparseABI')
  expect_match(source_text, 'cpp-eigen-direct-path-fail-fast')
  for (symbol in c(
    "_Pando_condition_product_matrix_cpp",
    "_Pando_condition_fit_multitask_path_cpp",
    "_Pando_condition_refit_path_cpp"
  )) {
    expect_match(source_text, symbol, fixed = TRUE)
  }
  expect_match(
    source_text,
    "missing registered native condition kernel",
    fixed = TRUE
  )

  tutorial <- paste(
    readLines("docs/tutorial-01-quick-start.md", warn = FALSE),
    collapse = "\n"
  )
  expect_match(tutorial, "Pando >= 1.5.2", fixed = TRUE)
  expect_match(tutorial, "native condition ABI 3", fixed = TRUE)
  expect_match(tutorial, "no R refit fallback", fixed = TRUE)
})

test_that("an installed compatible Pando exposes all native kernels", {
  skip_if_not_installed("Pando", minimum_version = "1.5.2")

  description <- utils::packageDescription("Pando")
  expect_identical(description[["Config/Pando/NativeSparseABI"]], "3")
  expect_identical(
    description[["Config/Pando/ConditionRefitBackend"]],
    "cpp-eigen-direct-path-fail-fast"
  )
  for (symbol in c(
    "_Pando_condition_product_matrix_cpp",
    "_Pando_condition_fit_multitask_path_cpp",
    "_Pando_condition_refit_path_cpp"
  )) {
    expect_true(is.loaded(symbol, PACKAGE = "Pando"), info = symbol)
  }
})
