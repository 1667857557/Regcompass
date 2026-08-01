test_that("RegCompass source requires the Pando hybrid ABI 5 contract", {
  description <- read.dcf("DESCRIPTION")
  expect_identical(unname(description[1L, "Version"]), "2.2.5")
  expect_match(
    unname(description[1L, "Suggests"]),
    "Pando \\(>= 1\\.6\\.1\\)"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoNativeSparseABI"]),
    "5"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionRefitBackend"]),
    "cpp-eigen-direct-path-fail-fast"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionTargetEngineBackend"]),
    "cpp-eigen-fused-hybrid-gram-nested-cv-path-refit-validation-stats-fail-fast"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionInnerCVBackend"]),
    "cpp-eigen-hybrid-gram-sufficient-statistics-fail-fast"
  )
  expect_identical(
    unname(description[1L, "Remotes"]),
    paste(
      "1667857557/SuperCell_Seurat_V4@agent/canonical-celltype-wnn-metacells,",
      "1667857557/Pando_regcompass"
    )
  )

  runtime_text <- paste(
    readLines("R/condition_grn_runtime.R", warn = FALSE),
    collapse = "\n"
  )
  bridge_text <- paste(
    readLines("R/condition_grn_contract.R", warn = FALSE),
    collapse = "\n"
  )
  expect_match(runtime_text, 'package_version\\("1\\.6\\.1"\\)')
  expect_match(runtime_text, 'Config/Pando/NativeSparseABI', fixed = TRUE)
  expect_match(
    runtime_text,
    "cpp-eigen-hybrid-gram-sufficient-statistics-fail-fast",
    fixed = TRUE
  )
  expect_match(
    runtime_text,
    "cpp-eigen-fused-hybrid-gram-nested-cv-path-refit-validation-stats-fail-fast",
    fixed = TRUE
  )
  expect_match(runtime_text, "cpp-eigen-direct-path-fail-fast", fixed = TRUE)
  for (symbol in c(
    "_Pando_condition_product_matrix_cpp",
    "_Pando_condition_fit_multitask_path_cpp",
    "_Pando_condition_refit_path_cpp",
    "_Pando_condition_fit_target_engine_cpp"
  )) {
    expect_match(runtime_text, symbol, fixed = TRUE)
  }
  expect_match(
    runtime_text,
    "missing registered native condition kernel",
    fixed = TRUE
  )
  expect_match(
    runtime_text,
    "No R fallback is permitted",
    fixed = TRUE
  )
  expect_match(
    bridge_text,
    ".rc_require_pando_hybrid_runtime()",
    fixed = TRUE
  )

  tutorial <- paste(
    readLines("docs/tutorial-01-quick-start.md", warn = FALSE),
    collapse = "\n"
  )
  expect_match(tutorial, "Pando >= 1.6.1", fixed = TRUE)
  expect_match(tutorial, "native condition ABI 5", fixed = TRUE)
  expect_match(tutorial, "hybrid Gram/sparse", fixed = TRUE)
  expect_match(tutorial, "no R fallback", fixed = TRUE)
})

test_that("an installed compatible Pando exposes the hybrid native runtime", {
  skip_if_not_installed("Pando", minimum_version = "1.6.1")

  description <- utils::packageDescription("Pando")
  expect_identical(description[["Config/Pando/NativeSparseABI"]], "5")
  expect_identical(
    description[["Config/Pando/ConditionRefitBackend"]],
    "cpp-eigen-direct-path-fail-fast"
  )
  expect_identical(
    description[["Config/Pando/ConditionTargetEngineBackend"]],
    "cpp-eigen-fused-hybrid-gram-nested-cv-path-refit-validation-stats-fail-fast"
  )
  expect_identical(
    description[["Config/Pando/ConditionInnerCVBackend"]],
    "cpp-eigen-hybrid-gram-sufficient-statistics-fail-fast"
  )
  for (symbol in c(
    "_Pando_condition_product_matrix_cpp",
    "_Pando_condition_fit_multitask_path_cpp",
    "_Pando_condition_refit_path_cpp",
    "_Pando_condition_fit_target_engine_cpp"
  )) {
    expect_true(is.loaded(symbol, PACKAGE = "Pando"), info = symbol)
  }
  expect_silent(RegCompassR:::.rc_require_pando_hybrid_runtime())
})
