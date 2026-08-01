test_that("RegCompass source requires the worker-safe Pando ABI 5 contract", {
  description <- read.dcf("DESCRIPTION")
  expect_identical(unname(description[1L, "Version"]), "2.2.6")
  expect_match(
    unname(description[1L, "Suggests"]),
    "Pando \\(>= 1\\.6\\.2\\)"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoNativeSparseABI"]),
    "5"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoNativeCallBinding"]),
    "registered-symbol-lookup-worker-safe-v1"
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
  expect_match(
    unname(description[1L, "Collate"]),
    "condition_grn_runtime_guard.R",
    fixed = TRUE
  )

  runtime_text <- paste(
    readLines("R/condition_grn_runtime.R", warn = FALSE),
    collapse = "\n"
  )
  guard_text <- paste(
    readLines("R/condition_grn_runtime_guard.R", warn = FALSE),
    collapse = "\n"
  )
  expect_match(runtime_text, 'package_version\\("1\\.6\\.2"\\)')
  expect_match(runtime_text, 'Config/Pando/NativeSparseABI', fixed = TRUE)
  expect_match(runtime_text, 'Config/Pando/NativeCallBinding', fixed = TRUE)
  expect_match(
    runtime_text,
    "registered-symbol-lookup-worker-safe-v1",
    fixed = TRUE
  )
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
  expect_match(runtime_text, "getNativeSymbolInfo", fixed = TRUE)
  expect_match(runtime_text, "BiocParallel::bplapply", fixed = TRUE)
  expect_match(runtime_text, ".pando_registered_call", fixed = TRUE)
  expect_match(
    runtime_text,
    "No R fallback is permitted",
    fixed = TRUE
  )
  expect_match(
    guard_text,
    ".rc_require_pando_hybrid_runtime(BPPARAM = BPPARAM)",
    fixed = TRUE
  )
  expect_match(
    guard_text,
    ".rc_fit_condition_grns_by_cell_type_unchecked",
    fixed = TRUE
  )
})

test_that("an installed compatible Pando exposes the worker-safe runtime", {
  skip_if_not_installed("Pando", minimum_version = "1.6.2")

  description <- utils::packageDescription("Pando")
  expect_identical(description[["Config/Pando/NativeSparseABI"]], "5")
  expect_identical(
    description[["Config/Pando/NativeCallBinding"]],
    "registered-symbol-lookup-worker-safe-v1"
  )
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
  namespace <- asNamespace("Pando")
  for (wrapper in c(
    ".condition_product_matrix_cpp",
    ".condition_fit_multitask_path_cpp",
    ".condition_refit_path_cpp",
    ".condition_fit_target_engine_cpp"
  )) {
    value <- get(wrapper, namespace, inherits = FALSE)
    expect_match(
      paste(deparse(body(value)), collapse = "\n"),
      ".pando_registered_call",
      fixed = TRUE,
      info = wrapper
    )
  }
  for (symbol in c(
    "_Pando_condition_product_matrix_cpp",
    "_Pando_condition_fit_multitask_path_cpp",
    "_Pando_condition_refit_path_cpp",
    "_Pando_condition_fit_target_engine_cpp"
  )) {
    info <- getNativeSymbolInfo(
      symbol,
      PACKAGE = "Pando",
      withRegistrationInfo = TRUE
    )
    expect_true(is.list(info) && !is.null(info$address), info = symbol)
  }
  expect_silent(RegCompassR:::.rc_require_pando_hybrid_runtime())
})
