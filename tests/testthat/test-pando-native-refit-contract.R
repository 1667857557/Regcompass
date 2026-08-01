test_that("RegCompass source requires the memory-bounded Pando ABI 6 contract", {
  package_root <- testthat::test_path("..", "..")
  description <- read.dcf(file.path(package_root, "DESCRIPTION"))
  expect_identical(unname(description[1L, "Version"]), "2.2.7")
  expect_match(
    unname(description[1L, "Suggests"]),
    "Pando \\(>= 1\\.6\\.3\\)"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoNativeSparseABI"]),
    "6"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoNativeCallBinding"]),
    "registered-symbol-lookup-worker-safe-v1"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionRefitBackend"]),
    "dense-direct-or-matrix-free-schur-pcg-v1"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionTargetEngineBackend"]),
    "cpp-eigen-memory-bounded-hybrid-target-v1"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionInnerCVBackend"]),
    "exact-refit-validation-sparse-residual-v1"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionMemoryContract"]),
    "no-full-p2-on-high-p-path-v1"
  )
  expect_identical(
    unname(description[1L, "Remotes"]),
    paste(
      "1667857557/SuperCell_Seurat_V4@agent/canonical-celltype-wnn-metacells,",
      "1667857557/Pando_regcompass@agent/high-p-memory-bounded-engine"
    )
  )
  expect_match(
    unname(description[1L, "Collate"]),
    "condition_grn_runtime_guard.R",
    fixed = TRUE
  )

  runtime_text <- paste(
    readLines(
      file.path(package_root, "R", "condition_grn_runtime.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  guard_text <- paste(
    readLines(
      file.path(package_root, "R", "condition_grn_runtime_guard.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  contract_text <- paste(
    readLines(
      file.path(package_root, "R", "condition_grn_contract.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_match(runtime_text, 'package_version\\("1\\.6\\.3"\\)')
  expect_match(runtime_text, 'Config/Pando/NativeSparseABI', fixed = TRUE)
  expect_match(runtime_text, 'Config/Pando/NativeCallBinding', fixed = TRUE)
  expect_match(
    runtime_text,
    "registered-symbol-lookup-worker-safe-v1",
    fixed = TRUE
  )
  expect_match(
    runtime_text,
    "exact-refit-validation-sparse-residual-v1",
    fixed = TRUE
  )
  expect_match(
    runtime_text,
    "cpp-eigen-memory-bounded-hybrid-target-v1",
    fixed = TRUE
  )
  expect_match(
    runtime_text, "dense-direct-or-matrix-free-schur-pcg-v1", fixed = TRUE
  )
  expect_match(runtime_text, "no-full-p2-on-high-p-path-v1", fixed = TRUE)
  for (symbol in c(
    "_Pando_condition_product_matrix_cpp",
    "_Pando_condition_fit_multitask_path_cpp",
    "_Pando_condition_refit_path_cpp",
    "_Pando_condition_fit_target_engine_cpp",
    "_Pando_condition_native_self_test_cpp"
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
    contract_text,
    ".rc_require_pando_hybrid_runtime(BPPARAM = BPPARAM)",
    fixed = TRUE
  )
  expect_match(
    guard_text,
    ".rc_fit_condition_grns_by_cell_type_unchecked",
    fixed = TRUE
  )
})

test_that("an installed compatible Pando exposes the memory-bounded runtime", {
  skip_if_not_installed("Pando", minimum_version = "1.6.3")

  description <- utils::packageDescription("Pando")
  expect_identical(description[["Config/Pando/NativeSparseABI"]], "6")
  expect_identical(
    description[["Config/Pando/NativeCallBinding"]],
    "registered-symbol-lookup-worker-safe-v1"
  )
  expect_identical(
    description[["Config/Pando/ConditionRefitBackend"]],
    "dense-direct-or-matrix-free-schur-pcg-v1"
  )
  expect_identical(
    description[["Config/Pando/ConditionTargetEngineBackend"]],
    "cpp-eigen-memory-bounded-hybrid-target-v1"
  )
  expect_identical(
    description[["Config/Pando/ConditionInnerCVBackend"]],
    "exact-refit-validation-sparse-residual-v1"
  )
  expect_identical(
    description[["Config/Pando/ConditionMemoryContract"]],
    "no-full-p2-on-high-p-path-v1"
  )
  namespace <- asNamespace("Pando")
  for (wrapper in c(
    ".condition_product_matrix_cpp",
    ".condition_fit_multitask_path_cpp",
    ".condition_refit_path_cpp",
    ".condition_fit_target_engine_cpp",
    ".condition_native_self_test_cpp"
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
    "_Pando_condition_fit_target_engine_cpp",
    "_Pando_condition_native_self_test_cpp"
  )) {
    info <- getNativeSymbolInfo(
      symbol,
      PACKAGE = "Pando",
      withRegistrationInfo = TRUE
    )
    expect_true(is.list(info) && !is.null(info$address), info = symbol)
  }
  self_test <- Pando:::.condition_native_self_test_cpp()
  expect_true(self_test$passed)
  expect_true(self_test$budget_guard_passed)
  expect_true(self_test$numerical$hybrid_preconditioner)
  expect_true(self_test$numerical$budget_guard_passed)
  expect_lt(self_test$numerical$schur_refit_relative_error, 1e-8)
  expect_false(self_test$execution_plan$full_predictor_square_allocated)
  expect_silent(RegCompassR:::.rc_require_pando_hybrid_runtime())
})
