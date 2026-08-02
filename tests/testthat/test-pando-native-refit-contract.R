test_that("RegCompass source requires the memory-bounded Pando ABI 6 contract", {
  package_root <- testthat::test_path("..", "..")
  description <- read.dcf(file.path(package_root, "DESCRIPTION"))
  expect_identical(unname(description[1L, "Version"]), "2.2.10")
  suggests <- unname(description[1L, "Suggests"])
  expect_match(suggests, "Pando", fixed = TRUE)
  expect_false(grepl("Pando (>=", suggests, fixed = TRUE))
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
  remotes <- gsub(
    "[[:space:]]+", " ",
    trimws(unname(description[1L, "Remotes"]))
  )
  expect_identical(
    remotes,
    paste(
      "1667857557/SuperCell_Seurat_V4,",
      "1667857557/Pando_regcompass"
    )
  )
  expect_false(grepl("Pando_regcompass@", remotes, fixed = TRUE))
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
  expect_false(grepl("RC_PANDO_MIN_VERSION", runtime_text, fixed = TRUE))
  expect_false(grepl('package_version("1.6.3")', runtime_text, fixed = TRUE))
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
  expect_match(
    runtime_text,
    "lightweight_metadata_symbol_api_v1",
    fixed = TRUE
  )
  for (symbol in c(
    "_Pando_condition_product_matrix_cpp",
    "_Pando_condition_fit_multitask_path_cpp",
    "_Pando_condition_refit_path_cpp",
    "_Pando_condition_fit_target_engine_cpp"
  )) {
    expect_match(runtime_text, symbol, fixed = TRUE)
  }
  expect_false(grepl("condition_native_self_test_cpp", runtime_text, fixed = TRUE))
  expect_false(grepl("schur_refit_relative_error", runtime_text, fixed = TRUE))
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

test_that("an installed compatible Pando passes lightweight runtime checks", {
  skip_if_not_installed("Pando")

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
  runtime <- RegCompassR:::.rc_require_pando_hybrid_runtime()
  expect_identical(
    runtime$runtime_check,
    "lightweight_metadata_symbol_api_v1"
  )
  expect_identical(length(runtime$native_symbols), 4L)
})
