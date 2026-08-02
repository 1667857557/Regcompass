.RC_PANDO_MIN_VERSION <- package_version("1.6.3")
.RC_PANDO_NATIVE_SPARSE_ABI <- "6"
.RC_PANDO_NATIVE_CALL_BINDING <-
  "registered-symbol-lookup-worker-safe-v1"
.RC_PANDO_TARGET_ENGINE_BACKEND <-
  "cpp-eigen-memory-bounded-hybrid-target-v1"
.RC_PANDO_INNER_CV_BACKEND <-
  "exact-refit-validation-sparse-residual-v1"
.RC_PANDO_REFIT_BACKEND <-
  "dense-direct-or-matrix-free-schur-pcg-v1"
.RC_PANDO_MEMORY_CONTRACT <- "no-full-p2-on-high-p-path-v1"
.RC_PANDO_NATIVE_SYMBOLS <- c(
  "_Pando_condition_product_matrix_cpp",
  "_Pando_condition_fit_multitask_path_cpp",
  "_Pando_condition_refit_path_cpp",
  "_Pando_condition_fit_target_engine_cpp",
  "_Pando_condition_native_self_test_cpp"
)
.RC_PANDO_NATIVE_WRAPPERS <- c(
  ".condition_product_matrix_cpp",
  ".condition_fit_multitask_path_cpp",
  ".condition_refit_path_cpp",
  ".condition_fit_target_engine_cpp",
  ".condition_native_self_test_cpp"
)

.rc_pando_registered_symbol_available <- function(symbol) {
  tryCatch({
    info <- getNativeSymbolInfo(
      symbol,
      PACKAGE = "Pando",
      withRegistrationInfo = TRUE
    )
    is.list(info) && !is.null(info$address)
  }, error = function(error) FALSE)
}

.rc_pando_worker_runtime_probe <- function(
    symbols, wrappers, expected, library_paths) {
  .libPaths(library_paths)
  loadNamespace("Pando")
  description <- utils::packageDescription("Pando")
  observed <- vapply(names(expected), function(name) {
    value <- description[[name]]
    if (is.null(value)) "" else as.character(value[[1L]])
  }, character(1))
  incompatible <- names(expected)[observed != expected]
  if (length(incompatible)) {
    stop(
      "Worker loaded incompatible Pando runtime metadata: ",
      paste0(
        incompatible, "=", observed[incompatible],
        " (required ", expected[incompatible], ")",
        collapse = "; "
      ),
      call. = FALSE
    )
  }
  namespace <- asNamespace("Pando")
  missing_wrappers <- wrappers[!vapply(
    wrappers,
    exists,
    logical(1),
    envir = namespace,
    inherits = FALSE
  )]
  if (length(missing_wrappers)) {
    stop(
      "Worker Pando namespace is missing native wrapper(s): ",
      paste(missing_wrappers, collapse = ", "),
      call. = FALSE
    )
  }
  unsafe_wrappers <- wrappers[!vapply(wrappers, function(name) {
    value <- get(name, envir = namespace, inherits = FALSE)
    grepl(
      ".pando_registered_call",
      paste(deparse(body(value)), collapse = "\n"),
      fixed = TRUE
    )
  }, logical(1))]
  if (length(unsafe_wrappers)) {
    stop(
      "Worker Pando uses namespace-bound Rcpp wrapper(s): ",
      paste(unsafe_wrappers, collapse = ", "),
      call. = FALSE
    )
  }
  missing_symbols <- symbols[!vapply(symbols, function(symbol) {
    tryCatch({
      info <- getNativeSymbolInfo(
        symbol,
        PACKAGE = "Pando",
        withRegistrationInfo = TRUE
      )
      is.list(info) && !is.null(info$address)
    }, error = function(error) FALSE)
  }, logical(1))]
  if (length(missing_symbols)) {
    stop(
      "Worker Pando DLL is missing registered native kernel(s): ",
      paste(missing_symbols, collapse = ", "),
      call. = FALSE
    )
  }
  self_test <- get(
    ".condition_native_self_test_cpp",
    envir = namespace,
    inherits = FALSE
  )()
  if (!is.list(self_test) || !isTRUE(self_test$passed) ||
      !isTRUE(self_test$budget_guard_passed) ||
      !isTRUE(self_test$numerical$passed) ||
      !is.numeric(self_test$numerical$schur_refit_relative_error) ||
      length(self_test$numerical$schur_refit_relative_error) != 1L ||
      !is.finite(self_test$numerical$schur_refit_relative_error) ||
      self_test$numerical$schur_refit_relative_error >= 1e-8 ||
      isTRUE(self_test$execution_plan$full_predictor_square_allocated) ||
      !identical(
        self_test$execution_plan$memory_contract,
        "no_full_p2_on_high_p_path_v1"
      )) {
    stop("Worker Pando native numerical self-test failed.", call. = FALSE)
  }
  TRUE
}

.rc_validate_pando_worker_runtime <- function(BPPARAM) {
  if (is.null(BPPARAM) || identical(BPPARAM, FALSE)) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop(
      "BiocParallel is required to validate the supplied `BPPARAM`.",
      call. = FALSE
    )
  }
  workers <- suppressWarnings(as.integer(BiocParallel::bpnworkers(BPPARAM)))
  if (!length(workers) || is.na(workers) || workers < 1L) workers <- 1L
  probes <- min(workers, 2L)
  result <- tryCatch(
    BiocParallel::bplapply(
      seq_len(probes),
      function(index, symbols, wrappers, expected, library_paths) {
        .rc_pando_worker_runtime_probe(
          symbols = symbols,
          wrappers = wrappers,
          expected = expected,
          library_paths = library_paths
        )
      },
      symbols = .RC_PANDO_NATIVE_SYMBOLS,
      wrappers = .RC_PANDO_NATIVE_WRAPPERS,
      expected = c(
        "Config/Pando/NativeSparseABI" = .RC_PANDO_NATIVE_SPARSE_ABI,
        "Config/Pando/NativeCallBinding" = .RC_PANDO_NATIVE_CALL_BINDING,
        "Config/Pando/ConditionRefitBackend" = .RC_PANDO_REFIT_BACKEND,
        "Config/Pando/ConditionTargetEngineBackend" =
          .RC_PANDO_TARGET_ENGINE_BACKEND,
        "Config/Pando/ConditionInnerCVBackend" =
          .RC_PANDO_INNER_CV_BACKEND,
        "Config/Pando/ConditionMemoryContract" =
          .RC_PANDO_MEMORY_CONTRACT
      ),
      library_paths = .libPaths(),
      BPPARAM = BPPARAM
    ),
    error = function(error) {
      stop(
        "Pando native runtime failed on a BiocParallel worker before GRN ",
        "fitting: ", conditionMessage(error),
        ". Reinstall 1667857557/Pando_regcompass >= 1.6.3 and restart R.",
        call. = FALSE
      )
    }
  )
  if (!all(vapply(result, isTRUE, logical(1)))) {
    stop(
      "Pando native runtime worker probe returned an invalid result.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_require_pando_hybrid_runtime <- function(BPPARAM = NULL) {
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop(
      "Pando >= 1.6.3 with native condition ABI 6 and the high-p memory contract is required. Install 1667857557/Pando_regcompass.",
      call. = FALSE
    )
  }
  installed <- utils::packageVersion("Pando")
  if (installed < .RC_PANDO_MIN_VERSION) {
    stop(
      sprintf(
        "Installed Pando %s is incompatible; RegCompass requires Pando >= 1.6.3 with worker-safe memory-bounded native calls.",
        as.character(installed)
      ),
      call. = FALSE
    )
  }
  description <- utils::packageDescription("Pando")
  required_fields <- c(
    "Config/Pando/NativeSparseABI",
    "Config/Pando/NativeCallBinding",
    "Config/Pando/ConditionRefitBackend",
    "Config/Pando/ConditionTargetEngineBackend",
    "Config/Pando/ConditionInnerCVBackend",
    "Config/Pando/ConditionMemoryContract"
  )
  missing_fields <- required_fields[
    !required_fields %in% names(description) |
      !nzchar(vapply(required_fields, function(name) {
        value <- description[[name]]
        if (is.null(value)) "" else as.character(value[[1L]])
      }, character(1)))
  ]
  if (length(missing_fields)) {
    stop(
      "Installed Pando is missing required native runtime metadata: ",
      paste(missing_fields, collapse = ", "),
      ". Install 1667857557/Pando_regcompass >= 1.6.3.",
      call. = FALSE
    )
  }
  expected <- c(
    "Config/Pando/NativeSparseABI" = .RC_PANDO_NATIVE_SPARSE_ABI,
    "Config/Pando/NativeCallBinding" = .RC_PANDO_NATIVE_CALL_BINDING,
    "Config/Pando/ConditionRefitBackend" = .RC_PANDO_REFIT_BACKEND,
    "Config/Pando/ConditionTargetEngineBackend" =
      .RC_PANDO_TARGET_ENGINE_BACKEND,
    "Config/Pando/ConditionInnerCVBackend" = .RC_PANDO_INNER_CV_BACKEND,
    "Config/Pando/ConditionMemoryContract" = .RC_PANDO_MEMORY_CONTRACT
  )
  observed <- vapply(names(expected), function(name) {
    as.character(description[[name]][[1L]])
  }, character(1))
  incompatible <- names(expected)[observed != expected]
  if (length(incompatible)) {
    detail <- paste0(
      incompatible, "=", observed[incompatible],
      " (required ", expected[incompatible], ")"
    )
    stop(
      "Installed Pando native runtime is incompatible: ",
      paste(detail, collapse = "; "),
      ". No R fallback is permitted.",
      call. = FALSE
    )
  }
  missing_symbols <- .RC_PANDO_NATIVE_SYMBOLS[
    !vapply(
      .RC_PANDO_NATIVE_SYMBOLS,
      .rc_pando_registered_symbol_available,
      logical(1)
    )
  ]
  if (length(missing_symbols)) {
    stop(
      "Installed Pando is missing registered native condition kernel(s): ",
      paste(missing_symbols, collapse = ", "),
      ". Reinstall 1667857557/Pando_regcompass >= 1.6.3 and restart R; no R fallback is permitted.",
      call. = FALSE
    )
  }
  namespace <- asNamespace("Pando")
  missing_wrappers <- .RC_PANDO_NATIVE_WRAPPERS[!vapply(
    .RC_PANDO_NATIVE_WRAPPERS,
    exists,
    logical(1),
    envir = namespace,
    inherits = FALSE
  )]
  if (length(missing_wrappers)) {
    stop(
      "Installed Pando is missing native wrapper(s): ",
      paste(missing_wrappers, collapse = ", "),
      call. = FALSE
    )
  }
  unsafe_wrappers <- .RC_PANDO_NATIVE_WRAPPERS[!vapply(
    .RC_PANDO_NATIVE_WRAPPERS,
    function(name) {
      value <- get(name, envir = namespace, inherits = FALSE)
      grepl(
        ".pando_registered_call",
        paste(deparse(body(value)), collapse = "\n"),
        fixed = TRUE
      )
    },
    logical(1)
  )]
  if (length(unsafe_wrappers)) {
    stop(
      "Installed Pando still uses namespace-bound Rcpp native wrappers: ",
      paste(unsafe_wrappers, collapse = ", "),
      ". Install Pando >= 1.6.3 and restart R.",
      call. = FALSE
    )
  }
  if (!exists(
        "infer_condition_grn", envir = namespace,
        inherits = FALSE
      )) {
    stop(
      "Installed Pando lacks infer_condition_grn(). Reinstall 1667857557/Pando_regcompass.",
      call. = FALSE
    )
  }
  self_test <- get(
    ".condition_native_self_test_cpp",
    envir = namespace,
    inherits = FALSE
  )()
  if (!is.list(self_test) || !isTRUE(self_test$passed) ||
      !isTRUE(self_test$budget_guard_passed) ||
      !isTRUE(self_test$numerical$passed) ||
      !is.numeric(self_test$numerical$schur_refit_relative_error) ||
      length(self_test$numerical$schur_refit_relative_error) != 1L ||
      !is.finite(self_test$numerical$schur_refit_relative_error) ||
      self_test$numerical$schur_refit_relative_error >= 1e-8 ||
      isTRUE(self_test$execution_plan$full_predictor_square_allocated) ||
      !identical(
        self_test$execution_plan$memory_contract,
        "no_full_p2_on_high_p_path_v1"
      )) {
    stop(
      "Installed Pando failed the native high-p numerical/allocation self-test.",
      call. = FALSE
    )
  }
  .rc_validate_pando_worker_runtime(BPPARAM)
  invisible(list(
    version = as.character(installed),
    native_sparse_abi = observed[["Config/Pando/NativeSparseABI"]],
    native_call_binding = observed[["Config/Pando/NativeCallBinding"]],
    target_engine_backend =
      observed[["Config/Pando/ConditionTargetEngineBackend"]],
    inner_cv_backend = observed[["Config/Pando/ConditionInnerCVBackend"]],
    refit_backend = observed[["Config/Pando/ConditionRefitBackend"]],
    memory_contract = observed[["Config/Pando/ConditionMemoryContract"]],
    native_self_test = self_test,
    native_symbols = .RC_PANDO_NATIVE_SYMBOLS,
    native_wrappers = .RC_PANDO_NATIVE_WRAPPERS
  ))
}
