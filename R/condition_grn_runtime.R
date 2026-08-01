.RC_PANDO_MIN_VERSION <- package_version("1.6.1")
.RC_PANDO_NATIVE_SPARSE_ABI <- "5"
.RC_PANDO_TARGET_ENGINE_BACKEND <-
  "cpp-eigen-fused-hybrid-gram-nested-cv-path-refit-validation-stats-fail-fast"
.RC_PANDO_INNER_CV_BACKEND <-
  "cpp-eigen-hybrid-gram-sufficient-statistics-fail-fast"
.RC_PANDO_REFIT_BACKEND <- "cpp-eigen-direct-path-fail-fast"
.RC_PANDO_NATIVE_SYMBOLS <- c(
  "_Pando_condition_product_matrix_cpp",
  "_Pando_condition_fit_multitask_path_cpp",
  "_Pando_condition_refit_path_cpp",
  "_Pando_condition_fit_target_engine_cpp"
)

.rc_require_pando_hybrid_runtime <- function() {
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop(
      "Pando >= 1.6.1 with native condition ABI 5 is required. Install 1667857557/Pando_regcompass.",
      call. = FALSE
    )
  }
  installed <- utils::packageVersion("Pando")
  if (installed < .RC_PANDO_MIN_VERSION) {
    stop(
      sprintf(
        "Installed Pando %s is incompatible; RegCompass requires Pando >= 1.6.1 with native condition ABI 5.",
        as.character(installed)
      ),
      call. = FALSE
    )
  }
  description <- utils::packageDescription("Pando")
  required_fields <- c(
    "Config/Pando/NativeSparseABI",
    "Config/Pando/ConditionRefitBackend",
    "Config/Pando/ConditionTargetEngineBackend",
    "Config/Pando/ConditionInnerCVBackend"
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
      ". Install 1667857557/Pando_regcompass.",
      call. = FALSE
    )
  }
  expected <- c(
    "Config/Pando/NativeSparseABI" = .RC_PANDO_NATIVE_SPARSE_ABI,
    "Config/Pando/ConditionRefitBackend" = .RC_PANDO_REFIT_BACKEND,
    "Config/Pando/ConditionTargetEngineBackend" =
      .RC_PANDO_TARGET_ENGINE_BACKEND,
    "Config/Pando/ConditionInnerCVBackend" = .RC_PANDO_INNER_CV_BACKEND
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
      is.loaded,
      logical(1),
      PACKAGE = "Pando"
    )
  ]
  if (length(missing_symbols)) {
    stop(
      "Installed Pando is missing registered native condition kernel(s): ",
      paste(missing_symbols, collapse = ", "),
      ". Reinstall 1667857557/Pando_regcompass; no R fallback is permitted.",
      call. = FALSE
    )
  }
  if (!exists(
        "infer_condition_grn", envir = asNamespace("Pando"),
        inherits = FALSE
      )) {
    stop(
      "Installed Pando lacks infer_condition_grn(). Reinstall 1667857557/Pando_regcompass.",
      call. = FALSE
    )
  }
  invisible(list(
    version = as.character(installed),
    native_sparse_abi = observed[["Config/Pando/NativeSparseABI"]],
    target_engine_backend =
      observed[["Config/Pando/ConditionTargetEngineBackend"]],
    inner_cv_backend = observed[["Config/Pando/ConditionInnerCVBackend"]],
    refit_backend = observed[["Config/Pando/ConditionRefitBackend"]],
    native_symbols = .RC_PANDO_NATIVE_SYMBOLS
  ))
}
