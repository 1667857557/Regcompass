# Consume the single canonical unversioned Pando condition-GRN fit schema.

.RC_PANDO_CONDITION_GRN_FIT_SCHEMA <- "pando_condition_grn_fit"

.rc_schema_transform <- local({
  old_schema <- c(
    "pando_condition_grn_fit_v4",
    "pando_condition_grn_fit_v5"
  )
  canonical_schema <- .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
  function(x) {
    if (is.character(x)) {
      x[x %in% old_schema] <- canonical_schema
      return(x)
    }
    if (is.call(x) || is.pairlist(x) || is.expression(x)) {
      for (i in seq_along(x)) {
        x[[i]] <- Recall(x[[i]])
      }
    }
    x
  }
})

.rc_schema_environment <- environment()
.rc_schema_functions <- setdiff(
  ls(envir = .rc_schema_environment, all.names = TRUE),
  c(".rc_schema_transform")
)
for (.rc_schema_name in .rc_schema_functions) {
  .rc_schema_value <- get(
    .rc_schema_name,
    envir = .rc_schema_environment,
    inherits = FALSE
  )
  if (!is.function(.rc_schema_value)) next
  body(.rc_schema_value) <- .rc_schema_transform(
    body(.rc_schema_value)
  )
  assign(
    .rc_schema_name,
    .rc_schema_value,
    envir = .rc_schema_environment
  )
}

.rc_require_pando_condition_grn_fit <- function(fit) {
  if (!inherits(fit, "ConditionGRNFit") ||
      !identical(
        fit$schema_version,
        .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
      )) {
    stop(
      "RegCompass requires the canonical pando_condition_grn_fit contract; ",
      "version-suffixed schemas are unsupported.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_extract_condition_grn_contract_canonical <-
  .rc_extract_condition_grn_contract
.rc_extract_condition_grn_contract <- function(...) {
  answer <- .rc_extract_condition_grn_contract_canonical(...)
  fits <- answer$fit_contracts %||% list()
  if (length(fits)) {
    invisible(lapply(fits, .rc_require_pando_condition_grn_fit))
  }
  answer$pando_fit_schema <- .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
  answer
}

rm(
  .rc_schema_environment,
  .rc_schema_functions,
  .rc_schema_name,
  .rc_schema_value,
  .rc_schema_transform
)
