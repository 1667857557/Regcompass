test_that("RegCompass accepts only the canonical Pando fit schema", {
  canonical <- structure(
    list(schema_version = "pando_condition_grn_fit"),
    class = c("ConditionGRNFit", "list")
  )
  expect_silent(
    RegCompassR:::.rc_require_pando_condition_grn_fit(canonical)
  )
  for (schema in c(
    "pando_condition_grn_fit_v4",
    "pando_condition_grn_fit_v5"
  )) {
    legacy <- canonical
    legacy$schema_version <- schema
    expect_error(
      RegCompassR:::.rc_require_pando_condition_grn_fit(legacy),
      "canonical pando_condition_grn_fit"
    )
  }
})

test_that("canonical consumers are direct definitions", {
  consumers <- c(
    ".rc_extract_condition_grn_contract",
    ".rc_fit_condition_grns_by_cell_type",
    ".rc_validate_layer1_stage"
  )
  for (name in consumers) {
    fun <- get(name, envir = asNamespace("RegCompassR"), inherits = FALSE)
    text <- paste(deparse(body(fun)), collapse = "\n")
    expect_false(grepl("pando_condition_grn_fit_v4", text, fixed = TRUE))
    expect_false(grepl("pando_condition_grn_fit_v5", text, fixed = TRUE))
    expect_false(grepl("body(", text, fixed = TRUE))
    expect_false(grepl("assign(", text, fixed = TRUE))
  }
})
