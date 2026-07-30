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
      "version-suffixed schemas are unsupported"
    )
  }
})

test_that("runtime consumers contain no version-suffixed schema contract", {
  consumers <- c(
    ".rc_extract_condition_grn_contract_canonical",
    ".rc_validate_layer1_stage"
  )
  for (name in consumers) {
    fun <- get(name, envir = asNamespace("RegCompassR"), inherits = FALSE)
    text <- paste(deparse(body(fun)), collapse = "\n")
    expect_false(grepl("pando_condition_grn_fit_v4", text, fixed = TRUE))
    expect_false(grepl("pando_condition_grn_fit_v5", text, fixed = TRUE))
  }
})
