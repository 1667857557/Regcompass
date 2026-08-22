test_that("recoverable Pando provenance is promoted to the fit contract", {
  dictionary <- data.frame(edge_id = "G||TF||chr1:1-2")
  attr(dictionary, "preprocessing_provenance_verified") <- TRUE
  fit <- structure(
    list(edge_dictionary = dictionary),
    class = c("ConditionGRNFit", "list")
  )

  completed <- .rc_complete_pando_condition_fit_contract(fit)

  expect_true(completed$dictionary_preprocessing_provenance_verified)
})

test_that("unverified dictionaries are not silently accepted", {
  dictionary <- data.frame(edge_id = "G||TF||chr1:1-2")
  fit <- structure(
    list(edge_dictionary = dictionary),
    class = c("ConditionGRNFit", "list")
  )

  completed <- .rc_complete_pando_condition_fit_contract(fit)

  expect_false(completed$dictionary_preprocessing_provenance_verified)
})

test_that("explicit provenance values remain authoritative", {
  dictionary <- data.frame(edge_id = "G||TF||chr1:1-2")
  attr(dictionary, "preprocessing_provenance_verified") <- TRUE
  fit <- structure(
    list(
      edge_dictionary = dictionary,
      dictionary_preprocessing_provenance_verified = FALSE
    ),
    class = c("ConditionGRNFit", "list")
  )

  completed <- .rc_complete_pando_condition_fit_contract(fit)

  expect_false(completed$dictionary_preprocessing_provenance_verified)
})


test_that("large condition contracts use indexed vectorized alignment", {
  validator <- paste(
    deparse(body(RegCompassR:::.rc_require_pando_condition_grn_fit)),
    collapse = "\n"
  )
  expect_match(validator, "coefficient_edge_index <- match", fixed = TRUE)
  expect_match(validator, "contrast_a_key <- paste", fixed = TRUE)
  expect_match(validator, "ia <- match", fixed = TRUE)
  expect_false(grepl(
    "for (i in seq_len(nrow(edge_inference)))", validator, fixed = TRUE
  ))
  expect_false(grepl(
    "for (i in seq_len(nrow(contrast)))", validator, fixed = TRUE
  ))
})
