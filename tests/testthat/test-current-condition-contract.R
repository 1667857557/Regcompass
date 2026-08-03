test_that("Layer 1 exposes one current projection route", {
  source <- paste(readLines("../../R/layer1_regulatory_support.R"),
                  collapse = "\n")
  expect_match(source, "gene_projection = projection\\$projection")
  expect_match(source, "reaction_expression = reaction_multiome")
  expect_false(grepl("condition_full_oof|common_oof|condition_unique", source))
})

test_that("Layer 2 exposes primary and RNA-only routes", {
  source <- paste(readLines("../../R/step_layer2.R"), collapse = "\n")
  expect_match(source, 'primary = "penalty"', fixed = TRUE)
  expect_match(source, 'rna_control = "penalty_rna_only"', fixed = TRUE)
  expect_false(grepl("penalty_common_oof|penalty_condition_unique", source))
})

test_that("current Pando condition API is used", {
  source <- paste(
    readLines("../../R/condition_grn_contract.R"), collapse = "\n"
  )
  expect_match(source, "Pando::condition_grn_fit(grn_object)", fixed = TRUE)
  expect_match(source, "Pando::project_condition_grn_cells", fixed = TRUE)
  expect_false(grepl("network_name = \\"regcompass_condition_grn\\"", source))
})

test_that("retired user arguments are absent", {
  functions <- list(
    rc_regcompass_step_layer1,
    rc_regcompass_step_grn,
    rc_extract_pando_tf_peak_gene
  )
  arguments <- unique(unlist(lapply(functions, function(fun) names(formals(fun)))))
  expect_false(any(c(
    "projection_component", "comparison_support", "regulatory_alpha",
    "min_abs_estimate", "min_model_rsq"
  ) %in% arguments))
})
