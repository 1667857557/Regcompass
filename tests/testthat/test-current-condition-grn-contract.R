test_that("condition bridge exposes only common-dictionary controls", {
  defaults <- eval(
    formals(.rc_fit_condition_grns_by_cell_type)$pando_infer_args
  )
  expect_identical(defaults$tf_cor, 0.1)
  expect_identical(defaults$peak_cor, 0)
  expect_identical(defaults$adjust_method, "BH")
  expect_identical(defaults$padj_threshold, 0.05)
  expect_identical(defaults$rank_action, "mark")
  expect_identical(defaults$min_residual_df, 1L)
  expect_false(any(c(
    "candidate_screen", "condition_mix", "condition_weight", "alpha",
    "nlambda", "outer_nfolds", "inner_nfolds", "scale",
    "engine_control"
  ) %in% names(defaults)))

  implementation <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  expect_match(implementation, "Pando::infer_condition_grn", fixed = TRUE)
  expect_match(implementation, "adjust_method = \"BH\"", fixed = TRUE)
  expect_match(implementation, "padj_threshold = 0.05", fixed = TRUE)
  expect_match(implementation, ".rc_extract_condition_grn_contract",
               fixed = TRUE)
})

test_that("contract extraction preserves all effects and filters penalty effects", {
  implementation <- paste(
    deparse(body(.rc_extract_condition_grn_contract)), collapse = "\n"
  )
  expect_match(implementation, "condition_estimate", fixed = TRUE)
  expect_match(implementation, "penalty_effect", fixed = TRUE)
  expect_match(implementation, "penalty_eligible", fixed = TRUE)
  expect_match(implementation, "min_model_rsq", fixed = TRUE)
  expect_match(
    implementation,
    "same_exact_edge_dictionary_unscaled_gaussian_glm",
    fixed = TRUE
  )
  expect_false(grepl("estimability_mask", implementation, fixed = TRUE))
  expect_false(grepl("projectable_structural_zero", implementation,
                     fixed = TRUE))
  expect_false(grepl("target_rsq_oof_pooled", implementation, fixed = TRUE))
})

test_that("condition Layer 1 is cell-first and BH-filtered", {
  implementation <- paste(
    deparse(body(.rc_condition_pando_projection)), collapse = "\n"
  )
  expect_match(implementation, "project_condition_grn_cells", fixed = TRUE)
  expect_match(implementation, "significant_only = TRUE", fixed = TRUE)
  expect_match(implementation, "aggregate_condition_grn_projection",
               fixed = TRUE)
  expect_match(implementation, "sqrt(pmin(1, pmax(0", fixed = TRUE)
  expect_false(grepl("support_policy", implementation, fixed = TRUE))
})

test_that("metacell builder retains the grouped WNN API", {
  builder <- paste(
    deparse(body(.rc_build_grouped_wnn_membership)), collapse = "\n"
  )
  expect_match(builder, "SCimplify_by_graph_group", fixed = TRUE)
  expect_match(builder, "cell.graph.group", fixed = TRUE)
  expect_match(builder, "cell.split.condition", fixed = TRUE)
  expect_match(builder, "return.group.results = FALSE", fixed = TRUE)
})

test_that("regulatory integration is bounded and zero preserving", {
  rna <- matrix(
    c(0, 0.5, 1), nrow = 1,
    dimnames = list("g1", c("u1", "u2", "u3"))
  )
  modifier <- matrix(
    c(1, 1, -1), nrow = 1, dimnames = dimnames(rna)
  )
  observed <- .rc_integrate_regulatory_support(rna, modifier, alpha = 1)
  expect_equal(observed[1, 1], 0)
  expect_equal(observed[1, 3], 1)
  expect_true(observed[1, 2] > 0.5)
  expect_true(all(observed >= 0 & observed <= 1))
})
