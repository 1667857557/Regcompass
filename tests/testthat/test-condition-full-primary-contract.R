test_that("Pando fit validator separates coefficient and projection support", {
  topology <- matrix(
    TRUE, 2, 2,
    dimnames = list(c("common", "closed"), c("A", "B"))
  )
  estimable <- matrix(
    c(TRUE, TRUE, FALSE, FALSE),
    nrow = 2,
    byrow = TRUE,
    dimnames = dimnames(topology)
  )
  structural <- topology & !estimable
  fit <- structure(list(
    schema_version = "pando_condition_grn_fit",
    topology_mask = topology,
    estimability_mask = estimable,
    coefficient_estimable_mask = estimable,
    projectable_structural_zero_mask = structural,
    projection_support_mask = estimable | structural,
    primary_projection = "projection_condition_full_oof",
    nonestimable_projection_policy = "structural_zero_by_condition",
    projection_condition_full_oof = matrix(
      0, 4, 1,
      dimnames = list(paste0("c", 1:4), "target")
    ),
    projection_common_oof = matrix(
      0, 4, 1,
      dimnames = list(paste0("c", 1:4), "target")
    )
  ), class = c("ConditionGRNFit", "list"))
  expect_invisible(.rc_require_pando_condition_grn_fit(fit))
  expect_true(all(fit$projection_support_mask))
  expect_true(all(fit$projectable_structural_zero_mask["closed", ]))

  invalid <- fit
  invalid$projectable_structural_zero_mask["closed", "A"] <- FALSE
  expect_error(
    .rc_require_pando_condition_grn_fit(invalid),
    "contract is inconsistent"
  )
})

test_that("Layer 1 source promotes condition-full and retains common decomposition", {
  helper <- paste(
    deparse(body(.rc_condition_pando_projection)), collapse = "\n"
  )
  layer1 <- paste(
    deparse(body(.rc_cell_first_projection_layer1)), collapse = "\n"
  )
  expect_match(helper, "project_condition_grn_primary_cells", fixed = TRUE)
  expect_match(helper, "condition_unique = primary - common", fixed = TRUE)
  expect_match(layer1, "reaction_expression = reaction_primary", fixed = TRUE)
  expect_match(layer1, "reaction_expression_condition_full_oof", fixed = TRUE)
  expect_match(layer1, "reaction_expression_common_oof", fixed = TRUE)
  expect_match(layer1, "gene_projection_condition_unique_oof", fixed = TRUE)
})

test_that("Layer 2 source contains only canonical evidence routes", {
  text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  required <- c(
    "penalty_condition_full_oof",
    "penalty_common_oof",
    "penalty_condition_unique_increment",
    "penalty_rna_only",
    "structural_zero_by_condition"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  removed <- c(
    "penalty_depth_matched_rna",
    "penalty_common_depth_interval_rna",
    "penalty_alpha_sensitivity",
    "zero_support_sensitive",
    "link_saturation_sensitive"
  )
  expect_false(any(vapply(removed, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("condition-unique penalty increment is an exact decomposition", {
  full <- matrix(
    c(1, 2, 3, 4), 2, 2,
    dimnames = list(c("r1", "r2"), c("u1", "u2"))
  )
  common <- matrix(
    c(0.5, 2, 2.5, 3), 2, 2,
    dimnames = dimnames(full)
  )
  layer2 <- list(
    penalty = full,
    penalty_condition_full_oof = full,
    penalty_common_oof = common,
    penalty_condition_unique_increment = full - common,
    penalty_rna_only = full * 0 + 1,
    vmax = full * 0 + 1,
    feasible = matrix(TRUE, 2, 2, dimnames = dimnames(full)),
    evaluated = matrix(TRUE, 2, 2, dimnames = dimnames(full)),
    score = full,
    unit_meta = data.frame(pool_id = c("u1", "u2")),
    schema_version = "regcompass_regulatory_layer2_v2",
    comparison_contract = list(
      primary = "penalty_condition_full_oof",
      common_component = "penalty_common_oof",
      nonestimable_edge_policy = "structural_zero_by_condition",
      removed_guardrails = c(
        "depth_matching",
        "common_depth_restriction",
        "alpha_sensitivity",
        "zero_support_sensitivity",
        "link_saturation_propagation"
      )
    )
  )
  class(layer2) <- c("regcompass_layer2_step", "list")
  expect_equal(
    layer2$penalty_condition_unique_increment,
    layer2$penalty_condition_full_oof - layer2$penalty_common_oof
  )
})
