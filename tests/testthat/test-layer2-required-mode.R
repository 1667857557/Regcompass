test_that("Layer 2 validation enforces a requested model mode", {
  value <- matrix(
    1,
    nrow = 1L,
    ncol = 1L,
    dimnames = list("R1", "U1")
  )
  zero <- value * 0
  flag <- matrix(
    TRUE,
    nrow = 1L,
    ncol = 1L,
    dimnames = dimnames(value)
  )
  layer2 <- structure(
    list(
      penalty = value,
      vmax = value,
      feasible = flag,
      evaluated = flag,
      score = value,
      unit_meta = data.frame(pool_id = "U1", stringsAsFactors = FALSE),
      model_mode = "meta_module_gem",
      penalty_condition_full_oof = value,
      penalty_common_oof = zero,
      penalty_rna_only = value,
      score_condition_full_oof = value,
      score_common_oof = zero,
      score_rna_only = value,
      penalty_condition_unique_increment = value,
      score_condition_unique_increment = value,
      comparison_contract = list(
        primary_route = "condition_full_oof",
        shared_structure = TRUE,
        shared_medium = TRUE,
        shared_directional_vmax = TRUE,
        normalization = "test",
        score_transform = "test",
        removed_guardrails = character()
      )
    ),
    class = "regcompass_layer2_step"
  )

  expect_identical(
    RegCompassR:::.rc_validate_layer2_stage(
      layer2,
      required_mode = "meta_module_gem"
    ),
    "U1"
  )
  expect_error(
    RegCompassR:::.rc_validate_layer2_stage(
      layer2,
      required_mode = "full_gem"
    ),
    "`full_gem` is required",
    fixed = TRUE
  )
  expect_error(
    RegCompassR:::.rc_validate_layer2_stage(
      layer2,
      required_mode = "unsupported"
    ),
    "`required_mode` must be one of",
    fixed = TRUE
  )
})
