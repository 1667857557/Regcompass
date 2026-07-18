test_that("Human-GEM MAR exchange reactions can be used by medium scenarios", {
  S <- Matrix::Matrix(c(-1, 0, 1), nrow = 3, ncol = 1, sparse = TRUE)
  rownames(S) <- c("MAM00001e", "MAM00001c", "MAM00002c")
  colnames(S) <- "MAR09034"
  gem <- list(
    S = S,
    lb = c(MAR09034 = -1000),
    ub = c(MAR09034 = 1000),
    reaction_meta = data.frame(
      reaction_id = "MAR09034",
      reaction_name = "Exchange of glucose",
      stringsAsFactors = FALSE
    )
  )
  gem <- rc_annotate_reaction_roles(gem)
  expect_equal(gem$reaction_meta$role, "exchange")
  medium <- rc_make_medium_scenarios(
    gem, scenario = "permissive_all_exchange"
  )
  expect_gt(nrow(medium), 0)
  expect_true(all(medium$assumption_level == "technical_sensitivity_baseline"))
  named_medium <- rc_make_medium_scenarios(
    gem,
    scenario = "normal_human_plasma",
    strict_preset_matching = FALSE
  )
  expect_true(all(named_medium$medium_scenario_id == "normal_human_plasma"))
  expect_true(all(
    named_medium$evidence_source == "literature_backed_medium_catalog"
  ))
})

test_that("medium scenarios fail informatively when no exchange reactions are annotated", {
  S <- Matrix::Matrix(c(-1), nrow = 1, ncol = 1, sparse = TRUE)
  rownames(S) <- "m1c"
  colnames(S) <- "MAR00001"
  gem <- list(S = S, lb = c(MAR00001 = 0), ub = c(MAR00001 = 1000))
  gem <- rc_annotate_reaction_roles(gem, infer_from_id = FALSE)
  expect_error(rc_make_medium_scenarios(gem), "No.*exchange reactions")
})

test_that("pure custom medium does not require a generic exchange baseline", {
  S <- Matrix::Matrix(c(-1), nrow = 1, ncol = 1, sparse = TRUE)
  rownames(S) <- "m1c"
  colnames(S) <- "R1"
  gem <- list(S = S, lb = c(R1 = 0), ub = c(R1 = 1000))
  custom <- data.frame(
    medium_scenario_id = "measured",
    exchange_reaction_id = "R1",
    lb = 0,
    ub = 10,
    available = TRUE,
    stringsAsFactors = FALSE
  )
  output <- rc_make_medium_scenarios(
    gem, scenario = "custom", custom_medium = custom
  )
  expect_equal(output$medium_scenario_id, "measured")
  expect_equal(output$evidence_source, "user_supplied_custom_medium")
})
