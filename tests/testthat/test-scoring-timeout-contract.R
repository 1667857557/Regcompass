test_that("scoring APIs expose no time-limit parameter", {
  expect_false("time_limit" %in% names(formals(rc_compass_two_step_lp_directional)))
  expect_false("time_limit" %in% names(formals(rc_run_microcompass)))
  expect_false("time_limit" %in% names(formals(.rc_run_microcompass_engine)))
  expect_false("time_limit" %in% names(formals(.rc_score_existing_union_cache)))
})

test_that("only union-GEM construction receives the completion limit", {
  engine <- paste(deparse(body(.rc_run_microcompass_engine)), collapse = "\n")
  scoring <- paste(
    deparse(body(rc_compass_two_step_lp_directional)), collapse = "\n"
  )
  target_scoring <- paste(
    deparse(body(.rc_score_existing_union_cache)), collapse = "\n"
  )

  expect_match(engine, "model_params$completion_time_limit", fixed = TRUE)
  expect_match(engine, ".rc_build_medium_specific_union_gem_cache", fixed = TRUE)
  expect_false(grepl("time_limit =", scoring, fixed = TRUE))
  expect_false(grepl("time_limit =", target_scoring, fixed = TRUE))
})

test_that("Stage 5 rejects retired timeout arguments", {
  body_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  expect_match(body_text, "Scoring `time_limit` has been removed", fixed = TRUE)
  expect_match(body_text, "completion_time_limit", fixed = TRUE)
  expect_match(body_text, "allowed_model_params", fixed = TRUE)
  expect_false("time_limit" %in% c(
    "model_params", "omega", "target_direction", "solver", "flux_threshold"
  ))
})

test_that("target-union results record unlimited scoring", {
  target_body <- paste(
    deparse(body(.rc_score_existing_union_cache)), collapse = "\n"
  )
  expect_match(target_body, 'scoring_time_limit = "none"', fixed = TRUE)
})
