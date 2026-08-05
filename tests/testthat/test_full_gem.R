.full_gem_test_solver <- function() {
  if (requireNamespace("highs", quietly = TRUE)) return("highs")
  if (requireNamespace("Rglpk", quietly = TRUE)) return("glpk")
  if (requireNamespace("gurobi", quietly = TRUE)) return("gurobi")
  NULL
}

.full_gem_test_model <- function() {
  S <- matrix(
    c(
      1, -1,  0,  0,  0,  0,
      0,  1, -1,  0,  0,  0,
      0,  0,  0,  1, -1,  0,
      0,  0,  0,  0,  1, -1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("A_c", "B_c", "C_c", "D_c"),
      c("EX_A", "R1", "EX_B", "EX_C", "R2", "EX_D")
    )
  )
  rc_make_gem(
    S,
    lb = stats::setNames(rep(0, 6), colnames(S)),
    ub = stats::setNames(rep(1000, 6), colnames(S)),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal", "exchange",
               "exchange", "internal", "exchange"),
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
}

.full_gem_test_medium <- function(open_c = FALSE) {
  data.frame(
    medium_scenario_id = "test",
    exchange_reaction_id = c("EX_A", "EX_B", "EX_C", "EX_D"),
    lb = 0,
    ub = c(1000, 1000, if (open_c) 1000 else 0, 1000),
    available = c(TRUE, TRUE, open_c, TRUE),
    condition = "all",
    stringsAsFactors = FALSE
  )
}

test_that("full GEM applies medium bounds without deleting reactions", {
  gem <- .full_gem_test_model()
  model <- rc_build_full_gem(
    gem,
    medium_table = .full_gem_test_medium()
  )

  expect_identical(colnames(model$S), colnames(gem$S))
  expect_identical(rownames(model$S), rownames(gem$S))
  expect_equal(model$S, gem$S)
  expect_equal(unname(model$ub[["EX_C"]]), 0)
  expect_identical(
    model$build_params$strategy,
    "compass_medium_constrained_full_gem"
  )
  expect_identical(model$build_params$n_input_reactions, 6L)
  expect_identical(model$build_params$n_reactions, 6L)
  expect_identical(model$build_params$n_medium_removed_reactions, 0L)
  expect_false(model$build_params$medium_direct_reaction_deletion)
  expect_false(model$build_params$flux_consistency_pruning)
  expect_identical(
    model$build_params$medium_handling,
    "exchange_bounds_only_no_reaction_deletion"
  )
  expect_identical(model$target_status, "not_prechecked")
})

test_that("full GEM cache retains medium-blocked targets", {
  gem <- .full_gem_test_model()
  dirs <- data.frame(
    reaction_id = c("R1", "R2"),
    target_direction = c("forward", "forward"),
    stringsAsFactors = FALSE
  )
  cache <- rc_build_full_gem_cache(
    gem = gem,
    dirs = dirs,
    medium_scenarios = .full_gem_test_medium(),
    solver = "highs",
    time_limit = 60,
    flux_consistency_epsilon = 1e-8
  )
  summary <- attr(cache, "summary")

  expect_equal(nrow(summary), 1L)
  expect_identical(
    summary$build_strategy,
    "compass_medium_constrained_full_gem"
  )
  expect_equal(summary$n_input_reactions, 6L)
  expect_equal(summary$n_reactions, 6L)
  expect_equal(summary$n_medium_removed_reactions, 0L)
  expect_true(summary$medium_applied)
  expect_setequal(
    vapply(cache, `[[`, character(1), "reaction_id"),
    c("R1", "R2")
  )
  expect_true(file.exists(summary$file[[1L]]))
  expect_identical(attr(cache, "completion_method"), "none")
  expect_false(attr(cache, "fastcore_executed"))
  expect_false(attr(cache, "corda2_executed"))
  expect_identical(
    attr(cache, "medium_handling"),
    "exchange_bounds_only_no_reaction_deletion"
  )
})

test_that("medium-blocked full-GEM target is handled by COMPASS vmax", {
  solver <- .full_gem_test_solver()
  skip_if(is.null(solver), "No LP solver is installed")
  model <- rc_build_full_gem(
    .full_gem_test_model(),
    medium_table = .full_gem_test_medium()
  )

  feasible <- rc_compass_vmax_directional(
    model$S, model$lb, model$ub,
    target_reaction = "R1",
    direction = "forward",
    solver = solver,
    flux_threshold = 1e-8
  )
  blocked <- rc_compass_vmax_directional(
    model$S, model$lb, model$ub,
    target_reaction = "R2",
    direction = "forward",
    solver = solver,
    flux_threshold = 1e-8
  )
  step2 <- .rc_compass_step2_from_vmax_directional(
    S = model$S,
    lb = model$lb,
    ub = model$ub,
    target_reaction = "R2",
    penalties = stats::setNames(rep(0, ncol(model$S)), colnames(model$S)),
    vmax_result = blocked,
    target_direction = "forward",
    solver = solver,
    flux_threshold = 1e-8
  )

  expect_true(feasible$feasible)
  expect_gt(feasible$vmax, 0)
  expect_false(blocked$feasible)
  expect_lt(blocked$vmax, 1e-8)
  expect_false(step2$feasible)
  expect_identical(step2$step2_status, "not_run")
})

test_that("all three routes reject direct medium reaction deletion", {
  solver <- .full_gem_test_solver()
  skip_if(is.null(solver), "No LP solver is installed")
  gem <- .full_gem_test_model()
  medium <- .full_gem_test_medium()
  reactions <- colnames(gem$S)

  full <- rc_build_full_gem(gem, medium_table = medium)
  fastcore_parent <- .rc_fastcore_parent(
    gem,
    medium_table = medium,
    forbidden_roles = character(),
    solver = solver,
    time_limit = 60,
    fastcore_epsilon = 1e-4
  )
  corda_parent <- .rc_corda_parent(
    gem,
    medium_table = medium,
    forbidden_roles = character(),
    solver = solver,
    time_limit = 60
  )

  expect_identical(colnames(full$S), reactions)
  expect_identical(colnames(fastcore_parent$S), reactions)
  expect_identical(colnames(corda_parent$S), reactions)
  expect_identical(corda_parent$corda_parent_prepruning, "none")
  expect_identical(corda_parent$corda_parent_role_blocking, "none")
})

test_that("full GEM completion contract is never labelled FASTCORE", {
  answer <- list(
    model_mode = "full_gem",
    params = list(),
    model_cache_summary = data.frame(
      medium_scenario = "base",
      condition = "all",
      n_input_reactions = 6L,
      n_reactions = 6L,
      n_medium_removed_reactions = 0L,
      n_medium_bound_changes = 1L,
      medium_applied = TRUE,
      medium_handling = "exchange_bounds_only_no_reaction_deletion",
      build_strategy = "compass_medium_constrained_full_gem",
      stringsAsFactors = FALSE
    )
  )

  result <- .rc_layer2_finalize_completion(
    answer,
    corda_options = list(),
    is_corda2 = FALSE,
    solver = "highs"
  )

  expect_identical(result$params$model_completion, "none")
  expect_identical(result$params$structural_completion, "none")
  expect_identical(
    result$params$structural_completion_algorithm,
    "compass_medium_bounds_only"
  )
  expect_false(result$params$medium_direct_reaction_deletion)
  expect_identical(result$completion_contract$model_completion, "none")
  expect_false(result$completion_contract$flux_consistency_pruning)
  expect_false(result$completion_contract$fastcore_executed)
  expect_false(result$completion_contract$corda2_executed)
})

test_that("full GEM rejects FASTCORE and CORDA2 completion controls", {
  expect_error(
    .rc_layer2_prepare_completion(
      layer1 = list(),
      meta_modules = list(),
      model_mode = "full_gem",
      layer2_args = list(
        model_params = list(model_completion = "fastcore")
      )
    ),
    "automatically skips FASTCORE"
  )
  expect_error(
    .rc_layer2_prepare_completion(
      layer1 = list(),
      meta_modules = list(),
      model_mode = "full_gem",
      layer2_args = list(
        model_params = list(model_completion = "corda2")
      )
    ),
    "automatically skips FASTCORE"
  )
  expect_error(
    .rc_layer2_prepare_completion(
      layer1 = list(),
      meta_modules = list(),
      model_mode = "full_gem",
      layer2_args = list(
        model_params = list(fastcore_epsilon = 1e-4)
      )
    ),
    "does not accept FASTCORE or CORDA2 controls"
  )
})
