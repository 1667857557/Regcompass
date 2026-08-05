.full_gem_test_solver <- function() {
  if (requireNamespace("highs", quietly = TRUE)) return("highs")
  if (requireNamespace("Rglpk", quietly = TRUE)) return("glpk")
  if (requireNamespace("gurobi", quietly = TRUE)) return("gurobi")
  NULL
}

.full_gem_test_model <- function() {
  S <- matrix(
    c(
      1, -1,  0,  0,
      0,  1, -1,  0,
      0,  0,  0, -1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("A_c", "B_c", "C_c"),
      c("EX_A", "R1", "EX_B", "R_BLOCKED")
    )
  )
  rc_make_gem(
    S,
    lb = stats::setNames(rep(0, 4), colnames(S)),
    ub = stats::setNames(rep(1000, 4), colnames(S)),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal", "exchange", "internal"),
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
}

test_that("full GEM removes only medium flux-inconsistent reactions", {
  solver <- .full_gem_test_solver()
  skip_if(is.null(solver), "No LP solver is installed")
  gem <- .full_gem_test_model()

  model <- rc_build_full_gem(
    gem,
    prune_flux_inconsistent = TRUE,
    solver = solver,
    time_limit = 60,
    flux_consistency_epsilon = 1e-8
  )

  expect_setequal(colnames(model$S), c("EX_A", "R1", "EX_B"))
  expect_false("R_BLOCKED" %in% colnames(model$S))
  expect_identical(
    model$build_params$strategy,
    "medium_flux_consistency_pruned_full_gem"
  )
  expect_false(model$build_params$fastcore_executed)
  expect_false(model$build_params$corda2_executed)
  expect_equal(model$build_params$n_flux_inconsistent_reactions, 1L)
})

test_that("full GEM cache excludes pruned targets and records the builder", {
  solver <- .full_gem_test_solver()
  skip_if(is.null(solver), "No LP solver is installed")
  gem <- .full_gem_test_model()
  dirs <- data.frame(
    reaction_id = c("R1", "R_BLOCKED"),
    target_direction = c("forward", "forward"),
    stringsAsFactors = FALSE
  )
  medium <- data.frame(
    medium_scenario_id = "base",
    exchange_reaction_id = NA_character_,
    lb = NA_real_, ub = NA_real_, available = FALSE,
    .no_constraints = TRUE,
    stringsAsFactors = FALSE
  )

  cache <- rc_build_full_gem_cache(
    gem = gem,
    dirs = dirs,
    medium_scenarios = medium,
    solver = solver,
    time_limit = 60,
    flux_consistency_epsilon = 1e-8
  )
  summary <- attr(cache, "summary")

  expect_equal(nrow(summary), 1L)
  expect_identical(
    summary$build_strategy,
    "medium_flux_consistency_pruned_full_gem"
  )
  expect_equal(summary$n_input_reactions, 4L)
  expect_equal(summary$n_reactions, 3L)
  expect_equal(summary$n_flux_inconsistent_reactions, 1L)
  expect_identical(summary$solver, solver)
  expect_equal(summary$completion_time_limit, 60)
  expect_true(file.exists(summary$file[[1L]]))
  expect_true(all(vapply(
    cache,
    function(entry) identical(entry$reaction_id, "R1"),
    logical(1)
  )))
  expect_identical(attr(cache, "completion_method"), "none")
  expect_false(attr(cache, "fastcore_executed"))
  expect_false(attr(cache, "corda2_executed"))
})

test_that("full GEM completion contract is never labelled FASTCORE", {
  answer <- list(
    model_mode = "full_gem",
    params = list(),
    model_cache_summary = data.frame(
      medium_scenario = "base",
      condition = "all",
      n_input_reactions = 4L,
      n_reactions = 3L,
      n_flux_inconsistent_reactions = 1L,
      flux_consistency_epsilon = 1e-8,
      build_strategy = "medium_flux_consistency_pruned_full_gem",
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
  expect_identical(
    result$params$structural_completion,
    "medium_flux_consistency"
  )
  expect_identical(result$completion_contract$model_completion, "none")
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

test_that("unpruned full GEM remains available to reconstruction parents", {
  gem <- .full_gem_test_model()
  model <- rc_build_full_gem(gem)
  expect_setequal(colnames(model$S), colnames(gem$S))
  expect_identical(model$build_params$strategy, "full_gem")
  expect_false(model$build_params$flux_consistency_pruning)
})
