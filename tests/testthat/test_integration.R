test_that("global meta-module model scores every metacell", {
  skip_if_not(
    requireNamespace("highs", quietly = TRUE) ||
      requireNamespace("Rglpk", quietly = TRUE) ||
      requireNamespace("gurobi", quietly = TRUE)
  )
  solver <- if (requireNamespace("highs", quietly = TRUE)) {
    "highs"
  } else if (requireNamespace("Rglpk", quietly = TRUE)) {
    "glpk"
  } else {
    "gurobi"
  }
  S <- matrix(
    c(
       1, -1,  0,
       0,  1, -1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("A_c", "B_c"),
      c("EX_A", "R1", "EX_B")
    )
  )
  gem <- rc_make_gem(
    S,
    lb = c(EX_A = 0, R1 = 0, EX_B = 0),
    ub = c(EX_A = 1000, R1 = 1000, EX_B = 1000),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal", "exchange"),
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
  membership <- data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    reaction_id = "R1",
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  layer1 <- list(
    reaction_expression = matrix(
      1,
      nrow = 3,
      ncol = 2,
      dimnames = list(colnames(S), c("u1", "u2"))
    ),
    unit_meta = data.frame(
      pool_id = c("u1", "u2"),
      unit_id = c("u1", "u2"),
      sample_id = c("S1", "S2"),
      condition = c("ctrl", "ctrl"),
      cell_type = c("T", "T"),
      stringsAsFactors = FALSE
    )
  )
  result <- NULL
  expect_warning(
    result <- rc_run_microcompass(
      layer1 = layer1,
      gem = gem,
      target_reactions = "R1",
      mode = "meta_module_gem",
      reaction_membership = membership,
      core_reactions = membership,
      unit = "metacell",
      target_direction = "forward",
      parallel = FALSE,
      solver = solver
    ),
    "descriptive pseudo-observations"
  )
  expect_true(all(result$evaluated[1L, ]))
  expect_true(all(result$feasible[1L, ]))
  expect_true(all(is.finite(result$penalty[1L, ])))
  expect_true(all(is.na(result$score[1L, ])))
  expect_true(result$noninformative_target[[1L]])
  expect_equal(nrow(result$model_cache_summary), 1L)
  expect_true(result$params$shared_gem)
  expect_equal(result$params$shared_gem_scope, "all_metacells")
})

test_that("condition-specific medium is rejected for shared-GEM scoring", {
  layer1 <- list(
    reaction_expression = matrix(
      1,
      nrow = 1,
      dimnames = list("R1", "u1")
    ),
    unit_meta = data.frame(
      pool_id = "u1", unit_id = "u1", sample_id = "S1",
      condition = "A", cell_type = "T", stringsAsFactors = FALSE
    )
  )
  gem <- rc_make_gem(
    matrix(0, nrow = 1, dimnames = list("m", "R1")),
    lb = c(R1 = 0),
    ub = c(R1 = 1)
  )
  medium <- data.frame(
    medium_scenario_id = "custom", exchange_reaction_id = "R1",
    lb = 0, ub = 1, available = TRUE, condition = "A"
  )
  expect_error(
    rc_run_microcompass(
      layer1,
      gem,
      "R1",
      medium_scenarios = medium,
      mode = "full_gem",
      parallel = FALSE
    ),
    "condition-invariant"
  )
})

test_that("v2 exporter writes model and LP diagnostics", {
  output_dir <- tempfile("regcompass_export_")
  result <- list(
    score = matrix(1, nrow = 1, dimnames = list("R1", "u1")),
    penalty = matrix(1, nrow = 1, dimnames = list("R1", "u1")),
    vmax = matrix(1, nrow = 1, dimnames = list("R1", "u1")),
    feasible = matrix(TRUE, nrow = 1, dimnames = list("R1", "u1")),
    evaluated = matrix(TRUE, nrow = 1, dimnames = list("R1", "u1")),
    penalty_components = list(),
    params = list(model_mode = "meta_module_gem", shared_gem = TRUE),
    medium_scenarios = data.frame(
      medium_scenario_id = "base",
      stringsAsFactors = FALSE
    ),
    model_cache_summary = data.frame(
      sample_id = "global",
      module_id = "GLOBAL_UNION",
      stringsAsFactors = FALSE
    ),
    model_diagnostics = data.frame(
      reaction_id = "R1",
      completion_status = "fastcore_completed",
      stringsAsFactors = FALSE
    ),
    lp_diagnostics = data.frame(
      reaction_id = "R1",
      solver_status = "optimal",
      stringsAsFactors = FALSE
    )
  )
  rc_export_microcompass(result, output_dir)
  expect_true(file.exists(file.path(
    output_dir, "03_models", "model_cache_summary.tsv.gz"
  )))
  expect_true(file.exists(file.path(
    output_dir, "03_models", "model_diagnostics.tsv.gz"
  )))
  expect_true(file.exists(file.path(
    output_dir, "04_microcompass", "lp_diagnostics.tsv.gz"
  )))
  expect_true(file.exists(file.path(
    output_dir, "04_microcompass", "evaluated_matrix.rds"
  )))
})
