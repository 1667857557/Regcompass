test_that("merged meta-module model scores all units on a shared GEM", {
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
    sample_id = "S1",
    module_id = "S1::M1",
    reaction_id = "R1",
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  layer1 <- list(
    C_rel = matrix(
      1,
      nrow = 3,
      ncol = 2,
      dimnames = list(colnames(S), c("u1", "u2"))
    ),
    reaction_confidence = matrix(
      1,
      nrow = 3,
      ncol = 2,
      dimnames = list(colnames(S), c("u1", "u2"))
    ),
    gpr_diagnostics = NULL,
    unit_meta = data.frame(
      pool_id = c("u1", "u2"),
      unit_id = c("u1", "u2"),
      sample_id = c("S1", "S2"),
      condition = c("ctrl", "ctrl"),
      cell_type = c("T", "T"),
      stringsAsFactors = FALSE
    )
  )
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
  )
  expect_true(result$params$merge_sample_modules)
  expect_true(result$evaluated[1L, "u1"])
  expect_true(result$evaluated[1L, "u2"])
  expect_true(result$feasible[1L, "u1"])
  expect_true(result$feasible[1L, "u2"])
  expect_equal(nrow(result$model_cache_summary), 1L)
  expect_equal(unique(result$model_cache_summary$sample_id), "__merged_meta_module__")
})

test_that("sample-specific meta-module matching remains available by opt-out", {
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
  S <- matrix(c(1, -1, 0, 0, 1, -1), nrow = 2, byrow = TRUE,
              dimnames = list(c("A_c", "B_c"), c("EX_A", "R1", "EX_B")))
  gem <- rc_make_gem(
    S, lb = c(EX_A = 0, R1 = 0, EX_B = 0),
    ub = c(EX_A = 1000, R1 = 1000, EX_B = 1000),
    reaction_meta = data.frame(reaction_id = colnames(S),
                               role = c("exchange", "internal", "exchange"),
                               role_source = "curated", stringsAsFactors = FALSE)
  )
  membership <- data.frame(sample_id = "S1", module_id = "S1::M1",
                           reaction_id = "R1", is_core = TRUE,
                           stringsAsFactors = FALSE)
  layer1 <- list(
    C_rel = matrix(1, nrow = 3, ncol = 2, dimnames = list(colnames(S), c("u1", "u2"))),
    reaction_confidence = matrix(1, nrow = 3, ncol = 2, dimnames = list(colnames(S), c("u1", "u2"))),
    gpr_diagnostics = NULL,
    unit_meta = data.frame(pool_id = c("u1", "u2"), unit_id = c("u1", "u2"),
                           sample_id = c("S1", "S2"), condition = c("ctrl", "ctrl"),
                           cell_type = c("T", "T"), stringsAsFactors = FALSE)
  )
  result <- rc_run_microcompass(
    layer1 = layer1, gem = gem, target_reactions = "R1",
    mode = "meta_module_gem", reaction_membership = membership,
    core_reactions = membership, unit = "metacell",
    target_direction = "forward", parallel = FALSE, solver = solver,
    model_params = list(merge_sample_modules = FALSE)
  )
  expect_false(result$params$merge_sample_modules)
  expect_true(result$evaluated[1L, "u1"])
  expect_false(result$evaluated[1L, "u2"])
})

test_that("v1.3 exporter writes model and LP diagnostics", {
  output_dir <- tempfile("regcompass_export_")
  result <- list(
    score = matrix(1, nrow = 1, dimnames = list("R1", "u1")),
    penalty = matrix(1, nrow = 1, dimnames = list("R1", "u1")),
    vmax = matrix(1, nrow = 1, dimnames = list("R1", "u1")),
    feasible = matrix(TRUE, nrow = 1, dimnames = list("R1", "u1")),
    evaluated = matrix(TRUE, nrow = 1, dimnames = list("R1", "u1")),
    penalty_components = list(),
    params = list(model_mode = "meta_module_gem"),
    medium_scenarios = data.frame(
      medium_scenario_id = "base",
      stringsAsFactors = FALSE
    ),
    model_cache_summary = data.frame(
      sample_id = "S1",
      module_id = "M1",
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
