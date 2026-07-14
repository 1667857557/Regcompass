rc_fastcore_solver_available <- function() {
  requireNamespace("highs", quietly = TRUE) ||
    requireNamespace("Rglpk", quietly = TRUE) ||
    requireNamespace("gurobi", quietly = TRUE)
}

rc_fastcore_test_solver <- function() {
  if (requireNamespace("highs", quietly = TRUE)) return("highs")
  if (requireNamespace("Rglpk", quietly = TRUE)) return("glpk")
  if (requireNamespace("gurobi", quietly = TRUE)) return("gurobi")
  "highs"
}

rc_fastcore_forward_toy <- function() {
  S <- matrix(
    c(
       1, -1, -1,  0,
       0,  1,  1, -1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("A_c", "B_c"),
      c("EX_A", "Rcore", "Ralt", "EX_B")
    )
  )
  rc_make_gem(
    S,
    lb = c(EX_A = 0, Rcore = 0, Ralt = 0, EX_B = 0),
    ub = c(EX_A = 1000, Rcore = 1000, Ralt = 1000, EX_B = 1000),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal", "internal", "exchange"),
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
}

rc_fastcore_reverse_toy <- function() {
  S <- matrix(
    c(
       0, -1, -1,
       1,  1,  0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("A_c", "B_c"),
      c("EX_B", "Rrev", "EX_A")
    )
  )
  rc_make_gem(
    S,
    lb = c(EX_B = 0, Rrev = -1000, EX_A = 0),
    ub = c(EX_B = 1000, Rrev = 0, EX_A = 1000),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal", "exchange"),
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
}

test_that("reaction directions are derived only from signed bounds", {
  S <- matrix(
    0,
    nrow = 1,
    ncol = 4,
    dimnames = list("dummy", c("F", "R", "B", "Z"))
  )
  gem <- rc_make_gem(
    S,
    lb = c(F = 0, R = -10, B = -10, Z = 0),
    ub = c(F = 10, R = 0, B = 10, Z = 0)
  )
  directions <- rc_prepare_directional_targets(
    gem,
    colnames(S),
    "both"
  )
  expect_setequal(
    paste(
      directions$reaction_id,
      directions$target_direction,
      sep = ":"
    ),
    c("F:forward", "R:reverse", "B:forward", "B:reverse")
  )
  expect_false("Z" %in% directions$reaction_id)
})

test_that("signed absolute-penalty LP preserves forced non-zero bounds", {
  skip_if_not(rc_fastcore_solver_available())
  solver <- rc_fastcore_test_solver()
  S <- matrix(
    0,
    nrow = 1,
    ncol = 1,
    dimnames = list("dummy", "R")
  )

  forward <- rc_compass_two_step_lp_directional(
    S,
    lb = c(R = 5),
    ub = c(R = 7),
    target_reaction = "R",
    penalties = c(R = 1),
    target_direction = "forward",
    omega = 0.5,
    solver = solver
  )
  expect_true(forward$feasible)
  expect_equal(forward$vmax, 7, tolerance = 1e-7)
  expect_equal(forward$penalty, 5, tolerance = 1e-7)

  reverse <- rc_compass_two_step_lp_directional(
    S,
    lb = c(R = -9),
    ub = c(R = -5),
    target_reaction = "R",
    penalties = c(R = 1),
    target_direction = "reverse",
    omega = 0.5,
    solver = solver
  )
  expect_true(reverse$feasible)
  expect_equal(reverse$vmax, 9, tolerance = 1e-7)
  expect_equal(reverse$penalty, 5, tolerance = 1e-7)
})

test_that("add-only FASTCORE preserves the full biological set", {
  skip_if_not(rc_fastcore_solver_available())
  gem <- rc_fastcore_forward_toy()
  membership <- data.frame(
    sample_id = "S1",
    module_id = "S1::GRN0001",
    reaction_id = c("Rcore", "Ralt"),
    is_core = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  core <- membership[membership$is_core, , drop = FALSE]
  model <- rc_build_meta_module_gem(
    gem,
    membership,
    core,
    sample_id = "S1",
    module_id = "S1::GRN0001",
    solver = rc_fastcore_test_solver(),
    strict = TRUE
  )

  expect_true(all(c("Rcore", "Ralt") %in% colnames(model$S)))
  expect_setequal(
    model$reaction_meta$reaction_id[
      model$reaction_meta$fastcore_support
    ],
    c("EX_A", "EX_B")
  )
  expect_true(all(
    model$closure_diagnostics$final_feasible[
      model$closure_diagnostics$feasible
    ]
  ))
  expect_equal(
    model$build_params$algorithm,
    "add_only_fastcore_lp7_lp10"
  )
})

test_that("reverse-only core reactions are completed by orientation", {
  skip_if_not(rc_fastcore_solver_available())
  gem <- rc_fastcore_reverse_toy()
  membership <- data.frame(
    sample_id = "S1",
    module_id = "S1::GRN0002",
    reaction_id = "Rrev",
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  model <- rc_build_meta_module_gem(
    gem,
    membership,
    membership,
    sample_id = "S1",
    module_id = "S1::GRN0002",
    solver = rc_fastcore_test_solver(),
    strict = TRUE
  )
  expect_equal(
    model$target_directions$target_direction,
    "reverse"
  )
  expect_true(model$closure_diagnostics$final_feasible)
  expect_setequal(
    model$reaction_meta$reaction_id[
      model$reaction_meta$fastcore_support
    ],
    c("EX_A", "EX_B")
  )
})

test_that("parent-blocked core reactions are not gap-filled", {
  skip_if_not(rc_fastcore_solver_available())
  S <- matrix(
    1,
    nrow = 1,
    dimnames = list("A_c", "Rblocked")
  )
  gem <- rc_make_gem(
    S,
    lb = c(Rblocked = 0),
    ub = c(Rblocked = 1000),
    reaction_meta = data.frame(
      reaction_id = "Rblocked",
      role = "internal",
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
  membership <- data.frame(
    sample_id = "S1",
    module_id = "M1",
    reaction_id = "Rblocked",
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  model <- rc_build_meta_module_gem(
    gem,
    membership,
    membership,
    sample_id = "S1",
    module_id = "M1",
    solver = rc_fastcore_test_solver(),
    strict = TRUE
  )
  expect_equal(
    model$closure_diagnostics$completion_status,
    "parent_blocked"
  )
  expect_false(any(model$reaction_meta$fastcore_support))
})

test_that("LP-10 scaling retains smaller stoichiometric support flux", {
  skip_if_not(rc_fastcore_solver_available())
  S <- matrix(
    c(
       2, -1,  0,
       0,  1, -1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("A_c", "B_c"),
      c("EX_2A", "Rcore", "EX_B")
    )
  )
  gem <- rc_make_gem(
    S,
    lb = c(EX_2A = 0, Rcore = 0, EX_B = 0),
    ub = c(EX_2A = 1000, Rcore = 1000, EX_B = 1000),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal", "exchange"),
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
  membership <- data.frame(
    sample_id = "S1",
    module_id = "Mscale",
    reaction_id = "Rcore",
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  model <- rc_build_meta_module_gem(
    gem,
    membership,
    membership,
    sample_id = "S1",
    module_id = "Mscale",
    solver = rc_fastcore_test_solver(),
    fastcore_epsilon = 1e-4,
    strict = TRUE
  )
  expect_setequal(
    model$reaction_meta$reaction_id[
      model$reaction_meta$fastcore_support
    ],
    c("EX_2A", "EX_B")
  )
  expect_true(all(
    model$completion_iterations$lp10_scaling_factor == 1e5
  ))
})

test_that("FASTCC does not create a false reversible split cycle", {
  skip_if_not(rc_fastcore_solver_available())
  S <- matrix(
    1,
    nrow = 1,
    ncol = 1,
    dimnames = list("A_c", "Rrev_blocked")
  )
  gem <- rc_make_gem(
    S,
    lb = c(Rrev_blocked = -1000),
    ub = c(Rrev_blocked = 1000),
    reaction_meta = data.frame(
      reaction_id = "Rrev_blocked",
      role = "internal",
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
  consistent <- .rc_fastcc_consistent_reactions(
    gem,
    solver = rc_fastcore_test_solver(),
    epsilon = 1e-4
  )
  expect_false("Rrev_blocked" %in% consistent)
})

test_that("one completed GEM is cached per sample module and medium", {
  skip_if_not(rc_fastcore_solver_available())
  gem <- rc_fastcore_forward_toy()
  membership <- data.frame(
    sample_id = "S1",
    module_id = "S1::GRN0001",
    reaction_id = "Rcore",
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  cache <- rc_build_meta_module_gem_cache(
    gem = gem,
    reaction_membership = membership,
    core_reactions = membership,
    medium_scenarios = NULL,
    solver = rc_fastcore_test_solver()
  )
  summary <- attr(cache, "summary")
  expect_equal(nrow(summary), 1)
  expect_equal(summary$build_strategy, "meta_module_gem")
  expect_true(file.exists(summary$file[[1L]]))
  expect_equal(
    length(unique(vapply(cache, `[[`, character(1), "file"))),
    1
  )
})

test_that("the Layer 2 API exposes exactly two structural modes", {
  expect_identical(
    eval(formals(rc_run_microcompass)$mode),
    c("full_gem", "meta_module_gem")
  )
  expect_error(
    rc_run_microcompass(
      layer1 = list(),
      gem = list(),
      mode = "module_meso_gem"
    ),
    "should be one of"
  )
  exports <- getNamespaceExports("RegCompassR")
  expect_false(any(c(
    "rc_build_target_microgem",
    "rc_build_microgem_cache",
    "rc_build_module_meso_gem",
    "rc_build_module_gem_cache",
    "rc_run_regcompass_v12"
  ) %in% exports))
})

test_that("v1.3 labeled row IDs retain sample and module", {
  parsed <- rc_parse_microcompass_row_id(
    "sample=S1::module=S1%3A%3AM1::reaction=R1::direction=forward::medium=base"
  )
  expect_equal(parsed$sample_id, "S1")
  expect_equal(parsed$module_id, "S1::M1")
  expect_equal(parsed$reaction_id, "R1")
  expect_equal(parsed$target_direction, "forward")
  expect_equal(parsed$medium_scenario, "base")
})
