test_that("CORDA2 options match the pinned Python constructor and constants", {
  defaults <- RegCompassR:::.rc_layer2_corda_options(list())
  expect_identical(defaults$model_completion, "fastcore")

  options <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda2"
  ))
  expect_identical(options$model_completion, "corda2")
  expect_identical(options$requested_model_completion, "corda2")
  expect_null(options$met_prod)
  expect_identical(options$n, 3L)
  expect_identical(options$redundancies, 3L)
  expect_equal(options$penalty_factor, 100)
  expect_identical(options$support, 5L)
  expect_identical(options$cost_increase, 1.01)
  expect_identical(options$target_flux, 1)
  expect_identical(options$upper_bound, 1e6)
  expect_identical(
    options$algorithm,
    "resendislab_python_CORDA2_c02e06d_exact_semantics"
  )
  expect_identical(options$supported_scope, "met_prod = NULL")

  exact <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda2",
    corda2_args = list(
      met_prod = NULL, n = 2L, penalty_factor = 10, support = 4L
    )
  ))
  expect_identical(exact$n, 2L)
  expect_equal(exact$penalty_factor, 10)
  expect_identical(exact$support, 4L)

  for (alias in c("corda", "corda_like")) {
    value <- RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = alias
    ))
    expect_identical(value$requested_model_completion, "corda2")
    expect_identical(value$algorithm, options$algorithm)
  }

  for (parameter in c(
    "corda2_cost_increase", "corda2_target_flux",
    "corda2_flux_tolerance", "corda_seed"
  )) {
    args <- list(model_completion = "corda2")
    args[[parameter]] <- 2
    expect_error(
      RegCompassR:::.rc_layer2_corda_options(args),
      "does not expose"
    )
  }
  expect_error(
    RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = "corda2",
      corda2_args = list(met_prod = "atp_c")
    )),
    "met_prod = NULL"
  )
})

test_that("CORDA2 solver tolerance follows the selected solver", {
  expect_identical(
    RegCompassR:::.rc_corda2_solver_feasibility_tolerance("highs"),
    1e-7
  )
  expect_identical(
    RegCompassR:::.rc_corda2_solver_feasibility_tolerance("glpk"),
    1e-7
  )
  expect_identical(
    RegCompassR:::.rc_corda2_solver_feasibility_tolerance("gurobi"),
    1e-6
  )
})

test_that("CORDA2 creates both COBRA direction variables and normalizes bounds", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1, 1, -1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), c("IRR", "REV"))
  )
  split <- RegCompassR:::.rc_corda_split_model(list(
    S = S,
    lb = c(IRR = 0, REV = -5),
    ub = c(IRR = 10, REV = 7)
  ), tolerance = 1e-7)
  expect_identical(
    split$direction_table$variable_id,
    c(
      "IRR::forward", "IRR::reverse",
      "REV::forward", "REV::reverse"
    )
  )
  expect_equal(split$ub[["IRR::forward"]], 1e6)
  expect_equal(split$ub[["IRR::reverse"]], 0)
  expect_equal(split$ub[["REV::forward"]], 1e6)
  expect_equal(split$ub[["REV::reverse"]], 1e6)
})

test_that("Layer 2 directly prepares and finalizes CORDA2", {
  implementation <- paste(
    deparse(body(rc_regcompass_step_layer2)), collapse = "\n"
  )
  expect_match(implementation, ".rc_prepare_corda_worker_pool", fixed = TRUE)
  expect_match(implementation, ".rc_layer2_prepare_completion", fixed = TRUE)
  expect_match(implementation, ".rc_layer2_finalize_completion", fixed = TRUE)
  expect_false(grepl("_base", implementation, fixed = TRUE))
  expect_false(grepl("before_", implementation, fixed = TRUE))
  expect_identical(
    names(formals(rc_regcompass_step_layer2)),
    c(
      "layer1", "meta_modules", "gem", "medium_scenarios", "outdir",
      "model_mode", "layer2_args", "parallel", "BPPARAM", "progress"
    )
  )
})

test_that("original union cache directly dispatches to CORDA2 helper", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_build_celltype_medium_union_gem_cache)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    ".rc_build_celltype_medium_corda_cache",
    fixed = TRUE
  )
  expect_match(implementation, "context$model_completion", fixed = TRUE)
  expect_false(grepl("_base", implementation, fixed = TRUE))

  helper <- paste(
    deparse(body(RegCompassR:::.rc_build_celltype_medium_corda_cache)),
    collapse = "\n"
  )
  expect_match(
    helper,
    "serial_within_each_python_corda2_instance",
    fixed = TRUE
  )
  expect_match(helper, "celltype_medium_python_corda2_exact", fixed = TRUE)
})
