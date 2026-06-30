test_that("rc_toy_gem validates a minimal GEM", {
  model <- rc_toy_gem()
  expect_equal(dim(model$S), c(2L, 4L))
  expect_equal(model$reaction_id, c("EX_glc", "ATPM", "BIOMASS", "DM_lac"))
  expect_true(all(model$lb <= model$ub))
})

test_that("rc_build_qp creates OSQP matrices with mass balance and bounds", {
  model <- rc_toy_gem()
  qp <- rc_build_baseline_qp(model, penalty = c(1, 2, 3, 4), lambda = 1e-4, atpm_rxn = "ATPM", atpm_min = 1)

  expect_equal(dim(qp$P), c(4L, 4L))
  expect_equal(length(qp$q), 4L)
  expect_equal(dim(qp$A), c(2L + 4L + 1L, 4L))
  expect_equal(length(qp$l), nrow(qp$A))
  expect_equal(length(qp$u), nrow(qp$A))
  expect_equal(qp$l[seq_len(model$S@Dim[1])], c(0, 0))
  expect_equal(qp$reaction_id, model$reaction_id)
})

test_that("rc_demand_qp adds a selected reaction demand constraint", {
  model <- rc_toy_gem()
  qp <- rc_build_baseline_qp(model, penalty = rep(1, 4))
  demanded <- rc_demand_qp(qp, "BIOMASS", delta = 2)

  expect_equal(nrow(demanded$A), nrow(qp$A) + 1L)
  expect_equal(tail(demanded$l, 1), 2)
  expect_true(is.infinite(tail(demanded$u, 1)))
  expect_equal(demanded$demand_reaction, "BIOMASS")
})

test_that("rc_solve_qp returns an OSQP status for toy baseline when rosqp is installed", {
  skip_if_not_installed("rosqp")

  model <- rc_toy_gem()
  qp <- rc_build_baseline_qp(model, penalty = rep(1, 4), atpm_rxn = "ATPM", atpm_min = 1)
  sol <- rc_solve_qp(qp, settings = list(verbose = FALSE))
  expect_true(nzchar(rc_osqp_status(sol)))
  expect_equal(names(sol$x), model$reaction_id)
})

test_that("rc_solve_selected_demand_qp solves selected toy demands when rosqp is installed", {
  skip_if_not_installed("rosqp")

  model <- rc_toy_gem()
  qp <- rc_build_baseline_qp(model, penalty = rep(1, 4))
  out <- rc_solve_selected_demand_qp(qp, reactions = c("ATPM", "BIOMASS"), delta = c(1, 1), settings = list(verbose = FALSE))
  expect_equal(out$reaction_id, c("ATPM", "BIOMASS"))
  expect_true(all(nzchar(out$osqp_status)))
})
