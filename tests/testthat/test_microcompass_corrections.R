test_that("reverse vmax uses minimization objective convention", {
  S <- matrix(0, nrow = 0, ncol = 1, dimnames = list(NULL, "Rrev"))
  ans <- rc_compass_vmax_directional(S, lb = c(Rrev = -3), ub = c(Rrev = 2), target_reaction = "Rrev", direction = "reverse")
  if (identical(ans$status, "error")) skip("No LP solver available")
  expect_true(ans$feasible)
  expect_equal(ans$vmax, 3, tolerance = 1e-6)
})

test_that("curated support penalty overrides instead of adding to evidence penalty", {
  C <- matrix(0.01, nrow = 2, ncol = 1, dimnames = list(c("EX_curated", "EX_inferred"), "u1"))
  Conf <- C
  roles <- data.frame(
    reaction_id = c("EX_curated", "EX_inferred"),
    role = "exchange",
    role_source = c("curated", "id_pattern"),
    stringsAsFactors = FALSE
  )
  out <- rc_compute_multiome_penalty(C, Conf, reaction_roles = roles)
  expect_equal(out$penalty["EX_curated", "u1"], 0.05)
  expect_gt(out$penalty["EX_inferred", "u1"], 0.05)
  expect_true(out$components$role_override_flag["EX_curated"])
  expect_false(out$components$role_override_flag["EX_inferred"])
})

test_that("single-stoichiometry reactions are boundary_like, not exchange", {
  S <- matrix(c(1, 0, 1, -1), nrow = 2, dimnames = list(c("m1", "m2"), c("R_boundary", "R_internal")))
  gem <- rc_annotate_reaction_roles(rc_make_gem(S, lb = c(0, 0), ub = c(1000, 1000)), infer_from_id = FALSE, infer_from_compartment = FALSE)
  roles <- stats::setNames(gem$reaction_roles$role, gem$reaction_roles$reaction_id)
  expect_equal(roles[["R_boundary"]], "boundary_like")
})

test_that("microCOMPASS row IDs parse medium scenario explicitly", {
  parsed <- rc_parse_microcompass_row_id(c("R1::forward::medium=blood_like", "R2::reverse"))
  expect_equal(parsed$reaction_id, c("R1", "R2"))
  expect_equal(parsed$target_direction, c("forward", "reverse"))
  expect_equal(parsed$medium_scenario[1], "blood_like")
  expect_true(is.na(parsed$medium_scenario[2]))
})

test_that("relaxed LP uses fallback target minimum when strict infeasible", {
  S <- matrix(c(1, 0), nrow = 1, dimnames = list("m1", c("Rtarget", "Rdummy")))
  mg <- rc_make_gem(S, lb = c(Rtarget = 0, Rdummy = 0), ub = c(Rtarget = 0, Rdummy = 0))
  out <- rc_run_relaxed_balance_lp(mg, penalties = c(Rtarget = 1, Rdummy = 1), target_reaction = "Rtarget", target_min_if_strict_infeasible = 0.01)
  if (identical(out$solver_status, "error")) skip("No LP solver available")
  expect_false(out$strict_feasible)
  expect_equal(out$target_min_used, 0.01)
  expect_equal(out$target_min_if_strict_infeasible, 0.01)
})
