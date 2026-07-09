rc_test_toy_gem <- function() {
  S <- matrix(c(1, -1, 0,
                0,  1, -1), nrow = 2, byrow = TRUE,
              dimnames = list(c("m1", "m2"), c("EX_m1", "Rtarget", "DM_m2")))
  meta <- data.frame(reaction_id = colnames(S), role = c("exchange", "internal", "demand"), role_source = "curated", stringsAsFactors = FALSE)
  rc_make_gem(S, lb = c(EX_m1 = -10, Rtarget = 0, DM_m2 = 0), ub = c(EX_m1 = 1000, Rtarget = 1000, DM_m2 = 1000), reaction_meta = meta)
}

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
  roles <- data.frame(reaction_id = c("EX_curated", "EX_inferred"), role = "exchange", role_source = c("curated", "id_pattern"), stringsAsFactors = FALSE)
  out <- rc_compute_multiome_penalty(C, Conf, reaction_roles = roles)
  expect_equal(out$evidence_policy, "penalty_only")
  expect_equal(out$penalty["EX_curated", "u1"], 0.05)
  expect_gt(out$penalty["EX_inferred", "u1"], 0.05)
})

test_that("single-stoichiometry reactions are boundary_like, not exchange", {
  S <- matrix(c(1, 0, 1, -1), nrow = 2, dimnames = list(c("m1", "m2"), c("R_boundary", "R_internal")))
  gem <- rc_annotate_reaction_roles(rc_make_gem(S, lb = c(0, 0), ub = c(1000, 1000)), infer_from_id = FALSE, infer_from_compartment = FALSE)
  roles <- stats::setNames(gem$reaction_roles$role, gem$reaction_roles$reaction_id)
  expect_equal(roles[["R_boundary"]], "boundary_like")
})

test_that("microCOMPASS row IDs parse medium scenario explicitly", {
  parsed <- rc_parse_microcompass_row_id(c("R1::forward::blood_like", "R2::reverse::medium=old"))
  expect_equal(parsed$reaction_id, c("R1", "R2"))
  expect_equal(parsed$target_direction, c("forward", "reverse"))
  expect_equal(parsed$medium_scenario, c("blood_like", "old"))
})

test_that("structural micro-GEM cache is evidence-independent", {
  gem <- rc_test_toy_gem()
  dirs <- data.frame(reaction_id = "Rtarget", target_direction = "forward", stringsAsFactors = FALSE)
  medium <- data.frame(medium_scenario_id = "base", exchange_reaction_id = character(), lb = numeric(), ub = numeric(), available = logical())
  cache1 <- rc_build_microgem_cache(gem, dirs, medium)
  cache2 <- rc_build_microgem_cache(gem, dirs, medium)
  expect_identical(lapply(cache1, function(x) colnames(x$S)), lapply(cache2, function(x) colnames(x$S)))
})

test_that("penalty is unit-specific", {
  C <- matrix(c(0.9, 0.1), nrow = 1, dimnames = list("Rtarget", c("u1", "u2")))
  Conf <- matrix(1, nrow = 1, ncol = 2, dimnames = dimnames(C))
  pen <- rc_compute_multiome_penalty(C, Conf)
  expect_false(identical(pen$penalty[, 1], pen$penalty[, 2]))
})

test_that("Human2 GEM requires model_info and gpr_table", {
  gem <- rc_test_toy_gem()
  expect_error(rc_validate_human2_gem(gem), "model_info")
  gem$model_info <- list(source = "Human-GEM", version = "2.0.0", commit = "abc", checksum = "sha", conversion_date = "2026-01-01")
  expect_error(rc_validate_human2_gem(gem), "gpr_table")
  gem$gpr_table <- data.frame(reaction_id = "Rtarget", and_group_id = 1L, gene = "G", stringsAsFactors = FALSE)
  expect_silent(rc_validate_human2_gem(gem))
})


test_that("parallel and serial microCOMPASS agree with structural cache", {
  gem <- rc_test_toy_gem()
  layer1 <- list(
    C_rel = matrix(c(0.9, 0.2, 0.8, 0.8, 0.7, 0.3), nrow = 3,
                   dimnames = list(c("EX_m1", "Rtarget", "DM_m2"), c("p1", "p2"))),
    reaction_confidence = matrix(1, nrow = 3, ncol = 2,
                                 dimnames = list(c("EX_m1", "Rtarget", "DM_m2"), c("p1", "p2"))),
    gpr_diagnostics = NULL,
    unit_meta = data.frame(pool_id = c("p1", "p2"), sample_id = c("s1", "s2"), condition = c("a", "b"), cell_type = "T", stringsAsFactors = FALSE)
  )
  medium <- data.frame(medium_scenario_id = "base", exchange_reaction_id = character(), lb = numeric(), ub = numeric(), available = logical())
  res1 <- rc_run_microcompass(layer1, gem, "Rtarget", medium_scenarios = medium, unit = "metacell", target_direction = "forward", parallel = FALSE)
  res2 <- rc_run_microcompass(layer1, gem, "Rtarget", medium_scenarios = medium, unit = "metacell", target_direction = "forward", parallel = TRUE, BPPARAM = FALSE)
  expect_equal(res1$score, res2$score, tolerance = 1e-8)
  expect_equal(res1$penalty, res2$penalty, tolerance = 1e-8)
  expect_equal(res1$feasible, res2$feasible)
})
