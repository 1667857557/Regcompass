test_that("micro-GEM auto strategy falls back from target khop to module meso-GEM", {
  S <- matrix(c(1, -1, 0,
                0,  1, -1), nrow = 2, byrow = TRUE,
              dimnames = list(c("m1", "m2"), c("EX_m1", "Rtarget", "DM_m2")))
  reaction_meta <- data.frame(
    reaction_id = colnames(S),
    role = c("exchange", "internal", "demand"),
    metabolic_module = "module_a",
    stringsAsFactors = FALSE
  )
  gem <- rc_make_gem(S, lb = c(EX_m1 = -10, Rtarget = 0, DM_m2 = 0), ub = c(EX_m1 = 1000, Rtarget = 1000, DM_m2 = 1000), reaction_meta = reaction_meta)
  dirs <- data.frame(reaction_id = "Rtarget", target_direction = "forward", stringsAsFactors = FALSE)
  medium <- data.frame(medium_scenario_id = "base", exchange_reaction_id = character(), lb = numeric(), ub = numeric(), available = logical(), stringsAsFactors = FALSE)

  cache <- rc_build_microgem_cache(gem, dirs, medium, microgem_params = list(strategy = "auto", k_hop = 0, module_col = "metabolic_module"))
  mg <- cache[[1]]
  expect_equal(mg$build_params$fallback_from, "target_khop")
  expect_equal(mg$target_status, "ok")
  expect_true(isTRUE(mg$closure_diagnostics$strict_target_feasible[[1]]))
})

test_that("strict_closure is enforced for structurally infeasible target micro-GEMs", {
  S <- matrix(c(1), nrow = 1, dimnames = list("m1", "Rtarget"))
  gem <- rc_make_gem(S, lb = c(Rtarget = 0), ub = c(Rtarget = 1000), reaction_meta = data.frame(reaction_id = "Rtarget", role = "internal", stringsAsFactors = FALSE))
  expect_error(rc_build_target_microgem(gem, "Rtarget", k_hop = 0, strict_closure = TRUE), "strict closure")
})

test_that("HiGHS numeric status 7 with valid optimal fields is parsed as optimal", {
  res <- list(status = 7L, model_status = "Optimal", objective_value = 1, primal_solution = c(1, 0))
  expect_equal(.rc_highs_status(res, n_variables = 2), "optimal")
})

test_that("infeasible COMPASS scores are NA rather than zero", {
  P <- matrix(c(1, 2), nrow = 1, dimnames = list("R1", c("u1", "u2")))
  feasible <- matrix(c(TRUE, FALSE), nrow = 1, dimnames = dimnames(P))
  score <- rc_compass_score_from_penalty(P, feasible)
  expect_true(is.na(score["R1", "u2"]))
})

test_that("single-sample differential output is descriptive-only when strict design is disabled", {
  result <- list(
    score = matrix(c(0.2, 0.8), nrow = 1, dimnames = list("R1::forward::base", c("u1", "u2"))),
    unit_meta = data.frame(unit_id = c("u1", "u2"), sample_id = "s1", condition = "ctrl", cell_type = "T", n_cells = c(10, 20), stringsAsFactors = FALSE)
  )
  expect_error(rc_test_microcompass_differential(result, strict_replicate_design = TRUE), "biological samples")
  out <- rc_test_microcompass_differential(result, strict_replicate_design = FALSE)
  expect_true(all(is.na(out$p_value)))
  expect_true(all(is.na(out$FDR)))
  expect_equal(unique(out$model_status), "descriptive_only")
  expect_equal(unique(out$n_biological_samples), 1L)
})

test_that("multi-class condition lm uses omnibus test instead of first coefficient", {
  result <- list(
    score = matrix(c(1, 2, 3, 1.2, 2.2, 3.2), nrow = 1, dimnames = list("R1::forward::base", paste0("u", 1:6))),
    unit_meta = data.frame(unit_id = paste0("u", 1:6), sample_id = paste0("s", 1:6), condition = rep(c("A", "B", "C"), each = 2), cell_type = "T", stringsAsFactors = FALSE)
  )
  out <- rc_test_microcompass_differential(result, method = "lm", min_samples_per_group = 2, strict_replicate_design = TRUE, test_type = "omnibus")
  expect_equal(out$contrast, "condition_omnibus")
  expect_true(is.finite(out$p_value))
})
