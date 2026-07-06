test_that("Layer 2 penalty mapping penalizes missing evidence and exempts support reactions", {
  C <- matrix(c(1, NA, 0.1), nrow = 3, dimnames = list(c("R_high", "R_missing", "EX_a"), "u1"))
  Conf <- matrix(c(1, NA, 0.1), nrow = 3, dimnames = list(c("R_high", "R_missing", "EX_a"), "u1"))

  out <- rc_layer2_penalty(C, Conf, support_reactions = "EX_a")

  expect_lt(out$penalty["R_high", "u1"], out$penalty["R_missing", "u1"])
  expect_equal(out$penalty["EX_a", "u1"], 0)
  expect_true(out$components$missing_evidence["R_missing", "u1"])
})

test_that("Layer 2 accepts long-form Layer 1 reaction confidence", {
  C <- matrix(c(0.7, 0.2), nrow = 2, dimnames = list(c("R1", "R2"), "pool1"))
  conf <- data.frame(
    reaction_id = c("R1", "R2"),
    pool_id = "pool1",
    reaction_confidence = c(0.8, 0.3),
    stringsAsFactors = FALSE
  )

  M <- rc_layer2_confidence_matrix(conf, C)

  expect_equal(M["R1", "pool1"], 0.8)
  expect_equal(M["R2", "pool1"], 0.3)
})

test_that("Layer 2 reaction selection excludes all-missing, unsupported, and very-low Q95 reactions", {
  layer1 <- list(
    C_rel = matrix(c(0.9, NA, 0.8, 0.7), nrow = 4, dimnames = list(c("R_keep", "R_missing", "R_unsupported", "R_lowq"), "pool1")),
    reaction_confidence = data.frame(
      reaction_id = c("R_keep", "R_missing", "R_unsupported", "R_lowq"),
      pool_id = "pool1",
      reaction_confidence = c(0.9, NA, NA, 0.9),
      reaction_unsupported_by_complete_gpr_flag = c(FALSE, FALSE, TRUE, FALSE),
      stringsAsFactors = FALSE
    ),
    q95_diagnostics = data.frame(
      reaction_id = c("R_keep", "R_missing", "R_unsupported", "R_lowq"),
      all_missing_reaction_flag = c(FALSE, TRUE, FALSE, FALSE),
      q95_power_class = c("adequate", "adequate", "adequate", "very_low"),
      stringsAsFactors = FALSE
    )
  )
  S <- matrix(1, nrow = 1, dimnames = list("m1", rownames(layer1$C_rel)))
  gem <- rc_make_gem(S, lb = rep(0, 4), ub = rep(1000, 4))

  sel <- rc_select_layer2_reactions(layer1, gem, neighbor_depth = 0)

  expect_equal(sel$reaction_id, "R_keep")
})

test_that("absolute penalty LP builder represents target minimum and variable bounds", {
  S <- matrix(c(1, -1), nrow = 1, dimnames = list("m1", c("R1", "R2")))
  lp <- rc_build_abs_penalty_lp(S, lb = c(R1 = -10, R2 = 0), ub = c(R1 = 5, R2 = 7), penalties = c(R1 = 2, R2 = 3), target_index = 1, target_min = 4)

  expect_equal(length(lp$obj), 4)
  expect_equal(lp$lhs[2], 4)
  expect_true(is.infinite(lp$rhs[2]))
  expect_equal(lp$ub, c(5, 7, 10, 0))
})
