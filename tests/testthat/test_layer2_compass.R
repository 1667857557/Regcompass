test_that("Layer 2 penalty mapping penalizes missing evidence and exempts support reactions", {
  C <- matrix(c(1, NA, 0.1), nrow = 3, dimnames = list(c("R_high", "R_missing", "EX_a"), "u1"))
  Conf <- matrix(c(1, NA, 0.1), nrow = 3, dimnames = list(c("R_high", "R_missing", "EX_a"), "u1"))

  out <- rc_layer2_penalty(C, Conf, support_reactions = "EX_a")

  expect_lt(out$penalty["R_high", "u1"], out$penalty["R_missing", "u1"])
  expect_equal(out$penalty["EX_a", "u1"], 0.05)
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


test_that("transport reactions are not support-exempt by default", {
  meta <- data.frame(
    reaction_id = c("EX_a", "DM_b", "T_gpr", "T_nogpr"),
    type = c("exchange", "demand", "transport", "transport"),
    gpr = c("", "", "geneA", ""),
    stringsAsFactors = FALSE
  )
  gem <- list(reaction_meta = meta)
  C <- matrix(NA_real_, nrow = 4, ncol = 1, dimnames = list(meta$reaction_id, "u1"))
  Conf <- C

  normal <- rc_layer2_support_penalties(gem, meta$reaction_id, C, Conf)
  reduced <- rc_layer2_support_penalties(gem, meta$reaction_id, C, Conf, transport_penalty_mode = "reduced", transport_reduced_penalty = 1)

  expect_true("EX_a" %in% names(normal))
  expect_true("DM_b" %in% names(normal))
  expect_false("T_gpr" %in% names(normal))
  expect_false("T_nogpr" %in% names(normal))
  expect_false("T_gpr" %in% names(reduced))
  expect_equal(reduced["T_nogpr"], c(T_nogpr = 1))
})

test_that("custom selected reactions filter invalid entries unless explicitly overridden", {
  layer1 <- list(
    C_rel = matrix(c(0.9, NA), nrow = 2, dimnames = list(c("R_keep", "R_missing"), "pool1")),
    reaction_confidence = data.frame(
      reaction_id = c("R_keep", "R_missing"),
      pool_id = "pool1",
      reaction_confidence = c(0.9, NA),
      stringsAsFactors = FALSE
    ),
    q95_diagnostics = data.frame(
      reaction_id = c("R_keep", "R_missing"),
      all_missing_reaction_flag = c(FALSE, TRUE),
      stringsAsFactors = FALSE
    )
  )
  S <- matrix(1, nrow = 1, dimnames = list("m1", rownames(layer1$C_rel)))
  gem <- rc_make_gem(S, lb = rep(0, 2), ub = rep(1000, 2))

  expect_warning(sel <- rc_select_layer2_reactions(layer1, gem, selected_reactions = c("R_keep", "R_missing"), neighbor_depth = 0), "Filtered invalid")
  override <- rc_select_layer2_reactions(layer1, gem, selected_reactions = c("R_keep", "R_missing"), neighbor_depth = 0, override_invalid = TRUE)

  expect_equal(sel$reaction_id, "R_keep")
  expect_true(grepl("R_missing", attr(sel, "custom_invalid_reaction_warning")))
  expect_setequal(override$reaction_id, c("R_keep", "R_missing"))
})
