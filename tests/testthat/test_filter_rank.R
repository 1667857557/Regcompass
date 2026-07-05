test_that("valid reaction filtering excludes missing unsupported low-power reactions", {
  layer1 <- list(
    C_rel = rbind(R_good = c(0.5, 0.8), R_bad = c(NA_real_, NA_real_)),
    C_raw = rbind(R_good = c(1, 2), R_bad = c(NA_real_, NA_real_)),
    q95_diagnostics = data.frame(
      reaction_id = c("R_good", "R_bad"),
      all_missing_reaction_flag = c(FALSE, TRUE),
      q95_power_class = c("adequate", "very_low")
    ),
    reaction_confidence = data.frame(
      reaction_id = c("R_good", "R_bad"),
      reaction_confidence = c(0.8, NA_real_),
      reaction_unsupported_by_complete_gpr_flag = c(FALSE, TRUE)
    )
  )
  out <- rc_filter_valid_reactions(layer1)
  expect_equal(out$valid_reactions, "R_good")
  ranked <- rc_rank_reactions(layer1)
  expect_equal(ranked$reaction_id, "R_good")
})
