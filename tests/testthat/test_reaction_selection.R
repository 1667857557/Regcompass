test_that("target reaction selection honors require_complete_gpr", {
  layer1 <- list(
    C_rel = matrix(c(0.9, 0.8), nrow = 2, dimnames = list(c("R_supported", "R_unsupported"), "mc1")),
    reaction_confidence = data.frame(
      reaction_id = c("R_supported", "R_unsupported"),
      pool_id = "mc1",
      reaction_confidence = c(0.9, 0.9),
      reaction_unsupported_by_complete_gpr_flag = c(FALSE, TRUE),
      stringsAsFactors = FALSE
    ),
    unit_meta = data.frame(pool_id = "mc1", condition = "ctrl", cell_type = "T", stringsAsFactors = FALSE)
  )
  strict <- rc_select_target_reactions(layer1, require_complete_gpr = TRUE, min_units_per_group = 1, top_n = 10)
  permissive <- rc_select_target_reactions(layer1, require_complete_gpr = FALSE, min_units_per_group = 1, top_n = 10)
  expect_equal(strict$reaction_id, "R_supported")
  expect_true("R_unsupported" %in% permissive$reaction_id)
})
