test_that("rc_select_reactions selects variable, exchange, transport, and user reactions", {
  C_rel <- matrix(
    c(0.1, 0.1, 0.1,
      0.0, 0.5, 1.0,
      0.2, 0.3, 0.4,
      1.0, 0.0, 1.0),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(c("R_const", "R_var", "EX_glc", "T_lac"), paste0("pool", 1:3))
  )
  meta <- data.frame(
    reaction_id = rownames(C_rel),
    is_exchange = c(FALSE, FALSE, TRUE, FALSE),
    is_transport = c(FALSE, FALSE, FALSE, TRUE)
  )

  selected <- rc_select_reactions(C_rel, meta, top_n = 1, user_reactions = "R_user")

  expect_equal(selected, c("T_lac", "EX_glc", "R_user"))
})

test_that("rc_estimate_selected_demand_qp reports required workload plan", {
  plan <- rc_estimate_selected_demand_qp(
    n_pools = 10,
    selected_reactions = c("R1", "R2", "R2"),
    seconds_per_qp = 0.5,
    workers = 2,
    checkpoint_every = 7
  )

  expect_equal(plan$n_selected_reactions, 2L)
  expect_equal(plan$estimated_QP_count, 30L)
  expect_equal(plan$estimated_seconds_serial, 15)
  expect_equal(plan$estimated_seconds_parallel, 7.5)
  expect_equal(plan$expected_checkpoints, 5)
})
