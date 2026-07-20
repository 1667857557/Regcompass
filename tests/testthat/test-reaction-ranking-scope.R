test_that("reaction priority ranks are stratified by medium", {
  row_ids <- c(
    "reaction=R1::direction=forward::medium=base",
    "reaction=R2::direction=forward::medium=base",
    "reaction=R1::direction=forward::medium=stress",
    "reaction=R2::direction=forward::medium=stress"
  )
  microcompass <- list(
    penalty = matrix(
      c(0.2, 0.4, 5, 10),
      ncol = 1,
      dimnames = list(row_ids, "u1")
    ),
    vmax = matrix(
      1,
      nrow = 4,
      ncol = 1,
      dimnames = list(row_ids, "u1")
    ),
    unit_meta = data.frame(
      unit_id = "u1",
      condition = "A",
      cell_type = "T",
      stringsAsFactors = FALSE
    ),
    params = list(omega = 1)
  )

  answer <- .rc_condition_penalty_comparison(microcompass)
  expect_identical(answer$ranking_scope, "condition_x_celltype_x_medium")
  expect_equal(
    answer$ranking$priority_rank[
      answer$ranking$medium_scenario == "base"
    ],
    c(1L, 2L)
  )
  expect_equal(
    answer$ranking$priority_rank[
      answer$ranking$medium_scenario == "stress"
    ],
    c(1L, 2L)
  )
  expect_true(all(
    answer$ranking$ranking_scope == "condition_x_celltype_x_medium"
  ))
})
