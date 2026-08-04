test_that("condition and standard routes form a strict partition", {
  expect_true(.rc_validate_pando_route_partition(
    "epithelial_like", "stem_like"
  ))
  expect_error(
    .rc_validate_pando_route_partition(
      "epithelial_like", "epithelial_like"
    ),
    "both condition and standard"
  )
  expect_error(
    .rc_validate_pando_route_partition(
      c("epithelial_like", "epithelial_like"), character()
    ),
    "complete and unique"
  )
})

test_that("paired cells cannot be emitted by multiple Pando routes", {
  results <- list(
    list(paired_cell_metadata = data.frame(
      cell_id = c("cell1", "cell2"), stringsAsFactors = FALSE
    )),
    list(paired_cell_metadata = data.frame(
      cell_id = c("cell2", "cell3"), stringsAsFactors = FALSE
    ))
  )
  expect_error(
    .rc_validate_pando_result_cell_partition(results),
    "more than one Pando route"
  )
})

test_that("projection overlays reject all route overlap", {
  target <- matrix(
    c(1, NA_real_), nrow = 1,
    dimnames = list("gene", c("unit1", "unit2"))
  )
  incoming <- matrix(
    c(1, 2), nrow = 1,
    dimnames = dimnames(target)
  )
  expect_error(
    .rc_overlay_projection(target, incoming, "Standard-Pando"),
    "overlaps an existing route"
  )

  incoming[1, 1] <- NA_real_
  expect_identical(
    .rc_overlay_projection(target, incoming, "Standard-Pando"),
    matrix(c(1, 2), nrow = 1, dimnames = dimnames(target))
  )
})
