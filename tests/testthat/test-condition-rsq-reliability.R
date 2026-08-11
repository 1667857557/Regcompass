test_that("penalty q is one for valid active targets and ignores R-squared", {
  grn_result <- list(
    tf_peak_gene_condition = data.frame(
      target = c("G1", "G1", "G2", "G3"),
      condition = c("A", "B", "A", "A"),
      cell_type = c("T", "T", "T", "B"),
      rsq = c(0.01, 0.99, NA, -4),
      rsq_oof = c(-10, NA, 0.5, 2),
      stringsAsFactors = FALSE
    )
  )
  unit_meta <- data.frame(
    unit_id = c("u1", "u2", "u3", "u4"),
    condition = c("A", "B", "A", "A"),
    cell_type = c("T", "T", "B", "T"),
    stringsAsFactors = FALSE
  )
  template <- matrix(
    0, 4, 4,
    dimnames = list(c("g1", "g2", "g3", "g4"), unit_meta$unit_id)
  )

  q <- RegCompassR:::.rc_active_target_penalty_q(
    grn_result = grn_result,
    unit_meta = unit_meta,
    condition_col = "condition",
    celltype_col = "cell_type",
    template = template
  )

  expect_equal(q["g1", "u1"], 1)
  expect_equal(q["g1", "u2"], 1)
  expect_equal(q["g2", "u1"], 1)
  expect_equal(q["g3", "u3"], 1)
  expect_true(is.na(q["g1", "u3"]))
  expect_true(is.na(q["g2", "u2"]))
  expect_true(is.na(q["g4", "u1"]))
})

test_that("penalty q is determined only by the final active-edge table", {
  active <- list(
    tf_peak_gene_condition = data.frame(
      target = "G1", condition = "A", cell_type = "T",
      rsq = NA_real_, rsq_oof = NA_real_, stringsAsFactors = FALSE
    )
  )
  inactive <- list(tf_peak_gene_condition = data.frame())
  unit_meta <- data.frame(
    unit_id = "u1", condition = "A", cell_type = "T",
    stringsAsFactors = FALSE
  )
  template <- matrix(0, 1, 1, dimnames = list("g1", "u1"))

  q_active <- RegCompassR:::.rc_active_target_penalty_q(
    active, unit_meta, "condition", "cell_type", template
  )
  q_inactive <- RegCompassR:::.rc_active_target_penalty_q(
    inactive, unit_meta, "condition", "cell_type", template
  )

  expect_equal(q_active[[1L]], 1)
  expect_true(is.na(q_inactive[[1L]]))
})

test_that("combined Pando projection does not reuse route-level R-squared weights", {
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_project_pando_by_celltype)),
    collapse = "\n"
  )
  expect_match(body_text, ".rc_active_target_penalty_q", fixed = TRUE)
  expect_false(grepl("part$reliability", body_text, fixed = TRUE))
  expect_false(grepl("standard$reliability", body_text, fixed = TRUE))
  expect_match(body_text, "R2 and OOF R2 remain fit", fixed = TRUE)
})
