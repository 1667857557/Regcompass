test_that("condition Scheme E penalty q follows fitted continuous availability, not R2", {
  all_edges <- data.frame(
    target = c("G1", "G1", "G2", "G3", "G4"),
    condition = c("A", "B", "A", "A", "A"),
    cell_type = c("T", "T", "T", "B", "T"),
    analysis_mode = "condition_grn",
    fit_status = c("ok", "ok", "ok", "ok", "failed"),
    rsq = c(0.8, 0.01, NA_real_, 0.4, 0.9),
    penalty_effect = c(0.5, 0.4, -0.2, 0.3, 0.1),
    stringsAsFactors = FALSE
  )
  active_edges <- all_edges[all_edges$fit_status == "ok", , drop = FALSE]
  grn_result <- list(
    tf_peak_gene_condition_all = all_edges,
    tf_peak_gene_condition = active_edges
  )
  unit_meta <- data.frame(
    unit_id = c("u1", "u2", "u3", "u4"),
    condition = c("A", "B", "A", "C"),
    cell_type = c("T", "T", "B", "T"),
    stringsAsFactors = FALSE
  )
  template <- matrix(
    0, 5, 4,
    dimnames = list(c("g1", "g2", "g3", "g4", "g5"), unit_meta$unit_id)
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
  expect_true(is.na(q["g4", "u1"]))
  expect_true(is.na(q["g5", "u1"]))
  expect_true(is.na(q["g1", "u4"]))
})

test_that("q zero remains a neutral modifier for legacy rejected routes", {
  projection <- matrix(
    c(NA_real_, 2), nrow = 1,
    dimnames = list("g1", c("u_rejected", "u_active"))
  )
  reliability <- matrix(
    c(0, 1), nrow = 1, dimnames = dimnames(projection)
  )
  scale <- matrix(1, nrow = 1, ncol = 2, dimnames = dimnames(projection))

  modifier <- RegCompassR:::.rc_scaled_regulatory_modifier(
    projection, reliability, scale
  )
  expect_equal(modifier["g1", "u_rejected"], 0)
  expect_true(is.finite(modifier["g1", "u_active"]))
})

test_that("combined projection documents diagnostic-only R2 for condition Scheme E", {
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_project_pando_by_celltype)),
    collapse = "\n"
  )
  q_text <- paste(
    deparse(body(RegCompassR:::.rc_active_target_penalty_q)),
    collapse = "\n"
  )
  evaluated_text <- paste(
    deparse(body(RegCompassR:::.rc_penalty_evaluated_rows)),
    collapse = "\n"
  )
  expect_match(body_text, ".rc_active_target_penalty_q", fixed = TRUE)
  expect_match(q_text, ".rc_penalty_evaluated_rows", fixed = TRUE)
  expect_match(evaluated_text, "condition_evaluated <- fit_ok", fixed = TRUE)
  expect_match(
    body_text,
    "BH/R2 do not gate",
    fixed = TRUE
  )
})

test_that("standard Pando projection keeps its independent filtering route", {
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_standard_pando_projection)),
    collapse = "\n"
  )
  expect_false(grepl("Scheme-E", body_text, fixed = TRUE))
  expect_false(grepl("scheme_e", body_text, fixed = TRUE))
})
