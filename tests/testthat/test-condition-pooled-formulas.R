test_that("regulatory integration is bounded and zero preserving", {
  C <- matrix(c(0, 0.20, 0.20, 0.80), nrow = 4,
    dimnames = list(paste0("g", 1:4), "mc1"))
  R <- matrix(c(1, 1, -1, 1), nrow = 4, dimnames = dimnames(C))
  out <- .rc_integrate_regulatory_support(C, R, alpha = 1)
  expect_equal(out["g1", "mc1"], 0)
  expect_gt(out["g2", "mc1"], C["g2", "mc1"])
  expect_lt(out["g3", "mc1"], C["g3", "mc1"])
  expect_true(all(out >= 0 & out <= 1))
})

test_that("reaction penalty is positive and decreases with expression", {
  E <- matrix(c(0, 1, 3, NA_real_), nrow = 4,
    dimnames = list(paste0("R", 1:4), "mc1"))
  answer <- rc_compute_multiome_penalty(E)
  P <- answer$penalty[, "mc1"]
  expect_equal(P[["R1"]], 1)
  expect_gt(P[["R1"]], P[["R2"]])
  expect_gt(P[["R2"]], P[["R3"]])
  expect_equal(P[["R4"]], 1)
  expect_true(all(is.finite(P) & P > 0))
})

test_that("condition Pando remains independent by broad cell type", {
  implementation <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  expect_match(implementation, "Pando::infer_condition_grn", fixed = TRUE)
  expect_match(implementation, "cell_type = cell_type", fixed = TRUE)
  expect_match(implementation, "exact-edge union", fixed = TRUE)
  expect_false("sample_col" %in%
                 names(formals(.rc_fit_condition_grns_by_cell_type)))
  expect_false("strict_biological_defaults" %in% names(formals(rc_run_regcompass)))
})

test_that("condition Layer 1 delegates paired-cell projection to Pando", {
  helper <- paste(
    deparse(body(.rc_condition_pando_projection)), collapse = "\n"
  )
  layer1 <- paste(
    deparse(body(.rc_cell_first_projection_layer1)), collapse = "\n"
  )
  expect_match(helper, "Pando::project_condition_grn_cells", fixed = TRUE)
  expect_match(helper, "significant_only = TRUE", fixed = TRUE)
  expect_match(helper, "Pando::aggregate_condition_grn_projection", fixed = TRUE)
  expect_match(helper, "common = primary", fixed = TRUE)
  expect_match(helper, "primary * 0", fixed = TRUE)
  expect_match(helper, "sqrt(pmin(1, pmax(0", fixed = TRUE)
  expect_false(grepl("project_condition_grn_primary_cells", helper, fixed = TRUE))
  expect_match(layer1, ".rc_scaled_oof_modifier", fixed = TRUE)
  expect_match(layer1, "gene_projection_condition_full_oof", fixed = TRUE)
})

test_that("single-condition scoring uses penalty per required target flux", {
  row_ids <- c(
    "reaction=R1::direction=forward::medium=base",
    "reaction=R2::direction=forward::medium=base"
  )
  units <- c("u1", "u2")
  microcompass <- list(
    penalty = matrix(c(0.2, 0.4, 0.3, 0.5), nrow = 2,
      dimnames = list(row_ids, units)),
    vmax = matrix(c(10, 1, 10, 1), nrow = 2,
      dimnames = list(row_ids, units)),
    unit_meta = data.frame(unit_id = units, condition = "A", cell_type = "T"),
    params = list(omega = 0.95)
  )
  answer <- .rc_condition_penalty_comparison(microcompass)
  expect_identical(answer$analysis_mode, "single_condition_reaction_ranking")
  expect_equal(nrow(answer$contrast), 0L)
  expect_true(all(is.finite(answer$ranking$median_penalty_per_target_flux)))
})
