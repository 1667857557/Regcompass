test_that("penalty q distinguishes eligible, evaluated-neutral, and unavailable targets", {
  grn_result <- list(
    tf_peak_gene_condition_all = data.frame(
      target = c("G1", "G2", "G3", "G4"),
      condition = c("A", "A", "A", "B"),
      cell_type = c("T", "T", "T", "T"),
      target_model_evaluated = c(TRUE, TRUE, FALSE, TRUE),
      penalty_eligible = c(TRUE, FALSE, FALSE, TRUE),
      stringsAsFactors = FALSE
    )
  )
  unit_meta <- data.frame(
    unit_id = c("uA", "uB"), condition = c("A", "B"),
    cell_type = c("T", "T"), stringsAsFactors = FALSE
  )
  template <- matrix(
    0, 5, 2,
    dimnames = list(c("g1", "g2", "g3", "g4", "g5"), unit_meta$unit_id)
  )
  q <- RegCompassR:::.rc_active_target_penalty_q(
    grn_result, unit_meta, "condition", "cell_type", template
  )
  expect_equal(q["g1", "uA"], 1)
  expect_equal(q["g2", "uA"], 0)
  expect_true(is.na(q["g3", "uA"]))
  expect_equal(q["g4", "uB"], 1)
  expect_true(is.na(q["g5", "uA"]))
  expect_true(is.na(q["g1", "uB"]))
})

test_that("condition target gate uses final full-data R2 and leaves Pando significance intact", {
  fit <- list(
    padj_threshold = 0.05,
    coefficients = data.frame(
      edge_id = c("G1||TF||P", "G2||TF||P"),
      target = c("G1", "G2"), condition = c("A", "A"),
      estimate = c(0.4, 0.5), estimable = TRUE,
      padj = c(0.01, 0.01), statistically_supported = TRUE,
      global_support = TRUE, local_support = FALSE,
      active = TRUE, significant = TRUE,
      penalty_effect = c(0.4, 0.5), stringsAsFactors = FALSE
    ),
    fit = data.frame(
      target = c("G1", "G2"), condition = c("A", "A"),
      fit_status = c("ok", "ok"),
      rsq = c(0.2, 0.01), rsq_in_sample = c(0.2, 0.01),
      rsq_oof = c(-0.3, 0.8), stringsAsFactors = FALSE
    )
  )
  diagnostics <- RegCompassR:::.rc_condition_fit_diagnostics_for_coefficients(
    fit, fit$coefficients, target_rsq_threshold = 0.05
  )
  expect_identical(diagnostics$target_model_supported, c(TRUE, FALSE))
  annotated <- cbind(fit$coefficients, diagnostics)
  annotated$padj_threshold <- 0.05
  gate <- RegCompassR:::.rc_condition_penalty_gate(
    annotated, padj_threshold = 0.05, target_rsq_threshold = 0.05
  )
  expect_identical(gate, c(TRUE, FALSE))
  expect_true(all(fit$coefficients$active))
  expect_true(all(fit$coefficients$significant))
})

test_that("q zero is neutral even when projection is unavailable", {
  projection <- matrix(c(NA_real_, 2), 1, 2,
                       dimnames = list("g", c("u0", "u1")))
  q <- matrix(c(0, 1), 1, 2, dimnames = dimnames(projection))
  scale <- matrix(1, 1, 2, dimnames = dimnames(projection))
  modifier <- RegCompassR:::.rc_scaled_regulatory_modifier(projection, q, scale)
  expect_equal(modifier[[1L]], 0)
  expect_true(is.finite(modifier[[2L]]))
})

test_that("combined Pando projection uses tri-state target eligibility, not OOF weights", {
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_project_pando_by_celltype)), collapse = "\n"
  )
  expect_match(body_text, ".rc_active_target_penalty_q", fixed = TRUE)
  expect_false(grepl("part$reliability", body_text, fixed = TRUE))
  expect_false(grepl("standard$reliability", body_text, fixed = TRUE))
  expect_match(body_text, "q=1 eligible; q=0 evaluated neutral", fixed = TRUE)
})
