test_that("rc_sample_aggregate uses sample-celltype medians", {
  score_mat <- matrix(
    c(1, 3, 10, 20,
      2, 4, 30, 50),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("R1", "R2"), paste0("P", 1:4))
  )
  pool_meta <- data.frame(
    pool_id = paste0("P", 1:4),
    sample_id = c("S1", "S1", "S2", "S2"),
    cell_type = c("T", "T", "T", "B"),
    stringsAsFactors = FALSE
  )

  out <- rc_sample_aggregate(score_mat, pool_meta)

  expect_equal(colnames(out), c("S1.T", "S2.B", "S2.T"))
  expect_equal(out[, "S1.T"], c(R1 = 2, R2 = 3))
  expect_equal(out[, "S2.T"], c(R1 = 10, R2 = 30))
})

test_that("rc_lm_by_reaction fits sample-level linear models and BH q-values", {
  Y <- matrix(
    c(1, 2, 3, 4,
      4, 3, 2, 1),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("R_up", "R_down"), paste0("S", 1:4))
  )
  sample_meta <- data.frame(
    condition = factor(c("ctrl", "ctrl", "stim", "stim"), levels = c("ctrl", "stim")),
    batch = c("A", "B", "A", "B"),
    row.names = paste0("S", 1:4)
  )

  out <- rc_lm_by_reaction(Y, sample_meta, "~ condition")

  stim <- out[out$term == "conditionstim", ]
  expect_equal(stim$reaction_id, c("R_up", "R_down"))
  expect_equal(stim$estimate, c(2, -2))
  expect_true(all(c("std_error", "statistic", "p_value", "q_value", "n") %in% colnames(out)))
})

test_that("rc_rank_regulators prioritizes supported candidate regulators", {
  evidence <- data.frame(
    regulator_id = c("TF1", "TF2", "TF3"),
    reaction_id = c("R1", "R1", "R1"),
    direct_association = c(0.9, 0.4, 0.6),
    adjusted_association = c(0.8, 0.5, 0.7),
    motif_support = c(1, 0, 1),
    enhancer_support = c(1, 0, 0)
  )

  out <- rc_rank_regulators(evidence)

  expect_equal(out$regulator_id[[1]], "TF1")
  expect_equal(out$evidence_tier[[1]], "motif-and-enhancer-supported")
  expect_true(all(c("rra_p_value", "q_value") %in% colnames(out)))
})

test_that("rc_rank_regulators ranks candidates within each reaction", {
  evidence <- data.frame(
    regulator_id = c("TF1", "TF2", "TF1", "TF2"),
    reaction_id = c("R1", "R1", "R2", "R2"),
    direct_association = c(0.9, 0.1, 0.1, 0.9),
    adjusted_association = c(0.8, 0.2, 0.2, 0.8)
  )

  out <- rc_rank_regulators(evidence)
  top_by_reaction <- out[out$candidate_rank == 1, c("reaction_id", "regulator_id")]

  expect_equal(top_by_reaction$reaction_id, c("R1", "R2"))
  expect_equal(top_by_reaction$regulator_id, c("TF1", "TF2"))
})

test_that("rc_sample_aggregate rejects missing sample or cell-type labels", {
  score_mat <- matrix(
    c(1, 2, 3, 4),
    nrow = 2L,
    dimnames = list(c("R1", "R2"), c("P1", "P2"))
  )
  pool_meta <- data.frame(
    pool_id = c("P1", "P2"),
    sample_id = c("S1", NA),
    cell_type = c("T", "T")
  )

  expect_error(rc_sample_aggregate(score_mat, pool_meta), "must not contain missing")
})

test_that("rc_rank_regulators excludes numeric identifier columns from default evidence", {
  evidence <- data.frame(
    regulator_id = c(101, 102),
    reaction_id = c(1, 1),
    direct_association = c(0.1, 0.9)
  )

  out <- rc_rank_regulators(evidence)

  expect_false("rank_regulator_id" %in% colnames(out))
  expect_false("rank_reaction_id" %in% colnames(out))
  expect_equal(out$regulator_id[[1]], 102)
})
