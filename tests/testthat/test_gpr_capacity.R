test_that("rc_reaction_capacity returns reaction by pool matrix", {
  gpr_table <- data.frame(
    reaction_id = c("R_HEX", "R_PFK"),
    gpr = c("HK1 or HK2", "PFKM and PFKL"),
    stringsAsFactors = FALSE
  )
  parsed <- rc_parse_gpr_table(gpr_table)
  gene_score <- matrix(
    c(0.8, 0.4, 0.7, 0.5, 0.2, 0.9, 0.6, 0.4),
    nrow = 4,
    dimnames = list(c("hk1", "hk2", "pfkm", "pfkl"), c("pool1", "pool2"))
  )

  out <- rc_reaction_capacity(parsed, gene_score, promiscuity_mode = "none", tau = 0.08)
  expect_equal(dim(out), c(2L, 2L))
  expect_identical(rownames(out), c("R_HEX", "R_PFK"))
  expect_identical(colnames(out), c("pool1", "pool2"))
  expect_equal(out["R_HEX", "pool1"], 0.8 + 0.4)
})

test_that("rc_run_layer1_capacity returns MVP v0.3 outputs", {
  gpr_table <- data.frame(
    reaction_id = c("R_HEX", "R_PFK"),
    gpr = c("HK1 or HK2", "PFKM and PFKL"),
    stringsAsFactors = FALSE
  )
  expr <- matrix(
    c(10, 2, 5, 5, 4, 8, 6, 2),
    nrow = 4,
    dimnames = list(c("HK1", "HK2", "PFKM", "PFKL"), c("pool1", "pool2"))
  )
  detect <- matrix(
    c(1, 0.5, 0.8, 0.9, 1, 0.2, 0.7, 0.6),
    nrow = 4,
    dimnames = list(c("HK1", "HK2", "PFKM", "PFKL"), c("pool1", "pool2"))
  )

  out <- rc_run_layer1_capacity(gpr_table, expr, pool_detection = detect, promiscuity_mode = "sqrt", min_direct = 10)
  expect_true(all(c("reaction_capacity_L1", "reaction_confidence", "q95_diagnostics") %in% names(out)))
  expect_equal(dim(out$reaction_capacity_L1), c(2L, 2L))
  expect_true("mean_gpr_detection_rate" %in% colnames(out$reaction_confidence))
})

test_that("safe scale uses sigma-consistent MAD/IQR and clips z-scores", {
  x <- c(0, 0, 1, 1)
  expect_equal(rc_safe_scale(x, min_scale = 0.05), max(stats::mad(x, constant = 1.4826), stats::IQR(x) / 1.349, 0.05))
  X <- rbind(g1 = c(0, 100), g2 = c(1, 1))
  z <- rc_gene_zscore(X, z_clip = 2)
  expect_true(all(z <= 2 & z >= -2))
})

test_that("reaction confidence aggregates gene confidence by GPR genes and pools", {
  gprs <- list(r1 = list(c("g1", "g2")), r2 = list("g3"))
  gene_conf <- matrix(c(0.2, 0.8, 0.9, 0.4, 0.6, 0.7), nrow = 3,
                      dimnames = list(c("g1", "g2", "g3"), c("p1", "p2")))
  out <- rc_reaction_confidence(gprs, gene_confidence = gene_conf)
  expect_true(all(c("reaction_id", "pool_id", "reaction_confidence") %in% colnames(out)))
  expect_equal(out$reaction_confidence[out$reaction_id == "r1" & out$pool_id == "p1"], stats::median(c(0.2, 0.8)))
})
