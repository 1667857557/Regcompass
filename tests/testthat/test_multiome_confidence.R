test_that("concordance null correction uses finite-n discrete baseline", {
  p1 <- matrix(c(0, 1, 0.5, 0.5), nrow = 1, dimnames = list("g1", paste0("p", 1:4)))
  p2 <- matrix(c(0, 1, 0.5, 0.5), nrow = 1, dimnames = list("g1", paste0("p", 1:4)))
  out <- rc_concordance_null_correct(p1, p2)
  expect_equal(as.numeric(out), rep(1, 4))
})

test_that("Fisher shrinkage separates positive support from negative discordance", {
  x <- rbind(g1 = 1:8, g2 = 1:8)
  y <- rbind(g1 = 1:8, g2 = 8:1)
  out <- rc_fisher_shrink(x, y, n0 = 30)
  expect_gt(out["g1", "rho_shrink"], 0)
  expect_lt(out["g2", "rho_shrink"], 0)
  expect_true(all(abs(out$rho_shrink) < 1))
})

test_that("gene confidence remains nonnegative and bounded", {
  mat <- matrix(c(0, 0.5, 1, 1), nrow = 2, dimnames = list(c("g1", "g2"), c("p1", "p2")))
  conf <- rc_gene_confidence(mat, rel_ra_pos = c(1, 0), det_rna = mat, qc = c(1, 0.5))
  expect_true(all(conf >= 0 & conf <= 1))
})

test_that("concordance correction is zero-power safe for n equals one", {
  p1 <- matrix(0.5, nrow = 1, dimnames = list("g1", "p1"))
  p2 <- matrix(0.5, nrow = 1, dimnames = list("g1", "p1"))
  expect_equal(as.numeric(rc_concordance_null_correct(p1, p2)), 1)
})

test_that("concordance correction can use stratum-specific n", {
  p1 <- matrix(c(0, 1, 0.5, 0.5), nrow = 1, dimnames = list("g1", paste0("p", 1:4)))
  p2 <- p1
  pool_meta <- data.frame(pool_id = colnames(p1), cell_type = c("A", "A", "B", "B"))
  out <- rc_concordance_null_correct(p1, p2, pool_meta = pool_meta, stratum_col = "cell_type")
  expect_equal(as.numeric(out), rep(1, 4))
})

test_that("gene confidence aligns named reliability vectors by gene", {
  mat <- matrix(c(0.2, 0.8, 0.4, 0.6), nrow = 2, dimnames = list(c("g1", "g2"), c("p1", "p2")))
  rel <- c(g2 = 0, g1 = 1)
  conf <- rc_gene_confidence(mat, rel_ra_pos = rel, det_rna = mat, qc = c(1, 1))
  expect_gt(conf["g1", "p1"], conf["g2", "p1"])
})
