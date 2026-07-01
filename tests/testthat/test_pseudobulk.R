test_that("rc_pool_mean computes sparse pool means", {
  skip_if_not_installed("Matrix")

  mat <- Matrix::Matrix(
    c(
      1, 0, 3, 5,
      0, 2, 4, 6,
      1, 1, 1, 1
    ),
    nrow = 3,
    ncol = 4,
    sparse = TRUE
  )
  rownames(mat) <- paste0("gene", seq_len(3))
  colnames(mat) <- paste0("cell", seq_len(4))
  pool_map <- data.frame(
    pool_id = c("pool_1", "pool_1", "pool_2", "pool_2"),
    cell_id = colnames(mat),
    stringsAsFactors = FALSE
  )

  out <- rc_pool_mean(mat, pool_map)
  expected <- cbind(
    pool_1 = Matrix::rowMeans(mat[, c("cell1", "cell2"), drop = FALSE]),
    pool_2 = Matrix::rowMeans(mat[, c("cell3", "cell4"), drop = FALSE])
  )
  expect_equal(out, expected)
  expect_identical(rownames(out), rownames(mat))
})

test_that("rc_pool_detection computes detection fractions from counts", {
  skip_if_not_installed("Matrix")

  counts <- Matrix::Matrix(
    c(
      1, 0, 3, 0,
      0, 0, 4, 6,
      1, 1, 0, 0
    ),
    nrow = 3,
    ncol = 4,
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(3))
  colnames(counts) <- paste0("cell", seq_len(4))
  pool_map <- data.frame(
    pool_id = c("pool_1", "pool_1", "pool_2", "pool_2"),
    cell_id = colnames(counts),
    stringsAsFactors = FALSE
  )

  out <- rc_pool_detection(counts, pool_map)
  expected <- cbind(
    pool_1 = Matrix::rowMeans(counts[, c("cell1", "cell2"), drop = FALSE] > 0),
    pool_2 = Matrix::rowMeans(counts[, c("cell3", "cell4"), drop = FALSE] > 0)
  )
  expect_equal(out, expected)
  expect_true(all(out >= 0 & out <= 1))
})

test_that("pool-level pseudobulk validates pool map cell IDs", {
  skip_if_not_installed("Matrix")

  mat <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
  rownames(mat) <- "gene1"
  colnames(mat) <- "cell1"
  pool_map <- data.frame(pool_id = "pool_1", cell_id = "missing_cell")

  expect_error(rc_pool_mean(mat, pool_map), "absent from matrix columns")
  expect_error(rc_pool_detection(mat, pool_map), "absent from matrix columns")
})

test_that("raw-count pseudobulk sum then logCPM works and empty pools are rejected", {
  skip_if_not_installed("Matrix")
  counts <- Matrix::Matrix(c(1, 0, 3, 0, 0, 0), nrow = 2, ncol = 3, sparse = TRUE)
  rownames(counts) <- c("g1", "g2"); colnames(counts) <- paste0("c", 1:3)
  pool_map <- data.frame(pool_id = c("p1", "p1", "p2"), cell_id = colnames(counts), skipped = FALSE)
  pb <- rc_pseudobulk_counts(counts, pool_map)
  expect_equal(as.numeric(pb[, "p1"]), c(1, 3))
  expect_error(rc_logcpm(pb), "Empty pools")
  filtered <- rc_filter_empty_pools(pb, data.frame(pool_id = c("p1", "p2")))
  expect_equal(colnames(filtered$counts), "p1")
  expect_true(all(is.finite(rc_logcpm(filtered$counts))))
})
