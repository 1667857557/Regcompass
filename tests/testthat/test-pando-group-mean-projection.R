test_that("Pando projection uses paired TF-ATAC interaction with 25% shrinkage", {
  rna <- matrix(c(1, 3), ncol = 1L,
                dimnames = list(c("c1", "c2"), "TF1"))
  atac <- matrix(c(4, 2), ncol = 1L,
                 dimnames = list(c("c1", "c2"), "chr1-10-20"))
  edge <- data.frame(
    tf = "TF1", target = "G", region = "chr1-10-20", estimate = 2,
    stringsAsFactors = FALSE
  )

  score <- RegCompassR:::.rc_pando_projection_from_group_means(
    rna, atac, list(u1 = edge), list(u1 = c("c1", "c2")), "G"
  )

  paired <- mean(c(1, 3) * c(4, 2))
  product_of_means <- mean(c(1, 3)) * mean(c(4, 2))
  expected <- 2 * (0.75 * paired + 0.25 * product_of_means)
  expect_equal(score["u1", "g"], expected)
})

test_that("Pando paired-shrinkage projection supports sparse Matrix inputs", {
  rna <- Matrix::Matrix(
    matrix(c(1, 3), ncol = 1L,
           dimnames = list(c("c1", "c2"), "TF1")),
    sparse = TRUE
  )
  atac <- Matrix::Matrix(
    matrix(c(4, 2), ncol = 1L,
           dimnames = list(c("c1", "c2"), "chr1-10-20")),
    sparse = TRUE
  )
  edge <- data.frame(
    tf = "TF1", target = "G", region = "chr1-10-20", estimate = 2,
    stringsAsFactors = FALSE
  )

  score <- RegCompassR:::.rc_pando_projection_from_group_means(
    rna, atac, list(u1 = edge), list(u1 = c("c1", "c2")), "G"
  )

  paired <- mean(c(1, 3) * c(4, 2))
  product_of_means <- mean(c(1, 3)) * mean(c(4, 2))
  expect_equal(
    score["u1", "g"],
    2 * (0.75 * paired + 0.25 * product_of_means)
  )
})

test_that("Pando projection weight has exact paired and legacy boundaries", {
  rna <- matrix(c(1, 3), ncol = 1L,
                dimnames = list(c("c1", "c2"), "TF1"))
  atac <- matrix(c(4, 2), ncol = 1L,
                 dimnames = list(c("c1", "c2"), "chr1-10-20"))
  edge <- data.frame(
    tf = "TF1", target = "G", region = "chr1-10-20", estimate = 2,
    stringsAsFactors = FALSE
  )
  args <- list(
    rna = rna,
    atac = atac,
    edges_by_group = list(u1 = edge),
    cells_by_group = list(u1 = c("c1", "c2")),
    targets = "G"
  )

  paired <- do.call(
    RegCompassR:::.rc_pando_projection_from_group_means,
    c(args, list(product_of_means_weight = 0))
  )
  legacy <- do.call(
    RegCompassR:::.rc_pando_projection_from_group_means,
    c(args, list(product_of_means_weight = 1))
  )

  expect_equal(paired["u1", "g"], mean(2 * c(1, 3) * c(4, 2)))
  expect_equal(
    legacy["u1", "g"],
    2 * mean(c(1, 3)) * mean(c(4, 2))
  )
})

test_that("Pando projection retains paired covariance when marginal means match", {
  rna <- matrix(
    c(0, 2, 0, 2), ncol = 1L,
    dimnames = list(c("a1", "a2", "b1", "b2"), "TF1")
  )
  atac <- matrix(
    c(0, 2, 2, 0), ncol = 1L,
    dimnames = list(c("a1", "a2", "b1", "b2"), "chr1-10-20")
  )
  edge <- data.frame(
    tf = "TF1", target = "G", region = "chr1-10-20", estimate = 1,
    stringsAsFactors = FALSE
  )

  score <- RegCompassR:::.rc_pando_projection_from_group_means(
    rna, atac,
    edges_by_group = list(u_pos = edge, u_discordant = edge),
    cells_by_group = list(
      u_pos = c("a1", "a2"),
      u_discordant = c("b1", "b2")
    ),
    targets = "G"
  )

  expect_equal(score["u_pos", "g"], 1.75)
  expect_equal(score["u_discordant", "g"], 0.25)
  expect_gt(score["u_pos", "g"], score["u_discordant", "g"])
})

test_that("Pando paired-shrinkage projection sums edges by target", {
  rna <- matrix(c(1, 3, 2, 4), nrow = 2L,
                dimnames = list(c("c1", "c2"), c("TF1", "TF2")))
  atac <- matrix(c(4, 2, 5, 1), nrow = 2L,
                 dimnames = list(
                   c("c1", "c2"), c("chr1-10-20", "chr1-30-40")
                 ))
  edge <- data.frame(
    tf = c("TF1", "TF2"), target = "G",
    region = c("chr1-10-20", "chr1-30-40"), estimate = c(2, -1),
    stringsAsFactors = FALSE
  )

  score <- RegCompassR:::.rc_pando_projection_from_group_means(
    rna, atac, list(u1 = edge), list(u1 = c("c1", "c2")), "G"
  )
  interaction1 <-
    0.75 * mean(rna[, "TF1"] * atac[, "chr1-10-20"]) +
    0.25 * mean(rna[, "TF1"]) * mean(atac[, "chr1-10-20"])
  interaction2 <-
    0.75 * mean(rna[, "TF2"] * atac[, "chr1-30-40"]) +
    0.25 * mean(rna[, "TF2"]) * mean(atac[, "chr1-30-40"])
  expected <- 2 * interaction1 - interaction2
  expect_equal(score["u1", "g"], expected)
})

test_that("Pando group-mean projection rejects broken stage hand-offs", {
  rna <- matrix(1, nrow = 1L, dimnames = list("c1", "TF1"))
  atac <- matrix(1, nrow = 1L, dimnames = list("c1", "chr1-10-20"))
  edge <- data.frame(
    tf = "TF1", target = "G", region = "chr1-10-20", estimate = 1,
    stringsAsFactors = FALSE
  )

  expect_error(
    RegCompassR:::.rc_pando_projection_from_group_means(
      rna, atac, list(u2 = edge), list(u1 = "c1"), "G"
    ),
    "same unique names"
  )
  expect_error(
    RegCompassR:::.rc_pando_projection_from_group_means(
      rna, atac, list(u1 = edge), list(u1 = c("c1", "c2")), "G"
    ),
    "missing RNA=1, missing ATAC=1"
  )
  expect_error(
    RegCompassR:::.rc_pando_projection_from_group_means(
      rna, atac, list(u1 = edge), list(u1 = "c1"), c("G", "g")
    ),
    "targets must be unique"
  )
  expect_error(
    RegCompassR:::.rc_pando_projection_from_group_means(
      rna, atac, list(u1 = edge), list(u1 = "c1"), "G",
      product_of_means_weight = -0.1
    ),
    "in \\[0, 1\\]"
  )
  expect_error(
    RegCompassR:::.rc_pando_projection_from_group_means(
      rna, atac, list(u1 = edge), list(u1 = "c1"), "G",
      product_of_means_weight = 1.1
    ),
    "in \\[0, 1\\]"
  )
})

test_that("Pando group names align by ID rather than list position", {
  rna <- matrix(c(1, 2), nrow = 2L,
                dimnames = list(c("c1", "c2"), "TF1"))
  atac <- matrix(c(3, 4), nrow = 2L,
                 dimnames = list(c("c1", "c2"), "chr1-10-20"))
  edge1 <- data.frame(
    tf = "TF1", target = "G", region = "chr1-10-20", estimate = 1
  )
  edge2 <- edge1
  edge2$estimate <- 2

  score <- RegCompassR:::.rc_pando_projection_from_group_means(
    rna, atac,
    edges_by_group = list(u2 = edge2, u1 = edge1),
    cells_by_group = list(u1 = "c1", u2 = "c2"),
    targets = "G"
  )

  expect_equal(score["u1", "g"], 3)
  expect_equal(score["u2", "g"], 16)
})