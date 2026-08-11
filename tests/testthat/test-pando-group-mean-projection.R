test_that("Pando projection multiplies separate metacell modality means", {
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

  expect_equal(score["u1", "g"], 2 * mean(c(1, 3)) * mean(c(4, 2)))
  expect_false(isTRUE(all.equal(
    score["u1", "g"], mean(2 * c(1, 3) * c(4, 2))
  )))
})

test_that("Pando group-mean projection supports sparse Matrix inputs", {
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

  expect_equal(score["u1", "g"], 2 * mean(c(1, 3)) * mean(c(4, 2)))
})

test_that("Pando group-mean projection sums edges by target", {
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
  expected <- 2 * mean(rna[, "TF1"]) * mean(atac[, "chr1-10-20"]) -
    mean(rna[, "TF2"]) * mean(atac[, "chr1-30-40"])
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