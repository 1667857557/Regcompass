test_that("Pando exposes the shared GRN design contract", {
  skip_if_not_installed("Pando")
  expect_true(utils::packageVersion("Pando") >= "1.1.2")
  exports <- getNamespaceExports("Pando")
  expect_true("prepare_grn_design" %in% exports)
  expect_true("validate_grn_design" %in% exports)
  expect_true(exists(
    "prepare_grn_design.GRNData",
    envir = asNamespace("Pando"),
    inherits = FALSE
  ))
})

test_that("Pando candidate builder preserves motif peak target provenance", {
  skip_if_not_installed("Pando")
  builder <- get(
    ".pando_candidate_edge_table",
    envir = asNamespace("Pando"),
    inherits = FALSE
  )
  peaks2gene <- Matrix::Matrix(
    matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE),
    sparse = TRUE,
    dimnames = list(c("G1", "G2"), c("R1", "R2"))
  )
  peaks2motif <- Matrix::Matrix(
    matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE),
    sparse = TRUE,
    dimnames = list(c("R1", "R2"), c("M1", "M2"))
  )
  motif2tf <- Matrix::Matrix(
    matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE),
    sparse = TRUE,
    dimnames = list(c("M1", "M2"), c("TF1", "TF2"))
  )
  edges <- builder(
    peaks2gene = peaks2gene,
    peaks2motif = peaks2motif,
    motif2tf = motif2tf,
    region_to_peak = c(R1 = "chr1-1-10", R2 = "chr1-20-30")
  )
  expect_identical(
    edges$edge_id,
    c("TF1::R1::G1", "TF2::R1::G1", "TF2::R2::G2")
  )
  expect_identical(
    edges$atac_feature_id,
    c("chr1-1-10", "chr1-1-10", "chr1-20-30")
  )
  expect_false(anyDuplicated(edges$edge_id))
})
