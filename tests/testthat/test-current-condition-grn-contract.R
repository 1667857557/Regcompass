test_that("canonical Pando bridge is independent by cell type", {
  defaults <- eval(
    formals(.rc_fit_condition_grns_by_cell_type)$pando_infer_args
  )
  expect_identical(defaults$candidate_screen, "motif_domain")
  expect_identical(defaults$peak_cor, 0)
  expect_identical(defaults$condition_mix, 0.5)
  expect_identical(defaults$condition_weight, "equal")
  expect_true(defaults$scale)
  expect_false(any(c(
    "method", "cv_block_col", "sample_col"
  ) %in% names(defaults)))

  implementation <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)),
    collapse = "\n"
  )
  expect_match(implementation, "Pando::infer_condition_grn", fixed = TRUE)
  expect_match(implementation, "cell_type = cell_type", fixed = TRUE)
  expect_match(implementation, ".rc_extract_condition_grn_contract",
               fixed = TRUE)
})

test_that("contract extraction requires within-cell-type cell OOF", {
  implementation <- paste(
    deparse(body(.rc_extract_condition_grn_contract)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    "condition_sparse_within_cell_type_oof_refit",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "within_cell_type_condition_stratified_cells",
    fixed = TRUE
  )
  expect_match(implementation, "fit_cell_type", fixed = TRUE)
  expect_false(grepl("sample_blocked", implementation, fixed = TRUE))
  expect_false(grepl("cell_type_blocked", implementation, fixed = TRUE))
})

test_that("Layer 1 uses cell-first Pando projection and bounded OOF reliability", {
  implementation <- paste(
    deparse(body(.rc_cell_first_projection_layer1)),
    collapse = "\n"
  )
  expect_match(implementation, "Pando::project_condition_grn_cells",
               fixed = TRUE)
  expect_match(implementation, "Pando::aggregate_condition_grn_projection",
               fixed = TRUE)
  expect_match(
    implementation,
    "as.character(unit_meta[[celltype_col]]) == fit$cell_type",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "sqrt(pmin(1, pmax(0",
    fixed = TRUE
  )
})

test_that("metacell builder name states its condition-by-cell-type scope", {
  expect_true(exists(
    ".rc_make_condition_celltype_metacells", inherits = TRUE
  ))
  expect_false(exists(
    ".rc_make_condition_pooled_metacells", inherits = TRUE
  ))
})

test_that("regulatory integration is bounded and zero preserving", {
  rna <- matrix(
    c(0, 0.5, 1), nrow = 1,
    dimnames = list("g1", c("u1", "u2", "u3"))
  )
  modifier <- matrix(
    c(1, 1, -1), nrow = 1, dimnames = dimnames(rna)
  )
  observed <- .rc_integrate_regulatory_support(rna, modifier, alpha = 1)

  expect_equal(observed[1, 1], 0)
  expect_equal(observed[1, 3], 1)
  expect_true(observed[1, 2] > 0.5)
  expect_true(all(observed >= 0 & observed <= 1))
})
