test_that("strict stratum helpers reject non-unique and missing metadata definitions", {
  expect_identical(.rc_strict_stratum_cols("sample_id", "condition", "cell_type"), c("condition", "sample_id", "cell_type"))
  expect_error(.rc_strict_stratum_cols("sample_id", "sample_id", "cell_type"), "three valid")

  meta <- data.frame(condition = c("ctrl", "ctrl"), sample_id = c("s1", ""), cell_type = c("T", "T"), row.names = c("c1", "c2"), stringsAsFactors = FALSE)
  expect_error(.rc_add_stratum_id(meta, c("condition", "sample_id", "cell_type")), "missing or empty")
})

test_that("post-metacell filter uses usable metacell boundaries and subsets bundle", {
  ids9 <- paste0("mc9_", seq_len(9))
  ids10 <- paste0("mc10_", seq_len(10))
  ids11 <- paste0("mc11_", seq_len(11))
  ids <- c(ids9, ids10, ids11, "mc10_low")
  meta <- data.frame(
    metacell_id = ids,
    condition = "ctrl",
    sample_id = rep(c("s9", "s10", "s11", "s10"), c(9, 10, 11, 1)),
    cell_type = "T",
    n_cells = c(rep(20L, 30), 5L),
    stringsAsFactors = FALSE
  )
  mc <- list(
    metacell_meta = meta,
    rna_counts = Matrix::Matrix(matrix(1, nrow = 2, ncol = length(ids), dimnames = list(c("g1", "g2"), ids)), sparse = TRUE),
    atac_counts = Matrix::Matrix(matrix(1, nrow = 2, ncol = length(ids), dimnames = list(c("p1", "p2"), ids)), sparse = TRUE),
    membership = data.frame(cell_id = paste0("c", seq_along(ids)), metacell_id = ids, stringsAsFactors = FALSE),
    fragment_manifest = data.frame(fragment_file = "frag.tsv.gz", object_cell = ids, fragment_barcode = ids, stringsAsFactors = FALSE)
  )

  out <- rc_filter_post_metacell_strata(mc, min_metacells = 10L, min_metacell_size = 20L)

  expect_false(any(grepl("^mc9_", out$metacell_meta$metacell_id)))
  expect_false("mc10_low" %in% out$metacell_meta$metacell_id)
  expect_true(all(c(ids10, ids11) %in% out$metacell_meta$metacell_id))
  expect_identical(colnames(out$rna_counts), out$metacell_meta$metacell_id)
  expect_identical(colnames(out$atac_counts), out$metacell_meta$metacell_id)
  expect_true(all(out$membership$metacell_id %in% out$metacell_meta$metacell_id))
  expect_setequal(out$post_filter_diagnostics$n_usable_metacells, c(9L, 10L, 11L))
})

test_that("README formal workflow example uses current metacell API", {
  readme <- paste(readLines(test_path("..", "..", "README.md"), warn = FALSE), collapse = "\n")
  expect_false(grepl('link_stratum_cols\\s*=\\s*"cell_type"', readme))
  example_args <- c("min_cells_pre_metacell", "min_metacells_post_metacell")
  expect_true(all(example_args %in% names(formals(rc_run_regcompass_multiome_metacell))))
})
