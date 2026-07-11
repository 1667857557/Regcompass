test_that("strict stratum helpers reject non-unique and missing metadata definitions", {
  expect_identical(.rc_strict_stratum_cols("sample_id", "condition", "cell_type"), c("condition", "sample_id", "cell_type"))
  expect_error(.rc_strict_stratum_cols("sample_id", "sample_id", "cell_type"), "three valid")

  meta <- data.frame(condition = c("ctrl", "ctrl"), sample_id = c("s1", ""), cell_type = c("T", "T"), row.names = c("c1", "c2"), stringsAsFactors = FALSE)
  expect_error(.rc_add_stratum_id(meta, c("condition", "sample_id", "cell_type")), "missing or empty")
})

test_that("pre-metacell filter uses only observed cell count", {
  skip_if_not_installed("SeuratObject")
  ids <- c(paste0("a", seq_len(99)), paste0("b", seq_len(100)), paste0("c", seq_len(146)))
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = length(ids), dimnames = list(c("g1", "g2"), ids)), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$condition <- "ctrl"
  object$sample_id <- c(rep("s99", 99), rep("s100", 100), rep("s146", 146))
  object$cell_type <- "gdT"

  out <- rc_filter_pre_metacell_strata(object, min_cells = 100L)

  expect_false("ctrl|s99|gdT" %in% out$retained_strata)
  expect_true("ctrl|s100|gdT" %in% out$retained_strata)
  expect_true("ctrl|s146|gdT" %in% out$retained_strata)
  expect_false(any(c("expected_metacells", "expected_post_metacell_failure") %in% colnames(out$diagnostics)))
  expect_equal(out$diagnostics$status[out$diagnostics$stratum_id == "ctrl|s100|gdT"], "retained_for_metacell")
})

test_that("post-metacell filter uses actual metacell count boundaries and subsets bundle", {
  ids9 <- paste0("mc9_", seq_len(9))
  ids10 <- paste0("mc10_", seq_len(10))
  ids11 <- paste0("mc11_", seq_len(11))
  ids10_low <- paste0("mc10low_", seq_len(10))
  ids <- c(ids9, ids10, ids11, ids10_low)
  meta <- data.frame(
    metacell_id = ids,
    condition = "ctrl",
    sample_id = rep(c("s9", "s10", "s11", "s10low"), c(9, 10, 11, 10)),
    cell_type = "T",
    n_cells = c(rep(20L, 30), rep(5L, 10)),
    stringsAsFactors = FALSE
  )
  mc <- list(
    metacell_meta = meta,
    rna_counts = Matrix::Matrix(matrix(1, nrow = 2, ncol = length(ids), dimnames = list(c("g1", "g2"), ids)), sparse = TRUE),
    atac_counts = Matrix::Matrix(matrix(1, nrow = 2, ncol = length(ids), dimnames = list(c("p1", "p2"), ids)), sparse = TRUE),
    membership = data.frame(cell_id = paste0("c", seq_along(ids)), metacell_id = ids, stringsAsFactors = FALSE),
    fragment_manifest = data.frame(fragment_file = "frag.tsv.gz", object_cell = ids, fragment_barcode = ids, stringsAsFactors = FALSE)
  )

  out <- rc_filter_post_metacell_strata(mc, min_metacells = 10L)

  expect_false(any(grepl("^mc9_", out$metacell_meta$metacell_id)))
  expect_true(all(c(ids10, ids11, ids10_low) %in% out$metacell_meta$metacell_id))
  expect_identical(colnames(out$rna_counts), out$metacell_meta$metacell_id)
  expect_identical(colnames(out$atac_counts), out$metacell_meta$metacell_id)
  expect_true(all(out$membership$metacell_id %in% out$metacell_meta$metacell_id))
  expect_true(all(out$fragment_manifest$object_cell %in% out$metacell_meta$metacell_id))
  expect_setequal(out$post_filter_diagnostics$n_metacells, c(9L, 10L, 11L))
  expect_false(any(c("n_usable_metacells", "n_low_power_metacells") %in% colnames(out$post_filter_diagnostics)))
})

test_that("README formal workflow example uses current metacell API", {
  readme <- paste(readLines(test_path("..", "..", "README.md"), warn = FALSE), collapse = "\n")
  expect_false(grepl('link_stratum_cols\\s*=\\s*"cell_type"', readme))
  example_args <- c("min_cells_pre_metacell", "min_metacells_post_metacell")
  expect_true(all(example_args %in% names(formals(rc_run_regcompass_multiome_metacell))))
})

test_that("formal workflow contains a LinkPeaks stratum invariant before relinking", {
  txt <- paste(deparse(body(rc_run_regcompass_multiome_metacell)), collapse = "\n")
  expect_match(txt, "bad_link_strata", fixed = TRUE)
  expect_match(txt, "post-filtered LinkPeaks strata contain fewer than", fixed = TRUE)
  expect_match(txt, "link_stratum_cols = strict_cols", fixed = TRUE)
})
