test_that("metacell metadata and validation enforce unique sample-aware units", {
  membership <- data.frame(
    cell_id = paste0("c", 1:4),
    metacell_id = c("mc1", "mc1", "mc2", "mc2"),
    sample_id = "s1",
    condition = "ctrl",
    cell_type = "T",
    stringsAsFactors = FALSE
  )
  meta <- rc_build_metacell_metadata(membership)
  expect_equal(meta$n_cells, c(2L, 2L))
  counts <- Matrix::Matrix(matrix(1:6, nrow = 3, dimnames = list(paste0("g", 1:3), c("mc1", "mc2"))), sparse = TRUE)
  expect_true(rc_validate_metacell_inputs(counts, meta))
  bad_meta <- rbind(meta, meta[1, ])
  expect_error(rc_validate_metacell_inputs(counts, bad_meta), "unique")
})

test_that("metacell layer1 runs from raw metacell counts", {
  counts <- Matrix::Matrix(
    matrix(c(10, 1, 5, 5, 1, 10), nrow = 3, dimnames = list(c("HK1", "HK2", "PFKM"), c("mc1", "mc2"))),
    sparse = TRUE
  )
  meta <- data.frame(metacell_id = c("mc1", "mc2"), sample_id = c("s1", "s1"), condition = c("ctrl", "ctrl"), cell_type = c("T", "T"), n_cells = c(30L, 25L), stringsAsFactors = FALSE)
  gpr <- data.frame(reaction_id = c("R_HEX", "R_HEX", "R_PFK"), and_group_id = c(1, 2, 1), gene = c("HK1", "HK2", "PFKM"), stringsAsFactors = FALSE)
  atac <- Matrix::Matrix(matrix(c(10, 1, 1, 10), nrow = 2, dimnames = list(c("p1", "p2"), c("mc1", "mc2"))), sparse = TRUE)
  links <- data.frame(peak_id = c("p1", "p2"), gene = c("HK1", "PFKM"), weight = 1, link_stratum = "T", stringsAsFactors = FALSE)
  out <- rc_run_layer1_from_metacells(gpr_table = gpr, rna_metacell_counts = counts, metacell_meta = meta, atac_metacell_counts = atac, metacell_seurat = NULL, peak_gene_links = links, allow_supplied_links = TRUE, force_metacell_relink = FALSE, bootstrap = FALSE)
  expect_equal(colnames(out$C_raw), c("mc1", "mc2"))
  expect_equal(out$layer1_unit, "metacell")
  expect_true(all(out$metacell_meta$metacell_id %in% c("mc1", "mc2")))
})

test_that("metacell sample summary reports metacell diagnostics", {
  score <- matrix(c(1, 3, 2, 4), nrow = 2, dimnames = list(c("r1", "r2"), c("mc1", "mc2")))
  meta <- data.frame(metacell_id = c("mc1", "mc2"), sample_id = "s1", condition = "ctrl", cell_type = "T", n_cells = c(20L, 30L), stringsAsFactors = FALSE)
  out <- rc_metacell_sample_summary(score, meta, condition_col = "condition")
  expect_true(all(c("n_metacells_used", "single_metacell_group_flag") %in% colnames(out)))
  expect_equal(unique(out$n_metacells_used), 2L)
  expect_false(unique(out$single_metacell_group_flag))
})

test_that("fragment_files are normalized for SuperCell", {
  x <- .rc_normalize_fragment_files("fragments.tsv.gz", atac_assay = "ATAC")
  expect_true(is.list(x))
  expect_identical(names(x), "ATAC")
  expect_identical(x[["ATAC"]], "fragments.tsv.gz")

  y <- .rc_normalize_fragment_files(c(ATAC = "a.tsv.gz", Peaks = "b.tsv.gz"), atac_assay = "ATAC")
  expect_true(is.list(y))
  expect_identical(names(y), c("ATAC", "Peaks"))

  expect_error(.rc_normalize_fragment_files(c("a.tsv.gz", "b.tsv.gz")), "named")
})

test_that("fragment registration validates one cell vector per fragment", {
  frag1 <- tempfile(fileext = ".tsv.gz")
  frag2 <- tempfile(fileext = ".tsv.gz")
  file.create(frag1, paste0(frag1, ".tbi"), frag2, paste0(frag2, ".tbi"))
  expect_error(
    .rc_register_signac_fragments(list(), fragment_files = c(frag1, frag2), cells_by_fragment = list("mc1")),
    "one cell vector per fragment file"
  )
})
