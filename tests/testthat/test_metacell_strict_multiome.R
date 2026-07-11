test_that("formal multiome mode rejects supplied links", {
  rna <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 2, dimnames = list(c("G1", "G2"), c("mc1", "mc2"))), sparse = TRUE)
  atac <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 2, dimnames = list(c("p1", "p2"), c("mc1", "mc2"))), sparse = TRUE)
  meta <- data.frame(metacell_id = c("mc1", "mc2"), sample_id = "s1", condition = "ctrl", cell_type = "T", stringsAsFactors = FALSE)
  gpr <- data.frame(reaction_id = "R1", and_group_id = 1, gene = "G1", stringsAsFactors = FALSE)
  supplied_links <- data.frame(peak_id = "p1", gene = "G1", weight = 1, stringsAsFactors = FALSE)
  expect_error(rc_run_layer1_from_metacells(gpr_table = gpr, rna_metacell_counts = rna, atac_metacell_counts = atac, metacell_meta = meta, metacell_seurat = NULL, peak_gene_links = supplied_links), "not accepted")
})

test_that("stratum-aware link confidence uses matching link strata only", {
  p_atac <- matrix(c(1, 0, 0, 1), nrow = 2, dimnames = list(c("p1", "p2"), c("mc1", "mc2")))
  links <- data.frame(peak_id = c("p1", "p2"), gene = c("G1", "G1"), weight = 1, link_stratum = c("T", "B"), stringsAsFactors = FALSE)
  meta <- data.frame(pool_id = c("mc1", "mc2"), cell_type = c("T", "B"), stringsAsFactors = FALSE)
  out <- rc_link_confidence_by_stratum(p_atac, links, meta, link_stratum_cols = "cell_type")
  expect_equal(as.numeric(out["G1", ]), c(1, 1))
})

test_that("rc_write_tsv_gz writes non-empty gzip", {
  f <- tempfile(fileext = ".tsv.gz")
  .rc_write_tsv_gz(data.frame(a = 1:3), f)
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 20)
})

test_that("rc_logcpm handles sparse dgCMatrix", {
  m <- Matrix::rsparsematrix(100, 5, density = 0.1)
  m@x <- abs(m@x)
  colnames(m) <- paste0("p", 1:5)
  rownames(m) <- paste0("g", 1:100)
  empty <- Matrix::colSums(m) <= 0
  if (any(empty)) m[1, empty] <- 1
  expect_s4_class(rc_logcpm(m), "dgCMatrix")
})

test_that("metacell LinkPeaks requires fragments", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 3, dimnames = list(c("p1", "p2"), paste0("mc", 1:3))), sparse = TRUE)
  chrom <- Signac::CreateChromatinAssay(counts = counts, ranges = GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = c(1, 100), width = 50)))
  obj <- SeuratObject::CreateSeuratObject(counts = counts, assay = "ATAC")
  obj[["ATAC"]] <- chrom
  gpr <- data.frame(reaction_id = "R1", and_group_id = 1, gene = "G1", stringsAsFactors = FALSE)
  expect_error(rc_recompute_metacell_peak_gene_links(obj, gpr_table = gpr), "requires fragment files")
})

test_that("LinkPeaks stratum requires enough metacells", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 2, dimnames = list(c("g1", "g2"), c("mc1", "mc2"))), sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  meta <- data.frame(metacell_id = c("mc1", "mc2"), cell_type = "T", sample_id = "s1", condition = "ctrl", stringsAsFactors = FALSE)
  gpr <- data.frame(reaction_id = "R1", and_group_id = 1, gene = "G1", stringsAsFactors = FALSE)
  expect_error(rc_recompute_metacell_peak_gene_links_by_stratum(metacell_object = obj, metacell_meta = meta, gpr_table = gpr, min_metacells_for_linkpeaks = 80, on_too_few_metacells = "stop"), "Internal invariant failed")
})


test_that("fragment registration errors before relinking if fragment paths are missing", {
  expect_error(.rc_register_signac_fragments(list(), fragment_files = tempfile(fileext = ".tsv.gz")), "fragment files are missing")
})

test_that("LinkPeaks gene matching preserves expression feature case", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 2, dimnames = list(c("HK1", "PFKM"), c("mc1", "mc2"))), sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")
  expect_identical(rc_match_linkpeaks_genes(c("hk1", "pfkm", "missing"), obj, "RNA"), c("HK1", "PFKM"))
})
