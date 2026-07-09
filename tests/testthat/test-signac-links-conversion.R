test_that("Signac link data frames are normalized", {
  links <- data.frame(
    seqnames = "chr1",
    start = 10,
    end = 20,
    gene_name = "hk1",
    correlation = 0.4,
    stringsAsFactors = FALSE
  )
  out <- rc_signac_links_to_peak_gene_table(links, score_col = "correlation")
  expect_equal(out$peak_id, "chr1-10-20")
  expect_equal(out$gene, "HK1")
  expect_equal(out$weight, 0.4)
})

test_that("empty Signac link data frames return an empty normalized table", {
  links <- data.frame(gene = character(), score = numeric(), peak = character())
  out <- rc_signac_links_to_peak_gene_table(links)
  expect_equal(colnames(out), c("peak_id", "gene", "weight", "pvalue"))
  expect_equal(nrow(out), 0)
})
