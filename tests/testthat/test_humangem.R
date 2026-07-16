test_that("Human-GEM preparation returns symbol GPR tables", {
  repo <- file.path(tempdir(), paste0("Human-GEM-mini-", Sys.getpid()))
  model <- file.path(repo, "model")
  dir.create(model, recursive = TRUE)
  writeLines(c(
    "genes\tgeneSymbols",
    "ENSG00000156515\tHK1",
    "ENSG00000003436\tPFKM",
    "ENSG00000134333\tLDHA"
  ), file.path(model, "genes.tsv"))
  writeLines(c(
    "rxns\trxnKEGGID",
    "MAR00001\tR00001",
    "MAR00002\tR00002"
  ), file.path(model, "reactions.tsv"))
  writeLines(c(
    "reactions:",
    "- id: MAR00001",
    "  name: hexokinase",
    "  - gene_reaction_rule: ENSG00000156515 or ENSG00000003436",
    "- id: MAR00002",
    "  name: lactate dehydrogenase",
    "  - gene_reaction_rule: ENSG00000134333"
  ), file.path(model, "Human-GEM.yml"))

  out <- rc_prepare_humangem_gpr_table(repo)
  expect_true(all(c("gpr_table", "metabolic_genes", "reaction_rules", "genes", "reactions") %in% names(out)))
  expect_setequal(out$gpr_table$gene, c("HK1", "PFKM", "LDHA"))
  expect_setequal(out$metabolic_genes, c("HK1", "PFKM", "LDHA"))
  expect_true(all(c("reaction_id", "and_group_id", "gene") %in% colnames(out$gpr_table)))
})

test_that("Human-GEM preparation can keep Ensembl gene IDs", {
  repo <- file.path(tempdir(), paste0("Human-GEM-mini-ensembl-", Sys.getpid()))
  model <- file.path(repo, "model")
  dir.create(model, recursive = TRUE)
  writeLines(c("genes\tgeneSymbols", "ENSG1\tHK1"), file.path(model, "genes.tsv"))
  writeLines(c("rxns", "MAR00001"), file.path(model, "reactions.tsv"))
  writeLines(c("reactions:", "- id: MAR00001", "  - gene_reaction_rule: ENSG1"), file.path(model, "Human-GEM.yml"))

  out <- rc_prepare_humangem_gpr_table(repo, gene_format = "ensembl")
  expect_identical(out$gpr_table$gene, "ENSG1")
})

test_that("Human-GEM YAML parser supports omap/list-prefixed gene_reaction_rule", {
  yml <- tempfile(fileext = ".yml")
  writeLines(c(
    "- id: MAR00001",
    "  - gene_reaction_rule: ENSG000001 or ENSG000002"
  ), yml)
  out <- RegCompassR:::rc_read_humangem_yml_rules(yml)
  expect_equal(out$reaction_id, "MAR00001")
  expect_equal(out$gpr, "ENSG000001 or ENSG000002")
})

test_that("Human-GEM downloader prefers tags for semver and validates archives", {
  skip_if(!nzchar(Sys.which("zip")), "zip command is required to build mock archive")
  repo <- file.path(tempdir(), paste0("Human-GEM-zip-", Sys.getpid()))
  model <- file.path(repo, "model")
  dir.create(model, recursive = TRUE)
  writeLines(c("genes\tgeneSymbols", "ENSG1\tHK1"), file.path(model, "genes.tsv"))
  writeLines(c("rxns", "MAR00001"), file.path(model, "reactions.tsv"))
  writeLines(c("reactions:", "- id: MAR00001", "  - gene_reaction_rule: ENSG1"), file.path(model, "Human-GEM.yml"))
  zipfile <- tempfile(fileext = ".zip")
  old <- setwd(dirname(repo)); on.exit(setwd(old), add = TRUE)
  utils::zip(zipfile, files = basename(repo), flags = "-r9Xq")
  calls <- character()
  mock_download <- function(url, destfile, mode, quiet) {
    calls <<- c(calls, url)
    if (grepl("/tags/", url)) {
      file.copy(zipfile, destfile, overwrite = TRUE)
      return(0L)
    }
    warning("404 Not Found")
    writeLines("<html>not found</html>", destfile)
    22L
  }

  out <- rc_download_humangem_gpr_table(destdir = tempfile("hg-dl-"), ref = "v2.0.0", overwrite = TRUE, download_fun = mock_download)
  expect_match(calls[[1]], "/tags/")
  expect_equal(attr(out, "download_diagnostics")$archive_validation[[1]], "ok")
})

test_that("Human-GEM downloader rejects status-zero HTML archives and cleans part files", {
  mock_download <- function(url, destfile, mode, quiet) {
    writeLines("<html>not a zip</html>", destfile)
    0L
  }
  dest <- tempfile("hg-dl-bad-")
  expect_error(
    rc_download_humangem_gpr_table(destdir = dest, ref = "main", overwrite = TRUE, download_fun = mock_download),
    "invalid_zip_magic"
  )
  expect_false(any(grepl("\\.part$", list.files(dest, full.names = TRUE))))
})

test_that("Human-GEM downloader uses heads first for branch refs", {
  calls <- character()
  mock_download <- function(url, destfile, mode, quiet) {
    calls <<- c(calls, url)
    writeLines("not a zip", destfile)
    1L
  }
  expect_error(
    rc_download_humangem_gpr_table(
      destdir = tempfile("hg-dl-branch-"), ref = "develop",
      overwrite = TRUE, download_fun = mock_download
    ),
    "Failed to download"
  )
  expect_match(calls[[1L]], "/heads/")
  expect_match(calls[[2L]], "/tags/")
})
