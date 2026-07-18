test_that("fragment manifest expansion records metacell barcode mappings", {
  manifest <- data.frame(
    fragment_file = c("f1.tsv.gz", "f2.tsv.gz"),
    index_file = c("f1.tsv.gz.tbi", "f2.tsv.gz.tbi"),
    stringsAsFactors = FALSE
  )

  out <- .rc_expand_fragment_manifest(manifest, c("mc1", "mc2"))
  file_counts <- table(out$fragment_file)

  expect_equal(nrow(out), 4L)
  expect_setequal(out$object_cell, c("mc1", "mc2"))
  expect_identical(out$object_cell, out$fragment_barcode)
  expect_setequal(names(file_counts), c("f1.tsv.gz", "f2.tsv.gz"))
  expect_equal(as.integer(file_counts), c(2L, 2L))
})

test_that("peak matrices are aligned to the requested feature and cell order", {
  x <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("p2", "p1"), c("mc2", "mc1"))
    ),
    sparse = TRUE
  )

  out <- .rc_align_peak_count_matrix(
    x,
    feature_ids = c("p1", "p2", "p3"),
    cell_ids = c("mc1", "mc2", "mc3")
  )

  expect_identical(rownames(out), c("p1", "p2", "p3"))
  expect_identical(colnames(out), c("mc1", "mc2", "mc3"))
  expect_equal(as.matrix(out[c("p1", "p2"), c("mc1", "mc2")]), matrix(
    c(4, 3, 2, 1),
    nrow = 2,
    dimnames = list(c("p1", "p2"), c("mc1", "mc2"))
  ))
  expect_equal(sum(out["p3", ]), 0)
  expect_equal(sum(out[, "mc3"]), 0)
})

test_that("fragment-derived de novo peaks replace the ATAC assay before Pando", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  rna <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("mc1", "mc2"))
    ),
    sparse = TRUE
  )
  atac <- Matrix::Matrix(
    matrix(
      0,
      nrow = 2,
      ncol = 2,
      dimnames = list(
        c("chr1-1-10", "chr1-20-30"),
        c("mc1", "mc2")
      )
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = rna)
  object[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = atac,
    sep = c("-", "-")
  )

  files <- c(tempfile(fileext = ".tsv.gz"), tempfile(fileext = ".tsv.gz"))
  file.create(files)
  file.create(paste0(files, ".tbi"))
  manifest <- do.call(rbind, lapply(files, function(path) {
    data.frame(
      fragment_file = path,
      object_cell = c("mc1", "mc2"),
      fragment_barcode = c("mc1", "mc2"),
      stringsAsFactors = FALSE
    )
  }))
  called_peaks <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(
      start = c(5L, 40L),
      end = c(15L, 50L)
    )
  )
  called_peak_ids <- c("chr1-5-15", "chr1-40-50")

  create_fragment_fun <- function(path, cells, validate.fragments) {
    list(path = path, cells = cells)
  }
  call_peaks_fun <- function(object, effective.genome.size, ...) {
    expect_setequal(object, files)
    expect_equal(effective.genome.size, 1.87e9)
    called_peaks
  }
  feature_matrix_fun <- function(
      fragments, features, keep_all_features, cells, process_n, verbose) {
    expect_identical(.rc_peak_ids(features), called_peak_ids)
    path <- fragments[[1L]]$path
    values <- if (identical(path, files[[1L]])) {
      matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(called_peak_ids, c("mc1", "mc2"))
      )
    } else {
      matrix(
        c(10, 20, 30, 40),
        nrow = 2,
        dimnames = list(called_peak_ids, c("mc1", "mc2"))
      )
    }
    Matrix::Matrix(values, sparse = TRUE)
  }

  out <- .rc_recount_atac_from_fragment_manifest(
    object = object,
    fragment_manifest = manifest,
    atac_assay = "ATAC",
    effective_genome_size = 1.87e9,
    create_fragment_fun = create_fragment_fun,
    feature_matrix_fun = feature_matrix_fun,
    call_peaks_fun = call_peaks_fun
  )
  counts <- .rc_get_assay_counts_safe(out, "ATAC")

  expect_identical(rownames(counts), called_peak_ids)
  expect_false(any(c("chr1-1-10", "chr1-20-30") %in% rownames(counts)))
  expect_equal(
    as.matrix(counts),
    matrix(
      c(11, 22, 33, 44),
      nrow = 2,
      dimnames = list(called_peak_ids, c("mc1", "mc2"))
    )
  )
  expect_identical(
    .rc_peak_ids(methods::slot(out[["ATAC"]], "ranges")),
    called_peak_ids
  )
  expect_identical(
    out@misc$atac_count_source,
    "recomputed_from_metacell_fragments"
  )
  expect_identical(
    out@misc$atac_peak_source,
    "de_novo_macs2_from_metacell_fragments"
  )
  expect_identical(
    out@misc$atac_fragment_recount$fragment_registration,
    "not_registered_overlapping_fragment_files"
  )
  expect_equal(out$nCount_ATAC, c(33, 77))
  expect_equal(out$nFeature_ATAC, c(2, 2))
})

test_that("metacell API exposes fragment peak-calling controls", {
  args <- names(formals(rc_make_supercell2_metacells))
  expect_true(all(c(
    "call_peaks_from_fragments",
    "macs2_path",
    "peak_calling_effective_genome_size",
    "peak_calling_args"
  ) %in% args))
})
