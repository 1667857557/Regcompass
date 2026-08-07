test_that("named fragment paths route prefixed cells without changing metacell groups", {
  files <- c(
    Control_1 = tempfile(fileext = ".tsv.gz"),
    Cre_1 = tempfile(fileext = ".tsv.gz")
  )
  file.create(unname(files))
  membership <- data.frame(
    cell_id = c(
      "Control_1_AAAC-1", "Control_1_AAAG-1",
      "Cre_1_AAAC-1", "Cre_1_AAAT-1"
    ),
    metacell_id = c("ctrl_mc1", "ctrl_mc1", "cre_mc1", "cre_mc2"),
    stringsAsFactors = FALSE
  )

  specs <- .rc_resolve_fragment_memberships(files, membership)

  expect_length(specs, 2L)
  expect_identical(specs[[1L]]$source, "named_fragment_prefix")
  expect_identical(specs[[2L]]$source, "named_fragment_prefix")
  expect_identical(names(specs[[1L]]$membership), c("AAAC-1", "AAAG-1"))
  expect_identical(unname(specs[[1L]]$membership), c("ctrl_mc1", "ctrl_mc1"))
  expect_identical(names(specs[[2L]]$membership), c("AAAC-1", "AAAT-1"))
  expect_identical(unname(specs[[2L]]$membership), c("cre_mc1", "cre_mc2"))
})

test_that("explicit fragment manifest permits the same raw barcode in different files", {
  files <- c(tempfile(fileext = ".tsv.gz"), tempfile(fileext = ".tsv.gz"))
  file.create(files)
  membership <- data.frame(
    cell_id = c("sampleA_cell", "sampleB_cell"),
    metacell_id = c("mcA", "mcB"),
    stringsAsFactors = FALSE
  )
  manifest <- data.frame(
    fragment_file = files,
    object_cell = membership$cell_id,
    fragment_barcode = c("AAAC-1", "AAAC-1"),
    stringsAsFactors = FALSE
  )

  specs <- .rc_resolve_fragment_memberships(manifest, membership)

  expect_length(specs, 2L)
  expect_identical(names(specs[[1L]]$membership), "AAAC-1")
  expect_identical(names(specs[[2L]]$membership), "AAAC-1")
  expect_setequal(
    unlist(lapply(specs, function(x) unname(x$membership))),
    c("mcA", "mcB")
  )
})

test_that("multiple fragment paths require an explicit routing contract", {
  files <- c(tempfile(fileext = ".tsv.gz"), tempfile(fileext = ".tsv.gz"))
  file.create(files)
  membership <- data.frame(
    cell_id = c("s1_AAAC-1", "s2_AAAC-1"),
    metacell_id = c("mc1", "mc2"),
    stringsAsFactors = FALSE
  )

  expect_error(
    .rc_resolve_fragment_memberships(files, membership),
    "named character vector"
  )
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

test_that("fragment-derived peaks replace the metacell ATAC assay", {
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
      fragments, features, cells, process_n, verbose) {
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
  counts <- .rc_get_assay_counts(out, "ATAC")

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
  expect_equal(as.numeric(out$nCount_ATAC), c(33, 77))
  expect_equal(as.numeric(out$nFeature_ATAC), c(2, 2))
})
