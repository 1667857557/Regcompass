test_that("fragment cell maps preserve object-to-file barcode direction", {
  map <- data.frame(
    object_cell = c("T_MC1", "T_MC2"),
    fragment_barcode = c("MC1", "MC2"),
    stringsAsFactors = FALSE
  )

  out <- .rc_normalize_fragment_cell_map(map, object_cells = c("T_MC1", "T_MC2"))

  expect_identical(names(out), c("T_MC1", "T_MC2"))
  expect_identical(unname(out), c("MC1", "MC2"))
})

test_that("fragment registration rejects overlapping object cells", {
  files <- c(tempfile(fileext = ".tsv.gz"), tempfile(fileext = ".tsv.gz"))
  file.create(files)
  file.create(paste0(files, ".tbi"))

  maps <- list(
    c(MC1 = "MC1", MC2 = "MC2"),
    c(MC2 = "MC2", MC3 = "MC3")
  )

  expect_error(
    .rc_validate_fragment_registration_plan(
      files,
      maps,
      object_cells = c("MC1", "MC2", "MC3"),
      require_complete = TRUE
    ),
    "multiple fragment files"
  )
})

test_that("fragment registration validates complete object coverage", {
  files <- c(tempfile(fileext = ".tsv.gz"), tempfile(fileext = ".tsv.gz"))
  file.create(files)
  file.create(paste0(files, ".tbi"))

  expect_error(
    .rc_validate_fragment_registration_plan(
      files,
      list(c(MC1 = "MC1"), c(MC2 = "MC2")),
      object_cells = c("MC1", "MC2", "MC3"),
      require_complete = TRUE
    ),
    "does not exactly cover"
  )

  expect_true(
    .rc_validate_fragment_registration_plan(
      files,
      list(c(MC1 = "MC1"), c(MC2 = "MC2", MC3 = "MC3")),
      object_cells = c("MC1", "MC2", "MC3"),
      require_complete = TRUE
    )
  )
})

test_that("fragment manifest registration merges rows for the same fragment path", {
  f <- tempfile(fileext = ".tsv.gz")
  manifest <- data.frame(
    fragment_file = c(f, f),
    object_cell = c("A_MC1", "A_MC2"),
    fragment_barcode = c("MC1", "MC2"),
    stringsAsFactors = FALSE
  )

  reg <- .rc_fragment_registration_from_manifest(
    manifest,
    object_cells = c("A_MC1", "A_MC2")
  )

  expect_identical(reg$fragment_files, f)
  expect_length(reg$cell_maps, 1L)
  expect_identical(names(reg$cell_maps[[1L]]), c("A_MC1", "A_MC2"))
  expect_identical(unname(reg$cell_maps[[1L]]), c("MC1", "MC2"))
})

test_that("metacell merge rejects duplicate IDs before Seurat renaming", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      1,
      nrow = 2,
      ncol = 1,
      dimnames = list(c("g1", "g2"), "MC1")
    ),
    sparse = TRUE
  )
  obj1 <- SeuratObject::CreateSeuratObject(counts = counts)
  obj2 <- SeuratObject::CreateSeuratObject(counts = counts)

  expect_error(
    rc_load_or_merge_metacell_objects(list(obj1, obj2)),
    "not globally unique"
  )
})

test_that("fragment manifest filters extra mappings but requires complete used mappings", {
  manifest <- data.frame(
    fragment_file = "frag.tsv.gz",
    object_cell = c("mc1", "mc2", "mc3"),
    fragment_barcode = c("mc1", "mc2", "mc3"),
    stringsAsFactors = FALSE
  )

  reg <- .rc_fragment_registration_from_manifest(
    manifest,
    object_cells = c("mc1", "mc2")
  )
  expect_identical(names(reg$cell_maps[[1L]]), c("mc1", "mc2"))

  expect_error(
    .rc_fragment_registration_from_manifest(
      manifest[1, , drop = FALSE],
      object_cells = c("mc1", "mc2")
    ),
    "missing mappings"
  )
})

test_that("merged objects are subset before strict ID validation", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      1,
      nrow = 2,
      ncol = 3,
      dimnames = list(c("g1", "g2"), c("mc1", "mc2", "mc3"))
    ),
    sparse = TRUE
  )
  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  meta <- data.frame(
    metacell_id = c("mc1", "mc2"),
    sample_id = "s1",
    condition = "ctrl",
    cell_type = "T",
    stringsAsFactors = FALSE
  )

  out <- rc_load_or_merge_metacell_objects(
    list(obj),
    metacell_meta = meta
  )

  expect_identical(colnames(out), c("mc1", "mc2"))
  expect_identical(attr(out, "removed_extra_metacell_ids"), "mc3")
})

test_that("missing expected metacell IDs remain fatal", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      1,
      nrow = 2,
      ncol = 1,
      dimnames = list(c("g1", "g2"), "mc1")
    ),
    sparse = TRUE
  )
  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  meta <- data.frame(
    metacell_id = c("mc1", "mc2"),
    sample_id = "s1",
    condition = "ctrl",
    cell_type = "T",
    stringsAsFactors = FALSE
  )

  expect_error(
    rc_load_or_merge_metacell_objects(
      list(obj),
      metacell_meta = meta
    ),
    "mc2"
  )
})


test_that("fragment processing is selected by an explicit pure contract", {
  expect_false(.rc_fragments_requested(FALSE))
  expect_true(.rc_fragments_requested(NULL))
  expect_false(.rc_fragment_processing_enabled(FALSE, TRUE, "regcompass"))
  expect_false(.rc_fragment_processing_enabled(NULL, FALSE, "regcompass"))
  expect_false(.rc_fragment_processing_enabled(NULL, TRUE, "none"))
  expect_true(.rc_fragment_processing_enabled(NULL, TRUE, "regcompass"))
  expect_false(.rc_fragment_processing_enabled(
    NULL, TRUE, "regcompass", require_fragments = FALSE
  ))
})

test_that("obsolete fragment-file merge argument is absent", {
  expect_false("fragment_files" %in%
                 names(formals(rc_load_or_merge_metacell_objects)))
})
