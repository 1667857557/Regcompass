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

  reg <- .rc_fragment_registration_from_manifest(manifest, object_cells = c("A_MC1", "A_MC2"))

  expect_identical(reg$fragment_files, f)
  expect_length(reg$cell_maps, 1L)
  expect_identical(names(reg$cell_maps[[1L]]), c("A_MC1", "A_MC2"))
  expect_identical(unname(reg$cell_maps[[1L]]), c("MC1", "MC2"))
})

test_that("metacell merge rejects duplicate IDs before Seurat renaming", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 1, dimnames = list(c("g1", "g2"), "MC1")), sparse = TRUE)
  obj1 <- SeuratObject::CreateSeuratObject(counts = counts)
  obj2 <- SeuratObject::CreateSeuratObject(counts = counts)

  expect_error(
    rc_load_or_merge_metacell_objects(list(obj1, obj2)),
    "not globally unique"
  )
})
