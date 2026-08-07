test_that("explicit fragment manifests may include filtered-out cells but must cover Stage 2", {
  files <- c(tempfile(fileext = ".tsv.gz"), tempfile(fileext = ".tsv.gz"))
  file.create(files)
  membership <- data.frame(
    cell_id = c("keep_A", "keep_B"),
    metacell_id = c("mc1", "mc2"),
    stringsAsFactors = FALSE
  )
  manifest <- data.frame(
    fragment_file = c(files[[1L]], files[[2L]], files[[1L]]),
    object_cell = c("keep_A", "keep_B", "filtered_out"),
    fragment_barcode = c("AAAC-1", "AAAC-1", "AAAG-1"),
    stringsAsFactors = FALSE
  )

  specs <- .rc_resolve_fragment_memberships(manifest, membership)
  mapped <- unlist(lapply(specs, function(x) unname(x$membership)),
                   use.names = FALSE)
  expect_setequal(mapped, c("mc1", "mc2"))

  expect_error(
    .rc_resolve_fragment_memberships(
      manifest[manifest$object_cell != "keep_B", , drop = FALSE],
      membership
    ),
    "does not cover Stage 2 cells"
  )
})

test_that("named fragment prefixes must assign every Stage 2 cell exactly once", {
  files <- c(
    A = tempfile(fileext = ".tsv.gz"),
    A_B = tempfile(fileext = ".tsv.gz")
  )
  file.create(unname(files))
  membership <- data.frame(
    cell_id = c("A_cell1", "A_B_cell2"),
    metacell_id = c("mc1", "mc2"),
    stringsAsFactors = FALSE
  )

  expect_error(
    .rc_resolve_fragment_memberships(files, membership),
    "exactly one source"
  )
})
