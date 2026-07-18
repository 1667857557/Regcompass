test_that("authoritative metacell metadata restores non-syntactic names", {
  object_meta <- data.frame(
    X_rc_condition = c("wrong", "wrong"),
    nCount_RNA = c(10, 20),
    row.names = c("MC1", "MC2"),
    check.names = FALSE
  )
  metacell_meta <- data.frame(
    metacell_id = c("MC2", "MC1"),
    `_rc_condition` = c("treated", "control"),
    sample_id = c("S2", "S1"),
    check.names = FALSE
  )
  restored <- .rc_restore_metacell_metadata(
    object_meta,
    metacell_meta,
    expected_ids = c("MC1", "MC2")
  )
  expect_identical(rownames(restored), c("MC1", "MC2"))
  expect_true("_rc_condition" %in% colnames(restored))
  expect_false("X_rc_condition" %in% colnames(restored))
  expect_identical(restored$`_rc_condition`, c("control", "treated"))
  expect_identical(restored$sample_id, c("S1", "S2"))
  expect_identical(restored$nCount_RNA, c(10, 20))
})

test_that("metacell metadata contract rejects duplicate IDs", {
  object_meta <- data.frame(row.names = c("MC1", "MC2"))
  metacell_meta <- data.frame(
    metacell_id = c("MC1", "MC1"),
    check.names = FALSE
  )
  expect_error(
    .rc_restore_metacell_metadata(
      object_meta, metacell_meta, c("MC1", "MC2")
    ),
    "non-missing and unique"
  )
})

test_that("Pando projection remapping avoids sample_id suffix collisions", {
  projection <- data.frame(
    sample_id = c("G2", "G1"),
    gene = c("B", "A"),
    module_id = c("G2::GRN0001", "G1::GRN0001"),
    stringsAsFactors = FALSE
  )
  group_meta <- data.frame(
    group_id = c("G1", "G2"),
    sample_id = c("S1", "S2"),
    `_rc_condition` = c("control", "treated"),
    check.names = FALSE
  )
  remapped <- .rc_remap_projection_metadata(
    projection,
    group_meta,
    sample_col = "sample_id",
    display_cols = c("group_id", "sample_id", "_rc_condition")
  )
  expect_identical(remapped$group_id, c("G2", "G1"))
  expect_identical(remapped$sample_id, c("S2", "S1"))
  expect_identical(remapped$`_rc_condition`, c("treated", "control"))
  expect_false(any(c("sample_id.x", "sample_id.y") %in% colnames(remapped)))
})

test_that("Pando projection remapping supports custom sample columns", {
  projection <- data.frame(
    sample_id = "condition|donor|celltype",
    gene = "A",
    module_id = "condition|donor|celltype::GRN0001",
    stringsAsFactors = FALSE
  )
  group_meta <- data.frame(
    group_id = "condition|donor|celltype",
    donor = "D1",
    condition = "control",
    stringsAsFactors = FALSE
  )
  remapped <- .rc_remap_projection_metadata(
    projection,
    group_meta,
    sample_col = "donor",
    display_cols = c("group_id", "donor", "condition")
  )
  expect_identical(remapped$sample_id, "D1")
  expect_identical(remapped$donor, "D1")
})

test_that("Pando projection remapping rejects non-unique group maps", {
  projection <- data.frame(
    sample_id = "G1", gene = "A", module_id = "G1::GRN0001",
    stringsAsFactors = FALSE
  )
  group_meta <- data.frame(
    group_id = c("G1", "G1"), sample_id = c("S1", "S2"),
    stringsAsFactors = FALSE
  )
  expect_error(
    .rc_remap_projection_metadata(
      projection, group_meta, "sample_id", c("group_id", "sample_id")
    ),
    "unique, non-empty group IDs"
  )
})
