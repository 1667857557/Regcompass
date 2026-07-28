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

test_that("ConditionGRNFit extraction writes metadata without sample remapping", {
  implementation <- paste(
    deparse(body(.rc_extract_condition_grn_contract_without_comparison_guard)),
    collapse = "\n"
  )
  bridge <- paste(
    deparse(body(.rc_extract_condition_grn_contract)), collapse = "\n"
  )
  expect_match(
    implementation,
    "tab[[condition_col]] <- condition",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "tab[[celltype_col]] <- fit$cell_type",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "tab$group_id <- rc_make_stratum_id(",
    fixed = TRUE
  )
  expect_match(
    bridge,
    ".rc_extract_condition_grn_contract_without_comparison_guard",
    fixed = TRUE
  )
  expect_match(
    bridge,
    ".rc_apply_condition_comparison_semantics",
    fixed = TRUE
  )
  expect_false(exists(".rc_remap_projection_metadata", inherits = TRUE))
})
