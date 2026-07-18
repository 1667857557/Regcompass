from pathlib import Path


def replace_between(path, start_marker, end_marker, replacement):
    text = path.read_text()
    start = text.find(start_marker)
    end = text.find(end_marker, start + len(start_marker))
    if start < 0 or end < 0:
        raise RuntimeError(f"Could not locate helper boundaries in {path}")
    path.write_text(text[:start] + replacement + text[end:])


metacell = Path("R/metacell_fragments.R")
replace_between(
    metacell,
    ".rc_restore_metacell_metadata <- function(\n",
    "rc_load_or_merge_metacell_objects <- function(\n",
    r'''.rc_restore_metacell_metadata <- function(
    object_meta, metacell_meta, expected_ids) {
  if (!is.data.frame(object_meta) || !is.data.frame(metacell_meta)) {
    stop("Metacell metadata inputs must be data.frames.", call. = FALSE)
  }
  if (!"metacell_id" %in% colnames(metacell_meta)) {
    stop("`metacell_meta` must contain `metacell_id`.", call. = FALSE)
  }
  if (anyDuplicated(colnames(metacell_meta))) {
    stop("`metacell_meta` column names must be unique.", call. = FALSE)
  }
  expected_ids <- as.character(expected_ids)
  supplied_ids <- as.character(metacell_meta$metacell_id)
  if (anyNA(expected_ids) || any(!nzchar(expected_ids)) ||
      anyDuplicated(expected_ids)) {
    stop("Expected metacell IDs must be non-missing and unique.", call. = FALSE)
  }
  if (anyNA(supplied_ids) || any(!nzchar(supplied_ids)) ||
      anyDuplicated(supplied_ids)) {
    stop("`metacell_meta$metacell_id` must be non-missing and unique.",
         call. = FALSE)
  }
  missing_metadata <- setdiff(expected_ids, supplied_ids)
  if (length(missing_metadata)) {
    stop("Metacell metadata are missing expected IDs: ",
         paste(utils::head(missing_metadata, 10L), collapse = ", "),
         call. = FALSE)
  }
  missing_object_rows <- setdiff(expected_ids, rownames(object_meta))
  if (length(missing_object_rows)) {
    stop("Seurat metadata are missing expected metacell IDs: ",
         paste(utils::head(missing_object_rows, 10L), collapse = ", "),
         call. = FALSE)
  }
  aligned <- metacell_meta[
    match(expected_ids, supplied_ids), , drop = FALSE
  ]
  rownames(aligned) <- expected_ids
  object_meta <- object_meta[expected_ids, , drop = FALSE]
  supplied_names <- colnames(aligned)
  extra_columns <- setdiff(colnames(object_meta), supplied_names)
  for (column in extra_columns) {
    originals <- supplied_names[
      supplied_names != column &
        make.names(supplied_names, unique = FALSE) == column
    ]
    redundant_alias <- any(vapply(originals, function(original) {
      left <- object_meta[[column]]
      right <- aligned[[original]]
      if (is.factor(left)) left <- as.character(left)
      if (is.factor(right)) right <- as.character(right)
      isTRUE(all.equal(unname(left), unname(right), check.attributes = FALSE))
    }, logical(1)))
    if (!redundant_alias) aligned[[column]] <- object_meta[[column]]
  }
  aligned
}

'''
)

pando = Path("R/pando_grn.R")
replace_between(
    pando,
    ".rc_remap_projection_metadata <- function(\n",
    "rc_run_pando_meta_modules <- function(\n",
    r'''.rc_remap_projection_metadata <- function(
    x, group_meta, sample_col, display_cols) {
  if (!is.data.frame(x) || !is.data.frame(group_meta)) {
    stop("Projection and group metadata must be data.frames.", call. = FALSE)
  }
  if (!nrow(x)) return(x)
  sample_col <- as.character(sample_col)
  if (length(sample_col) != 1L || is.na(sample_col) || !nzchar(sample_col)) {
    stop("`sample_col` must name one metadata column.", call. = FALSE)
  }
  if (!"sample_id" %in% colnames(x)) {
    stop("Projection table is missing column: sample_id", call. = FALSE)
  }
  required_group <- c("group_id", sample_col)
  missing_group <- setdiff(required_group, colnames(group_meta))
  if (length(missing_group)) {
    stop("Projection group metadata are missing columns: ",
         paste(missing_group, collapse = ", "), call. = FALSE)
  }
  group_meta$group_id <- as.character(group_meta$group_id)
  if (anyNA(group_meta$group_id) || any(!nzchar(group_meta$group_id)) ||
      anyDuplicated(group_meta$group_id)) {
    stop("Projection group metadata require unique, non-empty group IDs.",
         call. = FALSE)
  }
  group_ids <- as.character(x$sample_id)
  if (anyNA(group_ids) || any(!nzchar(group_ids))) {
    stop("Projection group IDs must be non-missing and non-empty.",
         call. = FALSE)
  }
  group_index <- match(group_ids, group_meta$group_id)
  if (anyNA(group_index)) {
    missing_groups <- unique(group_ids[is.na(group_index)])
    stop("Projection metadata are missing group IDs: ",
         paste(utils::head(missing_groups, 10L), collapse = ", "),
         call. = FALSE)
  }
  x$sample_id <- NULL
  x$group_id <- group_ids
  group_columns <- setdiff(colnames(group_meta), "group_id")
  collisions <- intersect(group_columns, colnames(x))
  if (length(collisions)) {
    stop("Projection and group metadata contain conflicting columns: ",
         paste(collisions, collapse = ", "), call. = FALSE)
  }
  for (column in group_columns) {
    x[[column]] <- group_meta[[column]][group_index]
  }
  sample_values <- as.character(x[[sample_col]])
  if (anyNA(sample_values) || any(!nzchar(sample_values))) {
    stop("Projection metadata produced missing or empty sample IDs.",
         call. = FALSE)
  }
  x$sample_id <- sample_values
  display_cols <- unique(as.character(display_cols))
  missing_display <- setdiff(display_cols, colnames(x))
  if (length(missing_display)) {
    stop("Projection display columns are missing: ",
         paste(missing_display, collapse = ", "), call. = FALSE)
  }
  rownames(x) <- NULL
  x[, c(display_cols, setdiff(colnames(x), display_cols)), drop = FALSE]
}

'''
)

tests = Path("tests/testthat/test-metadata-contracts.R")
text = tests.read_text()
text = text.replace(
    'X_rc_condition = c("wrong", "wrong")',
    'X_rc_condition = c("control", "treated")',
    1
)
append = r'''

test_that("distinct metadata columns are not discarded as aliases", {
  object_meta <- data.frame(
    X_rc_condition = c("legacy_a", "legacy_b"),
    row.names = c("MC1", "MC2"),
    check.names = FALSE
  )
  metacell_meta <- data.frame(
    metacell_id = c("MC1", "MC2"),
    `_rc_condition` = c("control", "treated"),
    check.names = FALSE
  )
  restored <- .rc_restore_metacell_metadata(
    object_meta, metacell_meta, c("MC1", "MC2")
  )
  expect_identical(restored$X_rc_condition, c("legacy_a", "legacy_b"))
  expect_identical(restored$`_rc_condition`, c("control", "treated"))
})

test_that("metacell merge restores exact names on the Seurat object", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(1, nrow = 2, ncol = 2,
           dimnames = list(c("g1", "g2"), c("MC1", "MC2"))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object@meta.data$X_rc_condition <- c("control", "treated")
  metadata <- data.frame(
    metacell_id = c("MC2", "MC1"),
    `_rc_condition` = c("treated", "control"),
    sample_id = c("S2", "S1"),
    check.names = FALSE
  )
  restored <- rc_load_or_merge_metacell_objects(
    list(object), metacell_meta = metadata
  )
  expect_identical(colnames(restored), c("MC2", "MC1"))
  expect_true("_rc_condition" %in% colnames(restored@meta.data))
  expect_false("X_rc_condition" %in% colnames(restored@meta.data))
  expect_identical(restored@meta.data$`_rc_condition`, c("treated", "control"))
  expect_identical(restored@meta.data$sample_id, c("S2", "S1"))
})

test_that("Pando projection remapping rejects payload collisions", {
  projection <- data.frame(
    sample_id = "G1", condition = "projection", gene = "A",
    stringsAsFactors = FALSE
  )
  group_meta <- data.frame(
    group_id = "G1", sample_id = "S1", condition = "control",
    stringsAsFactors = FALSE
  )
  expect_error(
    .rc_remap_projection_metadata(
      projection, group_meta, "sample_id",
      c("group_id", "sample_id", "condition")
    ),
    "conflicting columns: condition"
  )
})
'''
if 'distinct metadata columns are not discarded as aliases' not in text:
    text += append
tests.write_text(text)
print("Applied robust metadata contract refinements")
