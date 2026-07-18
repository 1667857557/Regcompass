from pathlib import Path


def replace_between(path, start_marker, end_marker, replacement):
    text = path.read_text()
    start = text.find(start_marker)
    if start < 0:
        raise RuntimeError(f"Missing start marker in {path}: {start_marker}")
    end = text.find(end_marker, start)
    if end < 0:
        raise RuntimeError(f"Missing end marker in {path}: {end_marker}")
    path.write_text(text[:start] + replacement + text[end:])


def replace_once(path, old, new):
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Expected one match in {path}, found {count}")
    path.write_text(text.replace(old, new, 1))


metacell_path = Path("R/metacell_fragments.R")
metacell_helper = r'''.rc_restore_metacell_metadata <- function(
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
    original_names <- supplied_names[
      supplied_names != column &
        make.names(supplied_names, unique = FALSE) == column
    ]
    redundant_alias <- any(vapply(original_names, function(original_name) {
      left <- object_meta[[column]]
      right <- aligned[[original_name]]
      if (is.factor(left)) left <- as.character(left)
      if (is.factor(right)) right <- as.character(right)
      isTRUE(all.equal(
        unname(left), unname(right),
        check.attributes = FALSE
      ))
    }, logical(1)))
    if (!redundant_alias) aligned[[column]] <- object_meta[[column]]
  }
  aligned
}

'''
replace_between(
    metacell_path,
    ".rc_restore_metacell_metadata <- function(\n",
    "rc_load_or_merge_metacell_objects <- function(\n",
    metacell_helper,
)

old_metacell_call = r'''    extra_in_object <- setdiff(observed, expected)
    obj <- subset(obj, cells = expected)
    if (!identical(colnames(obj), expected)) {
      stop("Merged object could not be subset and reordered to expected ",
           "metacell IDs.", call. = FALSE)
    }
    obj@meta.data <- .rc_restore_metacell_metadata(
      object_meta = obj@meta.data,
      metacell_meta = metacell_meta,
      expected_ids = expected
    )
    attr(obj, "removed_extra_metacell_ids") <- extra_in_object
'''
new_metacell_call = r'''    extra_in_object <- setdiff(observed, expected)
    restored_meta <- .rc_restore_metacell_metadata(
      object_meta = obj@meta.data,
      metacell_meta = metacell_meta,
      expected_ids = expected
    )
    obj <- subset(obj, cells = expected)
    if (!identical(colnames(obj), expected)) {
      stop("Merged object could not be subset and reordered to expected ",
           "metacell IDs.", call. = FALSE)
    }
    obj@meta.data <- restored_meta
    attr(obj, "removed_extra_metacell_ids") <- extra_in_object
'''
replace_once(metacell_path, old_metacell_call, new_metacell_call)

pando_path = Path("R/pando_grn.R")
pando_helper = r'''.rc_remap_projection_metadata <- function(
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
replace_between(
    pando_path,
    ".rc_remap_projection_metadata <- function(\n",
    "rc_run_pando_meta_modules <- function(\n",
    pando_helper,
)

Path("tests/testthat/test-metadata-contracts.R").write_text(r'''test_that("authoritative metacell metadata restores non-syntactic names", {
  object_meta <- data.frame(
    X_rc_condition = c("control", "treated"),
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

test_that("metacell merge restores exact names on the Seurat object", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      1,
      nrow = 2,
      ncol = 2,
      dimnames = list(c("g1", "g2"), c("MC1", "MC2"))
    ),
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
    list(object),
    metacell_meta = metadata
  )
  expect_identical(colnames(restored), c("MC2", "MC1"))
  expect_true("_rc_condition" %in% colnames(restored@meta.data))
  expect_false("X_rc_condition" %in% colnames(restored@meta.data))
  expect_identical(
    restored@meta.data$`_rc_condition`,
    c("treated", "control")
  )
  expect_identical(restored@meta.data$sample_id, c("S2", "S1"))
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
''')

print("Refined metadata contract fixes and tests")
