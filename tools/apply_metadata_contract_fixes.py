from pathlib import Path


def replace_once(path, old, new):
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Expected one match in {path}, found {count}")
    path.write_text(text.replace(old, new, 1))


metacell_path = Path("R/metacell_fragments.R")
metacell_marker = "rc_load_or_merge_metacell_objects <- function(\n"
metacell_helper = r'''.rc_restore_metacell_metadata <- function(
    object_meta, metacell_meta, expected_ids) {
  if (!is.data.frame(object_meta) || !is.data.frame(metacell_meta)) {
    stop("Metacell metadata inputs must be data.frames.", call. = FALSE)
  }
  if (!"metacell_id" %in% colnames(metacell_meta)) {
    stop("`metacell_meta` must contain `metacell_id`.", call. = FALSE)
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
  syntactic_aliases <- make.names(supplied_names, unique = FALSE)
  repaired_aliases <- unique(
    syntactic_aliases[syntactic_aliases != supplied_names]
  )
  extra_columns <- setdiff(
    colnames(object_meta),
    unique(c(supplied_names, repaired_aliases))
  )
  for (column in extra_columns) aligned[[column]] <- object_meta[[column]]
  aligned
}

'''
replace_once(
    metacell_path,
    metacell_marker,
    metacell_helper + metacell_marker,
)

old_metacell_block = r'''  if (!is.null(metacell_meta)) {
    metacell_meta$metacell_id <- as.character(metacell_meta$metacell_id)
    expected <- metacell_meta$metacell_id
    observed <- colnames(obj)
    missing_in_object <- setdiff(expected, observed)
    if (length(missing_in_object)) {
      stop("Merged metacell object is missing expected IDs: ",
           paste(utils::head(missing_in_object, 10L), collapse = ", "),
           call. = FALSE)
    }
    extra_in_object <- setdiff(observed, expected)
    obj <- subset(obj, cells = expected)
    if (!identical(colnames(obj), expected)) {
      stop("Merged object could not be subset and reordered to expected ",
           "metacell IDs.", call. = FALSE)
    }
    attr(obj, "removed_extra_metacell_ids") <- extra_in_object
  }
'''
new_metacell_block = r'''  if (!is.null(metacell_meta)) {
    if (!is.data.frame(metacell_meta) ||
        !"metacell_id" %in% colnames(metacell_meta)) {
      stop("`metacell_meta` must be a data.frame containing `metacell_id`.",
           call. = FALSE)
    }
    metacell_meta$metacell_id <- as.character(metacell_meta$metacell_id)
    expected <- metacell_meta$metacell_id
    observed <- colnames(obj)
    missing_in_object <- setdiff(expected, observed)
    if (length(missing_in_object)) {
      stop("Merged metacell object is missing expected IDs: ",
           paste(utils::head(missing_in_object, 10L), collapse = ", "),
           call. = FALSE)
    }
    extra_in_object <- setdiff(observed, expected)
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
  }
'''
replace_once(metacell_path, old_metacell_block, new_metacell_block)

pando_path = Path("R/pando_grn.R")
pando_marker = "rc_run_pando_meta_modules <- function(metacell_object,\n"
pando_helper = r'''.rc_remap_projection_metadata <- function(
    x, group_meta, sample_col, display_cols) {
  if (!is.data.frame(x) || !is.data.frame(group_meta)) {
    stop("Projection and group metadata must be data.frames.", call. = FALSE)
  }
  if (!nrow(x)) return(x)
  required_projection <- c("sample_id", "module_id")
  missing_projection <- setdiff(required_projection, colnames(x))
  if (length(missing_projection)) {
    stop("Projection table is missing columns: ",
         paste(missing_projection, collapse = ", "), call. = FALSE)
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
  if (".rc_projection_order" %in% colnames(x)) {
    stop("Projection table contains reserved column `.rc_projection_order`.",
         call. = FALSE)
  }
  x$group_id <- as.character(x$sample_id)
  missing_groups <- setdiff(unique(x$group_id), group_meta$group_id)
  if (length(missing_groups)) {
    stop("Projection metadata are missing group IDs: ",
         paste(utils::head(missing_groups, 10L), collapse = ", "),
         call. = FALSE)
  }
  x$sample_id <- NULL
  x$.rc_projection_order <- seq_len(nrow(x))
  out <- merge(
    x, group_meta,
    by = "group_id", all.x = TRUE, sort = FALSE
  )
  if (nrow(out) != nrow(x)) {
    stop("Projection metadata merge changed the number of rows.",
         call. = FALSE)
  }
  out <- out[order(out$.rc_projection_order), , drop = FALSE]
  out$.rc_projection_order <- NULL
  sample_values <- as.character(out[[sample_col]])
  if (anyNA(sample_values) || any(!nzchar(sample_values))) {
    stop("Projection metadata produced missing or empty sample IDs.",
         call. = FALSE)
  }
  out$sample_id <- sample_values
  rownames(out) <- NULL
  out[, c(display_cols, setdiff(colnames(out), display_cols)), drop = FALSE]
}

'''
replace_once(pando_path, pando_marker, pando_helper + pando_marker)

old_pando_block = r'''  group_meta <- unique(status_table[, c("group_id", group_cols), drop = FALSE])
  remap_projection <- function(x) {
    if (!nrow(x)) return(x)
    x$group_id <- as.character(x$sample_id)
    x <- merge(x, group_meta, by = "group_id", all.x = TRUE, sort = FALSE)
    x$sample_id <- as.character(x[[sample_col]])
    x[, c(display_cols, setdiff(colnames(x), display_cols)), drop = FALSE]
  }
  projection$nodes <- remap_projection(projection$nodes)
  projection$edges <- remap_projection(projection$edges)
'''
new_pando_block = r'''  group_meta <- unique(status_table[, c("group_id", group_cols), drop = FALSE])
  projection$nodes <- .rc_remap_projection_metadata(
    projection$nodes,
    group_meta = group_meta,
    sample_col = sample_col,
    display_cols = display_cols
  )
  projection$edges <- .rc_remap_projection_metadata(
    projection$edges,
    group_meta = group_meta,
    sample_col = sample_col,
    display_cols = display_cols
  )
'''
replace_once(pando_path, old_pando_block, new_pando_block)

test_path = Path("tests/testthat/test-metadata-contracts.R")
test_path.write_text(r'''test_that("authoritative metacell metadata restores non-syntactic names", {
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
''')

print("Applied metadata contract fixes and regression tests")
