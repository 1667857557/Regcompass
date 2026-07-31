#' Validate metacell-level RegCompass inputs
rc_validate_metacell_inputs <- function(rna_metacell_counts,
                                        metacell_meta,
                                        atac_metacell_counts = NULL,
                                        metacell_id_col = "metacell_id",
                                        sample_col = NULL,
                                        condition_col = "condition",
                                        celltype_col = "cell_type") {
  if (is.null(dim(rna_metacell_counts)) || length(dim(rna_metacell_counts)) != 2L) stop("`rna_metacell_counts` must be a feature-by-metacell matrix.", call. = FALSE)
  if (is.null(colnames(rna_metacell_counts)) || anyNA(colnames(rna_metacell_counts)) || any(!nzchar(colnames(rna_metacell_counts)))) stop("`rna_metacell_counts` must have metacell IDs in colnames().", call. = FALSE)
  if (!is.data.frame(metacell_meta)) stop("`metacell_meta` must be a data.frame.", call. = FALSE)
  required <- c(metacell_id_col, sample_col, condition_col, celltype_col)
  required <- required[!is.null(required) & !is.na(required) & nzchar(required)]
  missing <- setdiff(required, colnames(metacell_meta))
  if (length(missing) > 0L) stop("`metacell_meta` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyNA(metacell_meta[[metacell_id_col]]) || anyDuplicated(as.character(metacell_meta[[metacell_id_col]]))) stop("Metacell IDs must be non-missing and unique.", call. = FALSE)
  missing_mc <- setdiff(colnames(rna_metacell_counts), as.character(metacell_meta[[metacell_id_col]]))
  if (length(missing_mc) > 0L) stop("`metacell_meta` is missing metadata for metacells: ", paste(utils::head(missing_mc, 5L), collapse = ", "), call. = FALSE)
  if (!is.null(atac_metacell_counts)) {
    if (is.null(dim(atac_metacell_counts)) || length(dim(atac_metacell_counts)) != 2L) stop("`atac_metacell_counts` must be a feature-by-metacell matrix.", call. = FALSE)
    rna_ids <- as.character(colnames(rna_metacell_counts))
    atac_ids <- as.character(colnames(atac_metacell_counts))
    if (!setequal(rna_ids, atac_ids)) stop("RNA and ATAC metacell matrices contain different metacell IDs.", call. = FALSE)
    if (!identical(rna_ids, atac_ids)) stop("RNA and ATAC metacell matrices contain the same IDs but in different order.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Construct RegCompass stratum IDs
#'
#' This helper exposes the same `interaction(..., sep = "|", lex.order = TRUE)`
#' convention used internally for strict strata.
#'
#' @param meta Metadata data frame.
#' @param cols Metadata columns to combine.
#' @param sep Separator used between column values.
#' @return Character vector with one stratum ID per row of `meta`.
rc_make_stratum_id <- function(meta, cols, sep = "|") {
  if (!is.data.frame(meta)) stop("`meta` must be a data.frame.", call. = FALSE)
  cols <- cols[!is.null(cols) & !is.na(cols) & nzchar(cols)]
  if (!length(cols)) stop("`cols` must contain at least one metadata column.", call. = FALSE)
  missing <- setdiff(cols, colnames(meta))
  if (length(missing)) stop("Missing stratum columns: ", paste(missing, collapse = ", "), call. = FALSE)
  bad <- vapply(meta[, cols, drop = FALSE], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1))
  if (any(bad)) stop("Stratum columns contain missing or empty values: ", paste(cols[bad], collapse = ", "), call. = FALSE)
  as.character(interaction(meta[, cols, drop = FALSE], sep = sep, drop = TRUE, lex.order = TRUE))
}

#' Build one-row-per-metacell metadata from membership
rc_build_metacell_metadata <- function(membership,
                                       metacell_id_col = "metacell_id",
                                       cell_id_col = "cell_id") {
  if (!is.data.frame(membership)) {
    stop("`membership` must be a data.frame.", call. = FALSE)
  }
  missing <- setdiff(
    c(metacell_id_col, cell_id_col),
    colnames(membership)
  )
  if (length(missing)) {
    stop("`membership` is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  metacell_id <- trimws(as.character(membership[[metacell_id_col]]))
  keep <- !is.na(metacell_id) & nzchar(metacell_id)
  x <- membership[keep, , drop = FALSE]
  x[[metacell_id_col]] <- metacell_id[keep]
  if (!nrow(x)) {
    out <- x[, setdiff(colnames(x), cell_id_col), drop = FALSE]
    out$n_cells <- integer()
    return(out)
  }

  strict_columns <- intersect(c("condition", "cell_type"), colnames(x))
  split_rows <- split(seq_len(nrow(x)), x[[metacell_id_col]])
  for (metacell in names(split_rows)) {
    rows <- split_rows[[metacell]]
    for (column in strict_columns) {
      values <- trimws(as.character(x[[column]][rows]))
      values <- unique(values[!is.na(values) & nzchar(values)])
      if (length(values) != 1L || anyNA(x[[column]][rows]) ||
          any(!nzchar(trimws(as.character(x[[column]][rows]))))) {
        stop(
          "Metacell `", metacell, "` mixes metadata or contains missing values in `",
          column, "`.",
          call. = FALSE
        )
      }
    }
  }

  columns <- setdiff(colnames(x), cell_id_col)
  out <- x[!duplicated(x[[metacell_id_col]]), columns, drop = FALSE]
  out$n_cells <- as.integer(vapply(
    as.character(out[[metacell_id_col]]),
    function(id) length(split_rows[[id]]),
    integer(1)
  ))
  rownames(out) <- NULL
  out
}

.rc_safe_path_component <- function(x) {
  gsub("[^A-Za-z0-9_.=-]+", "_", as.character(x))
}

.rc_write_tsv_gz <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(file, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x, file = con, sep = "\t", quote = FALSE, row.names = FALSE
  )
  invisible(file)
}

.rc_get_assay_counts_safe <- function(object, assay) {
  .rc_get_assay_counts(object, assay)
}
