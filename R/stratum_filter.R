.rc_strict_stratum_cols <- function(sample_col, condition_col, celltype_col) {
  cols <- unique(c(condition_col, sample_col, celltype_col))
  cols <- cols[!is.null(cols) & !is.na(cols) & nzchar(cols)]
  if (length(cols) != 3L) {
    stop("`condition_col`, `sample_col`, and `celltype_col` must identify three valid metadata columns.", call. = FALSE)
  }
  cols
}

.rc_add_stratum_id <- function(meta, cols) {
  if (!is.data.frame(meta)) stop("`meta` must be a data.frame.", call. = FALSE)
  missing <- setdiff(cols, colnames(meta))
  if (length(missing)) stop("Missing strict stratum columns: ", paste(missing, collapse = ", "), call. = FALSE)
  bad <- vapply(meta[, cols, drop = FALSE], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1))
  if (any(bad)) stop("Strict stratum columns contain missing or empty values: ", paste(cols[bad], collapse = ", "), call. = FALSE)
  meta$.rc_stratum_id <- as.character(interaction(meta[, cols, drop = FALSE], sep = "|", drop = TRUE, lex.order = TRUE))
  meta
}

#' Filter strict strata before metacell construction
#' @export
rc_filter_pre_metacell_strata <- function(object,
                                          sample_col = "sample_id",
                                          condition_col = "condition",
                                          celltype_col = "cell_type",
                                          min_cells = 100L) {
  min_cells <- suppressWarnings(as.integer(min_cells[[1L]]))
  if (!is.finite(min_cells) || min_cells < 1L) stop("`min_cells` must be a positive integer.", call. = FALSE)
  cols <- .rc_strict_stratum_cols(sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col)
  meta <- .rc_add_stratum_id(object@meta.data, cols)
  strata <- split(rownames(meta), meta$.rc_stratum_id)
  diag <- do.call(rbind, lapply(names(strata), function(id) {
    cells <- strata[[id]]
    idx <- match(cells[[1L]], rownames(meta))
    vals <- meta[idx, cols, drop = FALSE]
    n <- length(cells)
    eligible <- n >= min_cells
    data.frame(vals,
               stratum_id = id,
               n_cells = n,
               min_required_cells = min_cells,
               eligible = eligible,
               status = if (eligible) "retained_for_metacell" else paste0("excluded_below_", min_cells, "_cells"),
               stringsAsFactors = FALSE,
               check.names = FALSE)
  }))
  keep_strata <- as.character(diag$stratum_id[diag$eligible])
  keep_cells <- rownames(meta)[meta$.rc_stratum_id %in% keep_strata]
  excluded_cells <- data.frame(cell_id = rownames(meta)[!meta$.rc_stratum_id %in% keep_strata],
                               meta[!meta$.rc_stratum_id %in% keep_strata, c(cols, ".rc_stratum_id"), drop = FALSE],
                               row.names = NULL,
                               check.names = FALSE)
  names(excluded_cells)[names(excluded_cells) == ".rc_stratum_id"] <- "stratum_id"
  if (!length(keep_cells)) stop("No condition × sample × cell-type strata contain at least ", min_cells, " cells.", call. = FALSE)
  list(object = subset(object, cells = keep_cells), diagnostics = diag, excluded_cells = excluded_cells, retained_strata = keep_strata, stratum_cols = cols)
}

.rc_subset_metacell_bundle <- function(mc, keep_metacell_ids) {
  keep_metacell_ids <- unique(as.character(keep_metacell_ids))
  if (!length(keep_metacell_ids)) stop("No metacells remain after filtering.", call. = FALSE)
  meta <- mc$metacell_meta
  meta$metacell_id <- as.character(meta$metacell_id)
  keep_metacell_ids <- keep_metacell_ids[keep_metacell_ids %in% meta$metacell_id]
  meta <- meta[match(keep_metacell_ids, meta$metacell_id), , drop = FALSE]
  out <- mc
  out$metacell_meta <- meta
  out$metacell_meta_used <- meta
  out$used_metacell_ids <- keep_metacell_ids
  if (!is.null(out$rna_counts)) out$rna_counts <- out$rna_counts[, keep_metacell_ids, drop = FALSE]
  if (!is.null(out$atac_counts)) out$atac_counts <- out$atac_counts[, keep_metacell_ids, drop = FALSE]
  if (is.data.frame(out$membership) && "metacell_id" %in% colnames(out$membership)) {
    out$membership <- out$membership[as.character(out$membership$metacell_id) %in% keep_metacell_ids, , drop = FALSE]
    out$membership_used <- out$membership
  }
  if (is.data.frame(out$fragment_manifest) && nrow(out$fragment_manifest) && "object_cell" %in% colnames(out$fragment_manifest)) {
    out$fragment_manifest <- out$fragment_manifest[as.character(out$fragment_manifest$object_cell) %in% keep_metacell_ids, , drop = FALSE]
    out$fragment_manifest_used <- out$fragment_manifest
    out$fragment_files <- unique(as.character(out$fragment_manifest$fragment_file))
  }
  if (!is.null(out$rna_counts) && !identical(colnames(out$rna_counts), out$metacell_meta$metacell_id)) stop("RNA counts are not aligned to retained metacell metadata.", call. = FALSE)
  if (!is.null(out$atac_counts) && !identical(colnames(out$atac_counts), out$metacell_meta$metacell_id)) stop("ATAC counts are not aligned to retained metacell metadata.", call. = FALSE)
  if (is.data.frame(out$membership) && "metacell_id" %in% colnames(out$membership) && any(!as.character(out$membership$metacell_id) %in% out$metacell_meta$metacell_id)) stop("Membership contains metacells absent from retained metadata.", call. = FALSE)
  if (is.data.frame(out$fragment_manifest) && nrow(out$fragment_manifest) && "object_cell" %in% colnames(out$fragment_manifest) && any(!as.character(out$fragment_manifest$object_cell) %in% out$metacell_meta$metacell_id)) stop("Fragment manifest contains metacells absent from retained metadata.", call. = FALSE)
  out
}

#' Filter strict strata after metacell construction
#' @export
rc_filter_post_metacell_strata <- function(mc,
                                           sample_col = "sample_id",
                                           condition_col = "condition",
                                           celltype_col = "cell_type",
                                           min_metacells = 10L) {
  min_metacells <- suppressWarnings(as.integer(min_metacells[[1L]]))
  if (!is.finite(min_metacells) || min_metacells < 1L) stop("`min_metacells` must be a positive integer.", call. = FALSE)
  cols <- .rc_strict_stratum_cols(sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col)
  meta <- .rc_add_stratum_id(mc$metacell_meta, cols)
  if (!"metacell_id" %in% colnames(meta)) stop("`mc$metacell_meta` must contain `metacell_id`.", call. = FALSE)
  meta$metacell_id <- as.character(meta$metacell_id)
  metacell_tab <- table(meta$.rc_stratum_id)
  diag <- unique(meta[, c(cols, ".rc_stratum_id"), drop = FALSE])
  diag$n_cells_input <- as.integer(tapply(meta$n_cells, meta$.rc_stratum_id, sum, na.rm = TRUE)[diag$.rc_stratum_id])
  diag$n_metacells <- as.integer(metacell_tab[diag$.rc_stratum_id])
  diag$min_required_metacells <- min_metacells
  diag$eligible <- diag$n_metacells >= min_metacells
  diag$status <- ifelse(diag$eligible, "retained_for_linkpeaks", paste0("excluded_below_", min_metacells, "_metacells"))
  if ("effective_gamma" %in% colnames(meta)) diag$effective_gamma <- as.integer(tapply(meta$effective_gamma, meta$.rc_stratum_id, function(x) x[which(!is.na(x))[1L]] %||% NA_integer_)[diag$.rc_stratum_id])
  names(diag)[names(diag) == ".rc_stratum_id"] <- "stratum_id"
  keep_strata <- as.character(diag$stratum_id[diag$eligible])
  keep_ids <- as.character(meta$metacell_id[meta$.rc_stratum_id %in% keep_strata])
  excluded <- meta[!meta$metacell_id %in% keep_ids, , drop = FALSE]
  names(excluded)[names(excluded) == ".rc_stratum_id"] <- "stratum_id"
  if (!length(keep_ids)) stop("No condition × sample × cell-type strata contain at least ", min_metacells, " metacells.", call. = FALSE)
  out <- .rc_subset_metacell_bundle(mc, keep_metacell_ids = keep_ids)
  out$post_filter_diagnostics <- diag
  out$excluded_post_metacells <- excluded
  out$retained_strata <- keep_strata
  out
}

rc_write_stratum_filter_reports <- function(diagnostics, excluded_units, stage, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(diagnostics, file.path(outdir, paste0(stage, "_strata.tsv.gz")))
  excluded_name <- if (identical(stage, "pre_metacell")) "excluded_pre_metacell_cells.tsv.gz" else "excluded_post_metacells.tsv.gz"
  .rc_write_tsv_gz(excluded_units, file.path(outdir, excluded_name))
  invisible(TRUE)
}
