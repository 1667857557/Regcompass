.rc_fragment_input_enabled <- function(fragment_files) {
  if (is.null(fragment_files) || identical(fragment_files, FALSE)) return(FALSE)
  if (is.logical(fragment_files)) {
    stop(
      "`fragment_files` must be FALSE/NULL, fragment path(s), or a fragment mapping data.frame; TRUE is not a file input.",
      call. = FALSE
    )
  }
  TRUE
}

.rc_fragment_aggregation_defaults <- function() {
  list(
    workers = NULL,
    rows_per_chunk = 10000000L,
    bgzip_path = NULL,
    tabix_path = NULL,
    process_n = 2000L,
    call_peaks = TRUE,
    macs2_path = NULL,
    effective_genome_size = NULL,
    peak_calling_args = list()
  )
}

.rc_resolve_fragment_aggregation_args <- function(fragment_args = list()) {
  if (!is.list(fragment_args)) {
    stop("`metacell_args$fragment_args` must be a list.", call. = FALSE)
  }
  unknown <- setdiff(names(fragment_args), names(.rc_fragment_aggregation_defaults()))
  if (length(unknown)) {
    stop(
      "Unknown `metacell_args$fragment_args`: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  out <- modifyList(.rc_fragment_aggregation_defaults(), fragment_args)
  if (!is.null(out$workers)) {
    out$workers <- as.integer(out$workers)
    if (length(out$workers) != 1L || is.na(out$workers) || out$workers < 1L) {
      stop("`fragment_args$workers` must be NULL or one positive integer.",
           call. = FALSE)
    }
  }
  out$rows_per_chunk <- as.integer(out$rows_per_chunk)
  out$process_n <- as.integer(out$process_n)
  if (length(out$rows_per_chunk) != 1L || is.na(out$rows_per_chunk) ||
      out$rows_per_chunk < 1L) {
    stop("`fragment_args$rows_per_chunk` must be one positive integer.",
         call. = FALSE)
  }
  if (length(out$process_n) != 1L || is.na(out$process_n) ||
      out$process_n < 1L) {
    stop("`fragment_args$process_n` must be one positive integer.",
         call. = FALSE)
  }
  if (!is.logical(out$call_peaks) || length(out$call_peaks) != 1L ||
      is.na(out$call_peaks)) {
    stop("`fragment_args$call_peaks` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(out$peak_calling_args)) {
    stop("`fragment_args$peak_calling_args` must be a list.", call. = FALSE)
  }
  out
}

.rc_validate_fragment_file_paths <- function(paths) {
  paths <- trimws(as.character(paths))
  if (!length(paths) || anyNA(paths) || any(!nzchar(paths))) {
    stop("Fragment paths must be non-empty character values.", call. = FALSE)
  }
  if (anyDuplicated(paths)) {
    stop("Duplicated fragment paths are not allowed.", call. = FALSE)
  }
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Fragment files are missing: ",
      paste(utils::head(missing, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  paths
}

.rc_fragment_membership_from_explicit_manifest <- function(
    fragment_files, membership) {
  required <- c("fragment_file", "object_cell", "fragment_barcode")
  missing <- setdiff(required, colnames(fragment_files))
  if (length(missing)) {
    stop(
      "Fragment mapping data.frame is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  manifest <- fragment_files[, required, drop = FALSE]
  manifest[] <- lapply(manifest, function(x) trimws(as.character(x)))
  invalid <- !stats::complete.cases(manifest) |
    !nzchar(manifest$fragment_file) |
    !nzchar(manifest$object_cell) |
    !nzchar(manifest$fragment_barcode)
  if (any(invalid)) {
    stop("Fragment mapping contains missing or empty values.", call. = FALSE)
  }
  manifest$fragment_file <- .rc_validate_fragment_file_paths(
    unique(manifest$fragment_file)
  )[match(manifest$fragment_file, unique(manifest$fragment_file))]
  unknown_cells <- setdiff(manifest$object_cell, membership$cell_id)
  if (length(unknown_cells)) {
    stop(
      "Fragment mapping contains cells absent from the Stage 2 cell set: ",
      paste(utils::head(unknown_cells, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  manifest$metacell_id <- membership$metacell_id[
    match(manifest$object_cell, membership$cell_id)
  ]
  split_rows <- split(seq_len(nrow(manifest)), manifest$fragment_file)
  lapply(names(split_rows), function(path) {
    rows <- split_rows[[path]]
    one <- manifest[rows, , drop = FALSE]
    barcode_groups <- split(one$metacell_id, one$fragment_barcode)
    ambiguous <- names(barcode_groups)[vapply(
      barcode_groups,
      function(x) length(unique(as.character(x))) != 1L,
      logical(1)
    )]
    if (length(ambiguous)) {
      stop(
        "Within one fragment file, raw barcodes map to multiple metacells: ",
        paste(utils::head(ambiguous, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    map <- stats::setNames(
      as.character(one$metacell_id),
      as.character(one$fragment_barcode)
    )
    map <- map[!duplicated(names(map))]
    list(
      fragment_file = path,
      membership = map,
      source = "explicit_fragment_manifest"
    )
  })
}

.rc_fragment_membership_from_paths <- function(fragment_files, membership) {
  paths <- .rc_validate_fragment_file_paths(fragment_files)
  labels <- names(fragment_files)
  if (length(paths) == 1L &&
      (is.null(labels) || is.na(labels[[1L]]) || !nzchar(labels[[1L]]))) {
    return(list(list(
      fragment_file = paths[[1L]],
      membership = stats::setNames(
        as.character(membership$metacell_id),
        as.character(membership$cell_id)
      ),
      source = "single_fragment_file"
    )))
  }
  if (is.null(labels) || length(labels) != length(paths) ||
      anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels)) {
    stop(
      "Multiple fragment files must be a named character vector whose names are the exact prefixes used before `_` in Seurat cell IDs, or use an explicit mapping data.frame with `fragment_file`, `object_cell`, and `fragment_barcode`.",
      call. = FALSE
    )
  }
  labels <- trimws(labels)
  object_cells <- as.character(membership$cell_id)
  lapply(seq_along(paths), function(i) {
    prefix <- labels[[i]]
    prefix_token <- paste0(prefix, "_")
    keep <- startsWith(object_cells, prefix_token)
    if (!any(keep)) {
      stop(
        "No Stage 2 cells use fragment prefix `", prefix,
        "_`; use an explicit fragment mapping data.frame if cell IDs use a different convention.",
        call. = FALSE
      )
    }
    raw_barcode <- substring(object_cells[keep], nchar(prefix_token) + 1L)
    if (any(!nzchar(raw_barcode)) || anyDuplicated(raw_barcode)) {
      stop(
        "Fragment prefix `", prefix,
        "` does not yield unique non-empty raw barcodes.",
        call. = FALSE
      )
    }
    list(
      fragment_file = paths[[i]],
      membership = stats::setNames(
        as.character(membership$metacell_id[keep]), raw_barcode
      ),
      source = "named_fragment_prefix"
    )
  })
}

.rc_resolve_fragment_memberships <- function(fragment_files, membership) {
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership))) {
    stop("Metacell membership must contain `cell_id` and `metacell_id`.",
         call. = FALSE)
  }
  if (anyNA(membership$cell_id) || any(!nzchar(as.character(membership$cell_id))) ||
      anyDuplicated(as.character(membership$cell_id)) ||
      anyNA(membership$metacell_id) ||
      any(!nzchar(as.character(membership$metacell_id)))) {
    stop("Metacell membership contains invalid cell or metacell IDs.",
         call. = FALSE)
  }
  if (is.data.frame(fragment_files)) {
    return(.rc_fragment_membership_from_explicit_manifest(
      fragment_files, membership
    ))
  }
  if (!is.character(fragment_files)) {
    stop(
      "`fragment_files` must be a character path vector or an explicit mapping data.frame.",
      call. = FALSE
    )
  }
  .rc_fragment_membership_from_paths(fragment_files, membership)
}

.rc_require_fragment_aggregation_api <- function() {
  if (!requireNamespace("SuperCell", quietly = TRUE)) {
    stop("Package 'SuperCell' is required for fragment aggregation.", call. = FALSE)
  }
  exports <- getNamespaceExports("SuperCell")
  if (!"AggregateFragmentFile" %in% exports) {
    stop(
      "Installed SuperCell lacks `AggregateFragmentFile()`. Install the current 1667857557/SuperCell_Seurat_V4 release.",
      call. = FALSE
    )
  }
  fun <- getExportedValue("SuperCell", "AggregateFragmentFile")
  required <- c(
    "input_file", "membership", "output_name", "output_path", "tmp_path",
    "nb_cl", "nb_row_split", "bgzip_path", "tabix_path",
    "return_details", "keep_unmatched"
  )
  missing <- setdiff(required, names(formals(fun)))
  if (length(missing)) {
    stop(
      "Installed SuperCell fragment API is incompatible; missing formals: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  fun
}

.rc_aggregate_single_cell_fragments <- function(
    fragment_files, membership, outdir, fragment_args = list()) {
  specs <- .rc_resolve_fragment_memberships(fragment_files, membership)
  args <- .rc_resolve_fragment_aggregation_args(fragment_args)
  aggregate_fun <- .rc_require_fragment_aggregation_api()
  fragment_dir <- file.path(outdir, "fragments")
  dir.create(fragment_dir, recursive = TRUE, showWarnings = FALSE)
  details <- vector("list", length(specs))
  manifests <- vector("list", length(specs))
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    output_name <- sprintf("metacell_fragments_%03d.tsv.gz", i)
    tmp_path <- tempfile(
      pattern = sprintf("fragment_aggregate_%03d_", i),
      tmpdir = fragment_dir
    )
    dir.create(tmp_path, recursive = TRUE, showWarnings = FALSE)
    one <- tryCatch(
      aggregate_fun(
        input_file = spec$fragment_file,
        membership = spec$membership,
        output_name = output_name,
        output_path = fragment_dir,
        tmp_path = tmp_path,
        nb_cl = args$workers,
        nb_row_split = args$rows_per_chunk,
        bgzip_path = args$bgzip_path,
        tabix_path = args$tabix_path,
        returnOutputFileName = TRUE,
        return_details = TRUE,
        keep_unmatched = FALSE
      ),
      finally = unlink(tmp_path, recursive = TRUE, force = TRUE)
    )
    if (!is.list(one) || is.null(one$fragment_file)) {
      stop("SuperCell fragment aggregation returned an incompatible result.",
           call. = FALSE)
    }
    output_file <- normalizePath(
      as.character(one$fragment_file[[1L]]), mustWork = TRUE
    )
    index_file <- paste0(output_file, ".tbi")
    if (!file.exists(index_file)) {
      stop("Aggregated fragment index was not created: ", index_file,
           call. = FALSE)
    }
    metacell_ids <- unique(as.character(one$metacell_ids %||% character()))
    metacell_ids <- metacell_ids[!is.na(metacell_ids) & nzchar(metacell_ids)]
    if (!length(metacell_ids)) {
      stop(
        "Fragment file `", spec$fragment_file,
        "` produced no metacell fragments; verify barcode mapping.",
        call. = FALSE
      )
    }
    manifests[[i]] <- data.frame(
      source_fragment_file = normalizePath(spec$fragment_file, mustWork = TRUE),
      fragment_file = output_file,
      index_file = index_file,
      object_cell = metacell_ids,
      fragment_barcode = metacell_ids,
      mapping_source = spec$source,
      stringsAsFactors = FALSE
    )
    details[[i]] <- c(
      one,
      list(
        source_fragment_file = normalizePath(spec$fragment_file, mustWork = TRUE),
        mapping_source = spec$source
      )
    )
  }
  manifest <- do.call(rbind, manifests)
  rownames(manifest) <- NULL
  list(
    fragment_manifest = manifest,
    fragment_files = unique(as.character(manifest$fragment_file)),
    details = details,
    args = args
  )
}

.rc_validate_fragment_recount_manifest <- function(
    fragment_manifest, object_cells, require_complete = TRUE) {
  if (!is.data.frame(fragment_manifest) || !nrow(fragment_manifest)) {
    if (isTRUE(require_complete)) {
      stop("No metacell fragment manifest was supplied for ATAC recounting.",
           call. = FALSE)
    }
    return(data.frame())
  }
  required <- c("fragment_file", "object_cell", "fragment_barcode")
  missing <- setdiff(required, colnames(fragment_manifest))
  if (length(missing)) {
    stop(
      "Fragment recount manifest is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  manifest <- fragment_manifest
  manifest$fragment_file <- as.character(manifest$fragment_file)
  manifest$object_cell <- as.character(manifest$object_cell)
  manifest$fragment_barcode <- as.character(manifest$fragment_barcode)
  manifest <- manifest[
    !is.na(manifest$fragment_file) & nzchar(manifest$fragment_file) &
      !is.na(manifest$object_cell) & nzchar(manifest$object_cell) &
      !is.na(manifest$fragment_barcode) & nzchar(manifest$fragment_barcode),
    , drop = FALSE
  ]
  manifest <- manifest[manifest$object_cell %in% object_cells, , drop = FALSE]
  if (!nrow(manifest)) {
    stop("Fragment recount manifest contains no mappings for the metacell object.",
         call. = FALSE)
  }
  key <- paste(manifest$fragment_file, manifest$object_cell, sep = "\001")
  barcode_n <- tapply(manifest$fragment_barcode, key,
                      function(x) length(unique(x)))
  if (any(barcode_n > 1L)) {
    stop(
      "Fragment recount manifest maps one metacell to multiple barcodes within the same fragment file.",
      call. = FALSE
    )
  }
  manifest <- unique(manifest)
  files <- unique(manifest$fragment_file)
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop(
      "Metacell fragment files are missing: ",
      paste(utils::head(missing_files, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  missing_indexes <- vapply(
    files,
    function(path) {
      !file.exists(paste0(path, ".tbi")) &&
        !file.exists(paste0(path, ".csi"))
    },
    logical(1)
  )
  if (any(missing_indexes)) {
    stop(
      "Metacell fragment indexes are missing: ",
      paste(utils::head(files[missing_indexes], 10L), collapse = ", "),
      call. = FALSE
    )
  }
  if (isTRUE(require_complete)) {
    missing_cells <- setdiff(object_cells, unique(manifest$object_cell))
    if (length(missing_cells)) {
      stop(
        "Fragment recount manifest does not cover metacells: ",
        paste(utils::head(missing_cells, 10L), collapse = ", "),
        call. = FALSE
      )
    }
  }
  rownames(manifest) <- NULL
  manifest
}

.rc_align_peak_count_matrix <- function(x, feature_ids, cell_ids) {
  if (is.null(dim(x)) || is.null(rownames(x)) || is.null(colnames(x))) {
    stop("Fragment-derived peak count matrices require row and column names.",
         call. = FALSE)
  }
  out <- Matrix::Matrix(
    0,
    nrow = length(feature_ids),
    ncol = length(cell_ids),
    sparse = TRUE,
    dimnames = list(feature_ids, cell_ids)
  )
  common_features <- intersect(feature_ids, rownames(x))
  common_cells <- intersect(cell_ids, colnames(x))
  if (length(common_features) && length(common_cells)) {
    out[common_features, common_cells] <- .rc_as_sparse(
      x[common_features, common_cells, drop = FALSE]
    )
  }
  out
}
