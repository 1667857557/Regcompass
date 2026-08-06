.rc_expand_fragment_manifest <- function(manifest, metacell_ids) {
  if (!is.data.frame(manifest) || !nrow(manifest)) return(manifest)
  metacell_ids <- unique(trimws(as.character(metacell_ids)))
  metacell_ids <- metacell_ids[!is.na(metacell_ids) & nzchar(metacell_ids)]
  if (!length(metacell_ids)) {
    stop("No metacell IDs were available for fragment manifest expansion.",
         call. = FALSE)
  }
  if (all(c("object_cell", "fragment_barcode") %in% colnames(manifest))) {
    return(manifest)
  }
  pieces <- lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[rep(i, length(metacell_ids)), , drop = FALSE]
    row$object_cell <- metacell_ids
    row$fragment_barcode <- metacell_ids
    row
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
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
    stop("Fragment recount manifest is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  manifest <- fragment_manifest[, unique(c(required, colnames(fragment_manifest))),
                                drop = FALSE]
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
  conflicts <- names(barcode_n)[barcode_n > 1L]
  if (length(conflicts)) {
    stop("Fragment recount manifest maps one metacell to multiple barcodes ",
         "within the same fragment file.", call. = FALSE)
  }
  manifest <- unique(manifest)
  files <- unique(manifest$fragment_file)
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop("Metacell fragment files are missing: ",
         paste(utils::head(missing_files, 10L), collapse = ", "),
         call. = FALSE)
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
    stop("Metacell fragment indexes are missing: ",
         paste(utils::head(files[missing_indexes], 10L), collapse = ", "),
         call. = FALSE)
  }
  if (isTRUE(require_complete)) {
    missing_cells <- setdiff(object_cells, unique(manifest$object_cell))
    if (length(missing_cells)) {
      stop("Fragment recount manifest does not cover metacells: ",
           paste(utils::head(missing_cells, 10L), collapse = ", "),
           call. = FALSE)
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

.rc_clear_signac_fragments <- function(object, atac_assay = "ATAC") {
  if (!inherits(object, "Seurat")) return(object)
  if (!requireNamespace("Signac", quietly = TRUE)) return(object)
  if (!atac_assay %in% names(object@assays)) return(object)
  if (!inherits(object[[atac_assay]], "ChromatinAssay")) return(object)
  fragment_setter <- get("Fragments<-", envir = asNamespace("Signac"))
  object[[atac_assay]] <- fragment_setter(
    object[[atac_assay]], value = list()
  )
  object
}

.rc_normalize_fragment_cell_map <- function(
    cell_map, object_cells, fragment_file = NULL) {
  if (is.data.frame(cell_map)) {
    required <- c("object_cell", "fragment_barcode")
    missing <- setdiff(required, colnames(cell_map))
    if (length(missing)) {
      stop("`cell_map` is missing columns: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    cell_map <- stats::setNames(
      as.character(cell_map$fragment_barcode),
      as.character(cell_map$object_cell)
    )
  } else {
    cell_map <- as.character(cell_map)
    if (is.null(names(cell_map))) names(cell_map) <- cell_map
  }
  if (!length(cell_map)) {
    stop(
      "Fragment cell map is empty",
      if (!is.null(fragment_file)) paste0(": ", fragment_file) else ".",
      call. = FALSE
    )
  }
  if (anyNA(cell_map) || any(!nzchar(cell_map)) ||
      anyNA(names(cell_map)) || any(!nzchar(names(cell_map)))) {
    stop("Fragment cell map contains missing or empty identifiers.",
         call. = FALSE)
  }
  if (anyDuplicated(names(cell_map))) {
    duplicated_cells <- unique(names(cell_map)[duplicated(names(cell_map))])
    stop("Duplicated object cells within one fragment mapping: ",
         paste(utils::head(duplicated_cells, 10L), collapse = ", "),
         call. = FALSE)
  }
  unknown <- setdiff(names(cell_map), object_cells)
  if (length(unknown)) {
    stop("Fragment mapping contains cells absent from the merged object: ",
         paste(utils::head(unknown, 10L), collapse = ", "), call. = FALSE)
  }
  cell_map
}

.rc_validate_fragment_registration_plan <- function(
    fragment_files, cell_maps, object_cells, require_complete = TRUE) {
  fragment_files <- as.character(fragment_files)
  if (length(fragment_files) != length(cell_maps)) {
    stop("`fragment_files` and `cell_maps` must have the same length.",
         call. = FALSE)
  }
  if (!length(fragment_files)) {
    stop("No fragment files were supplied.", call. = FALSE)
  }
  missing_files <- fragment_files[!file.exists(fragment_files)]
  if (length(missing_files)) {
    stop("Metacell fragment files are missing: ",
         paste(utils::head(missing_files, 10L), collapse = ", "),
         call. = FALSE)
  }
  missing_indexes <- vapply(
    fragment_files,
    function(path) {
      !file.exists(paste0(path, ".tbi")) &&
        !file.exists(paste0(path, ".csi"))
    },
    logical(1)
  )
  if (any(missing_indexes)) {
    stop("Metacell fragment tabix indexes are missing: ",
         paste(utils::head(fragment_files[missing_indexes], 10L),
               collapse = ", "), call. = FALSE)
  }
  registered <- unlist(lapply(cell_maps, names), use.names = FALSE)
  if (anyDuplicated(registered)) {
    duplicated_cells <- unique(registered[duplicated(registered)])
    stop("Object cells are assigned to multiple fragment files: ",
         paste(utils::head(duplicated_cells, 10L), collapse = ", "),
         call. = FALSE)
  }
  if (isTRUE(require_complete)) {
    missing_cells <- setdiff(object_cells, registered)
    extra_cells <- setdiff(registered, object_cells)
    if (length(missing_cells) || length(extra_cells)) {
      stop(
        "Fragment registration does not exactly cover the merged object. ",
        "Missing cells: ",
        paste(utils::head(missing_cells, 10L), collapse = ", "),
        "; extra cells: ",
        paste(utils::head(extra_cells, 10L), collapse = ", "),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.rc_fragment_registration_from_manifest <- function(
    fragment_manifest, object_cells) {
  if (!"fragment_file" %in% colnames(fragment_manifest)) {
    stop("`fragment_manifest` must contain `fragment_file`.", call. = FALSE)
  }
  manifest <- fragment_manifest
  manifest$fragment_file <- as.character(manifest$fragment_file)
  required <- c("object_cell", "fragment_barcode")
  if (!all(required %in% colnames(manifest))) {
    stop(
      "Fragment manifest entries must contain explicit `object_cell` and ",
      "`fragment_barcode` columns.",
      call. = FALSE
    )
  }
  manifest$object_cell <- as.character(manifest$object_cell)
  manifest$fragment_barcode <- as.character(manifest$fragment_barcode)
  manifest <- manifest[manifest$object_cell %in% object_cells, , drop = FALSE]
  missing_maps <- setdiff(object_cells, unique(manifest$object_cell))
  if (length(missing_maps)) {
    stop("Fragment manifest is missing mappings for metacells: ",
         paste(utils::head(missing_maps, 10L), collapse = ", "),
         call. = FALSE)
  }
  manifest <- unique(
    manifest[, c("fragment_file", "object_cell", "fragment_barcode"),
             drop = FALSE]
  )
  cell_path <- paste(
    manifest$fragment_file, manifest$object_cell, sep = "\001"
  )
  barcode_by_cell_path <- tapply(
    manifest$fragment_barcode,
    cell_path,
    function(x) length(unique(x))
  )
  conflicts <- names(barcode_by_cell_path)[barcode_by_cell_path > 1L]
  if (length(conflicts)) {
    conflict_cells <- sub("^.*\\001", "", conflicts)
    stop("Fragment manifest assigns one object cell to multiple barcodes: ",
         paste(utils::head(conflict_cells, 10L), collapse = ", "),
         call. = FALSE)
  }
  files <- unique(manifest$fragment_file)
  maps <- lapply(files, function(path) {
    x <- manifest[manifest$fragment_file == path, , drop = FALSE]
    .rc_normalize_fragment_cell_map(
      x[, c("object_cell", "fragment_barcode"), drop = FALSE],
      object_cells = object_cells,
      fragment_file = path
    )
  })
  list(fragment_files = files, cell_maps = maps)
}

