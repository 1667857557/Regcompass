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

.rc_register_signac_fragments <- function(
    object, fragment_files = NULL, cells_by_fragment = NULL,
    atac_assay = "ATAC", replace_existing = TRUE,
    require_complete = TRUE, validate_fragments = TRUE) {
  if (is.null(fragment_files) || length(fragment_files) == 0L) return(object)
  fragment_files <- as.character(fragment_files)
  if (is.null(cells_by_fragment)) {
    if (length(fragment_files) != 1L) {
      stop("`cells_by_fragment` is required when registering multiple files.",
           call. = FALSE)
    }
    ids <- colnames(object)
    cells_by_fragment <- list(stats::setNames(ids, ids))
  }
  if (length(cells_by_fragment) != length(fragment_files)) {
    stop("`cells_by_fragment` must have one cell vector per fragment file.",
         call. = FALSE)
  }
  cell_maps <- Map(
    function(cell_map, path) {
      .rc_normalize_fragment_cell_map(
        cell_map,
        object_cells = colnames(object),
        fragment_file = path
      )
    },
    cells_by_fragment,
    fragment_files
  )
  .rc_validate_fragment_registration_plan(
    fragment_files,
    cell_maps,
    object_cells = colnames(object),
    require_complete = require_complete
  )
  if (!requireNamespace("Signac", quietly = TRUE)) {
    stop("Package 'Signac' is required to register metacell fragment files.",
         call. = FALSE)
  }
  if (!inherits(object, "Seurat") ||
      !atac_assay %in% names(object@assays)) {
    stop("Metacell object is missing ATAC assay `", atac_assay, "`.",
         call. = FALSE)
  }
  fragment_files <- normalizePath(fragment_files, mustWork = TRUE)
  fragment_setter <- get("Fragments<-", envir = asNamespace("Signac"))
  if (isTRUE(replace_existing)) {
    object[[atac_assay]] <- fragment_setter(
      object[[atac_assay]], value = list()
    )
  }
  fragments <- Map(function(path, cell_map) {
    tryCatch(
      Signac::CreateFragmentObject(
        path = path,
        cells = cell_map,
        validate.fragments = validate_fragments
      ),
      error = function(e) {
        stop("Failed to register metacell fragment file `", path, "`: ",
             conditionMessage(e), call. = FALSE)
      }
    )
  }, fragment_files, cell_maps)
  object[[atac_assay]] <- fragment_setter(
    object[[atac_assay]], value = fragments
  )
  registered <- unlist(lapply(fragments, SeuratObject::Cells),
                       use.names = FALSE)
  if (anyDuplicated(registered)) {
    stop("Post-registration validation detected cells in multiple Fragment ",
         "objects.", call. = FALSE)
  }
  object
}

rc_load_or_merge_metacell_objects <- function(
    metacell_objects, fragment_manifest = NULL, metacell_meta = NULL,
    rna_assay = "RNA", atac_assay = "ATAC",
    require_complete_fragments = TRUE) {
  if (is.null(metacell_objects) || !length(metacell_objects)) {
    stop("No metacell Seurat objects supplied.", call. = FALSE)
  }
  objs <- lapply(
    metacell_objects,
    function(x) if (inherits(x, "Seurat")) x else readRDS(x)
  )
  input_cells <- unlist(lapply(objs, colnames), use.names = FALSE)
  duplicated_before_merge <- unique(input_cells[duplicated(input_cells)])
  if (length(duplicated_before_merge)) {
    stop("Metacell IDs are not globally unique before merge: ",
         paste(utils::head(duplicated_before_merge, 10L), collapse = ", "),
         call. = FALSE)
  }
  has_fragment_manifest <- is.data.frame(fragment_manifest) &&
    nrow(fragment_manifest) > 0L
  precounted <- has_fragment_manifest && all(vapply(
    objs,
    function(x) {
      identical(
        tryCatch(x@misc$atac_count_source, error = function(e) NULL),
        "recomputed_from_metacell_fragments"
      )
    },
    logical(1)
  ))
  objs <- lapply(objs, .rc_clear_signac_fragments, atac_assay = atac_assay)
  obj <- if (length(objs) == 1L) {
    objs[[1L]]
  } else {
    Reduce(function(a, b) merge(x = a, y = b, merge.data = FALSE), objs)
  }
  if (anyDuplicated(colnames(obj))) {
    stop("Merged metacell object contains duplicated cell names.",
         call. = FALSE)
  }
  if (!is.null(metacell_meta)) {
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
  if (has_fragment_manifest && !precounted) {
    obj <- .rc_recount_atac_from_fragment_manifest(
      object = obj,
      fragment_manifest = fragment_manifest,
      atac_assay = atac_assay,
      require_complete = require_complete_fragments
    )
  } else if (has_fragment_manifest) {
    manifest <- .rc_validate_fragment_recount_manifest(
      fragment_manifest,
      object_cells = colnames(obj),
      require_complete = require_complete_fragments
    )
    mapped <- unlist(
      split(manifest$object_cell, manifest$fragment_file),
      use.names = FALSE
    )
    if (!anyDuplicated(mapped)) {
      registration <- .rc_fragment_registration_from_manifest(
        manifest,
        object_cells = colnames(obj)
      )
      obj <- .rc_register_signac_fragments(
        obj,
        fragment_files = registration$fragment_files,
        cells_by_fragment = registration$cell_maps,
        atac_assay = atac_assay,
        replace_existing = TRUE,
        require_complete = require_complete_fragments
      )
      fragment_registration <- "registered"
    } else {
      fragment_registration <- "not_registered_overlapping_fragment_files"
    }
    obj@misc$atac_count_source <- "recomputed_from_metacell_fragments"
    obj@misc$atac_fragment_recount$fragment_registration <-
      fragment_registration
  }
  obj
}
