.rc_fragment_objects_and_counts <- function(
    object, fragment_manifest, atac_assay = "ATAC",
    require_complete = TRUE, process_n = 2000L,
    create_fragment_fun = NULL, feature_matrix_fun = NULL) {
  if (!requireNamespace("Signac", quietly = TRUE) &&
      (is.null(create_fragment_fun) || is.null(feature_matrix_fun))) {
    stop("Package 'Signac' is required to recount ATAC peaks from fragments.",
         call. = FALSE)
  }
  if (!atac_assay %in% names(object@assays) ||
      !inherits(object[[atac_assay]], "ChromatinAssay")) {
    stop("ATAC assay `", atac_assay,
         "` must be a Signac ChromatinAssay.", call. = FALSE)
  }
  manifest <- .rc_validate_fragment_recount_manifest(
    fragment_manifest,
    object_cells = colnames(object),
    require_complete = require_complete
  )
  if (!nrow(manifest)) {
    return(list(counts = NULL, fragments = list(), manifest = manifest))
  }
  if (is.null(create_fragment_fun)) {
    create_fragment_fun <- getExportedValue("Signac", "CreateFragmentObject")
  }
  if (is.null(feature_matrix_fun)) {
    feature_matrix_fun <- getExportedValue("Signac", "FeatureMatrix")
  }
  features <- methods::slot(object[[atac_assay]], "ranges")
  feature_ids <- rownames(.rc_get_assay_counts_safe(object, atac_assay))
  cell_ids <- colnames(object)
  if (length(features) != length(feature_ids)) {
    stop("ATAC peak ranges and count-matrix rows have different lengths.",
         call. = FALSE)
  }
  names(features) <- feature_ids
  files <- unique(manifest$fragment_file)
  fragment_objects <- vector("list", length(files))
  count_matrices <- vector("list", length(files))
  for (i in seq_along(files)) {
    path <- files[[i]]
    rows <- manifest[manifest$fragment_file == path, , drop = FALSE]
    cell_map <- stats::setNames(
      as.character(rows$fragment_barcode),
      as.character(rows$object_cell)
    )
    fragment_objects[[i]] <- create_fragment_fun(
      path = path,
      cells = cell_map,
      validate.fragments = FALSE
    )
    counts_i <- feature_matrix_fun(
      fragments = list(fragment_objects[[i]]),
      features = features,
      keep_all_features = TRUE,
      cells = names(cell_map),
      process_n = as.integer(process_n),
      verbose = FALSE
    )
    count_matrices[[i]] <- .rc_align_peak_count_matrix(
      counts_i,
      feature_ids = feature_ids,
      cell_ids = cell_ids
    )
  }
  counts <- Reduce(`+`, count_matrices)
  counts <- .rc_as_sparse(counts)
  if (any(!is.finite(counts@x)) || any(counts@x < 0)) {
    stop("Fragment-derived ATAC peak counts contain invalid values.",
         call. = FALSE)
  }
  list(
    counts = counts,
    fragments = fragment_objects,
    manifest = manifest
  )
}
