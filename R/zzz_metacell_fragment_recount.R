if (!exists(".rc_aggregate_fragments_by_membership_without_recount", inherits = FALSE)) {
  .rc_aggregate_fragments_by_membership_without_recount <-
    rc_aggregate_fragments_by_membership
}
if (!exists(".rc_make_supercell2_metacells_without_recount", inherits = FALSE)) {
  .rc_make_supercell2_metacells_without_recount <- rc_make_supercell2_metacells
}
if (!exists(".rc_load_or_merge_metacell_objects_without_recount", inherits = FALSE)) {
  .rc_load_or_merge_metacell_objects_without_recount <-
    rc_load_or_merge_metacell_objects
}

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

rc_aggregate_fragments_by_membership <- function(
    fragment_files, membership, outdir, tmp_root = tempdir(),
    bgzip_path = "bgzip", tabix_path = "tabix", nb_cl = 1L) {
  manifest <- .rc_aggregate_fragments_by_membership_without_recount(
    fragment_files = fragment_files,
    membership = membership,
    outdir = outdir,
    tmp_root = tmp_root,
    bgzip_path = bgzip_path,
    tabix_path = tabix_path,
    nb_cl = nb_cl
  )
  .rc_expand_fragment_manifest(
    manifest,
    unique(as.character(membership$metacell_id))
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
  features <- Signac::granges(object[[atac_assay]])
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

.rc_recount_atac_from_fragment_manifest <- function(
    object, fragment_manifest, atac_assay = "ATAC",
    require_complete = TRUE, process_n = 2000L,
    create_fragment_fun = NULL, feature_matrix_fun = NULL) {
  answer <- .rc_fragment_objects_and_counts(
    object = object,
    fragment_manifest = fragment_manifest,
    atac_assay = atac_assay,
    require_complete = require_complete,
    process_n = process_n,
    create_fragment_fun = create_fragment_fun,
    feature_matrix_fun = feature_matrix_fun
  )
  if (is.null(answer$counts)) return(object)
  object <- SeuratObject::SetAssayData(
    object = object,
    assay = atac_assay,
    slot = "counts",
    new.data = answer$counts
  )
  object[[paste0("nCount_", atac_assay)]] <-
    as.numeric(Matrix::colSums(answer$counts))
  object[[paste0("nFeature_", atac_assay)]] <-
    as.numeric(Matrix::colSums(answer$counts > 0))
  mapped <- unlist(lapply(
    split(answer$manifest, answer$manifest$fragment_file),
    function(x) unique(as.character(x$object_cell))
  ), use.names = FALSE)
  if (!anyDuplicated(mapped)) {
    fragment_setter <- get("Fragments<-", envir = asNamespace("Signac"))
    object[[atac_assay]] <- fragment_setter(
      object[[atac_assay]],
      value = answer$fragments
    )
    registration <- "registered"
  } else {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
    registration <- "not_registered_overlapping_fragment_files"
  }
  object@misc$atac_count_source <- "recomputed_from_metacell_fragments"
  object@misc$atac_fragment_recount <- list(
    n_fragment_files = length(unique(answer$manifest$fragment_file)),
    n_peaks = nrow(answer$counts),
    n_metacells = ncol(answer$counts),
    fragment_registration = registration
  )
  object
}

rc_make_supercell2_metacells <- function(
    object,
    outdir,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    state_col = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    rna_reduction = "pca",
    atac_reduction = "lsi",
    rna_dims = 1:30,
    atac_dims = 2:30,
    gamma = 100,
    seed = 12345L,
    min_cells_per_stratum = 100,
    min_metacell_size = 20,
    min_metacells_per_stratum = 2L,
    label_col = NULL,
    fragment_files = NULL,
    bgzip_path = "bgzip",
    tabix_path = "tabix",
    fragment_nb_cl = 1L,
    save_metacell_object = TRUE,
    save_counts = TRUE,
    save_fragments = TRUE,
    require_fragment_aggregation = TRUE,
    fragment_aggregation_backend = c("regcompass", "supercell", "none"),
    overwrite = FALSE,
    BPPARAM = NULL,
    on_stratum_error = c("record", "stop")) {
  fragment_aggregation_backend <- match.arg(fragment_aggregation_backend)
  on_stratum_error <- match.arg(on_stratum_error)
  out <- .rc_make_supercell2_metacells_without_recount(
    object = object,
    outdir = outdir,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    state_col = state_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    rna_reduction = rna_reduction,
    atac_reduction = atac_reduction,
    rna_dims = rna_dims,
    atac_dims = atac_dims,
    gamma = gamma,
    seed = seed,
    min_cells_per_stratum = min_cells_per_stratum,
    min_metacell_size = min_metacell_size,
    min_metacells_per_stratum = min_metacells_per_stratum,
    label_col = label_col,
    fragment_files = fragment_files,
    bgzip_path = bgzip_path,
    tabix_path = tabix_path,
    fragment_nb_cl = fragment_nb_cl,
    save_metacell_object = save_metacell_object,
    save_counts = save_counts,
    save_fragments = save_fragments,
    require_fragment_aggregation = require_fragment_aggregation,
    fragment_aggregation_backend = fragment_aggregation_backend,
    overwrite = overwrite,
    BPPARAM = BPPARAM,
    on_stratum_error = on_stratum_error
  )
  fragment_enabled <- !identical(fragment_files, FALSE) &&
    isTRUE(save_fragments) &&
    !identical(fragment_aggregation_backend, "none")
  if (!fragment_enabled) {
    out$atac_count_source <- "aggregated_object_peak_counts"
    return(out)
  }
  if (!is.data.frame(out$fragment_manifest) ||
      !nrow(out$fragment_manifest)) {
    if (isTRUE(require_fragment_aggregation)) {
      stop("Fragment aggregation completed without a usable fragment manifest.",
           call. = FALSE)
    }
    warning(
      "No metacell fragment manifest was produced; retaining SuperCell ATAC ",
      "peak counts because fragment aggregation was not required.",
      call. = FALSE
    )
    out$atac_count_source <- "aggregated_object_peak_counts"
    return(out)
  }
  object_files <- as.character(out$metacell_objects)
  for (object_file in object_files) {
    mc <- readRDS(object_file)
    stratum_dir <- dirname(object_file)
    manifest_i <- out$fragment_manifest
    if ("stratum_dir" %in% colnames(manifest_i)) {
      manifest_i <- manifest_i[
        normalizePath(manifest_i$stratum_dir, mustWork = FALSE) ==
          normalizePath(stratum_dir, mustWork = FALSE),
        , drop = FALSE
      ]
    } else {
      manifest_i <- manifest_i[
        as.character(manifest_i$object_cell) %in% colnames(mc),
        , drop = FALSE
      ]
    }
    manifest_i <- .rc_expand_fragment_manifest(manifest_i, colnames(mc))
    mc <- .rc_recount_atac_from_fragment_manifest(
      object = mc,
      fragment_manifest = manifest_i,
      atac_assay = atac_assay,
      require_complete = TRUE
    )
    saveRDS(mc, object_file)
    saveRDS(
      .rc_as_sparse(.rc_get_assay_counts_safe(mc, atac_assay)),
      file.path(stratum_dir, "atac_counts.rds")
    )
    .rc_write_tsv_gz(
      manifest_i,
      file.path(stratum_dir, "fragments", "fragment_manifest.tsv.gz")
    )
  }
  refreshed <- rc_import_supercell2_metacells(
    unique(dirname(object_files)),
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    require_fragments = TRUE
  )
  refreshed$stratum_status <- out$stratum_status
  refreshed$atac_count_source <- "recomputed_from_metacell_fragments"
  refreshed
}

rc_load_or_merge_metacell_objects <- function(
    metacell_objects, fragment_manifest = NULL, metacell_meta = NULL,
    fragment_files = NULL, rna_assay = "RNA", atac_assay = "ATAC",
    require_complete_fragments = TRUE) {
  if (is.null(fragment_manifest) || !is.data.frame(fragment_manifest) ||
      !nrow(fragment_manifest)) {
    return(.rc_load_or_merge_metacell_objects_without_recount(
      metacell_objects = metacell_objects,
      fragment_manifest = fragment_manifest,
      metacell_meta = metacell_meta,
      fragment_files = fragment_files,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      require_complete_fragments = require_complete_fragments
    ))
  }
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
  precounted <- all(vapply(
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
  if (!precounted) {
    obj <- .rc_recount_atac_from_fragment_manifest(
      object = obj,
      fragment_manifest = fragment_manifest,
      atac_assay = atac_assay,
      require_complete = require_complete_fragments
    )
  } else {
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
