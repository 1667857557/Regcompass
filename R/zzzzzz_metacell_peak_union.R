.rc_cbind_sparse_feature_union <- function(matrices) {
  matrices <- matrices[
    vapply(
      matrices,
      function(x) !is.null(dim(x)) && ncol(x) > 0L,
      logical(1)
    )
  ]
  if (!length(matrices)) return(NULL)
  feature_ids <- unique(unlist(lapply(matrices, rownames), use.names = FALSE))
  cell_ids <- unlist(lapply(matrices, colnames), use.names = FALSE)
  if (anyNA(feature_ids) || any(!nzchar(feature_ids)) ||
      anyNA(cell_ids) || any(!nzchar(cell_ids))) {
    stop("Metacell matrices require non-empty feature and cell names.",
         call. = FALSE)
  }
  if (anyDuplicated(cell_ids)) {
    duplicated_ids <- unique(cell_ids[duplicated(cell_ids)])
    stop("Duplicated metacell IDs across matrices: ",
         paste(utils::head(duplicated_ids, 10L), collapse = ", "),
         call. = FALSE)
  }
  aligned <- lapply(matrices, function(x) {
    x <- .rc_as_sparse(x)
    row_map <- Matrix::sparseMatrix(
      i = match(rownames(x), feature_ids),
      j = seq_len(nrow(x)),
      x = 1,
      dims = c(length(feature_ids), nrow(x)),
      dimnames = list(feature_ids, rownames(x))
    )
    row_map %*% x
  })
  out <- do.call(cbind, aligned)
  rownames(out) <- feature_ids
  colnames(out) <- cell_ids
  .rc_as_sparse(out)
}

rc_import_supercell2_metacells <- function(
    metacell_dirs,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    require_fragments = FALSE) {
  metacell_dirs <- metacell_dirs[dir.exists(metacell_dirs)]
  if (!length(metacell_dirs)) {
    stop("No valid metacell directories supplied.", call. = FALSE)
  }
  read_tsv <- function(path) {
    utils::read.delim(
      gzfile(path),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  metas <- memberships <- fragment_manifests <- list()
  rnas <- atacs <- list()
  objects <- fragments <- character()
  for (directory in metacell_dirs) {
    meta_file <- file.path(directory, "metacell_metadata.tsv.gz")
    if (!file.exists(meta_file)) next
    metas[[directory]] <- read_tsv(meta_file)
    membership_file <- file.path(directory, "membership.tsv.gz")
    if (file.exists(membership_file)) {
      memberships[[directory]] <- read_tsv(membership_file)
    }
    rna_file <- file.path(directory, "rna_counts.rds")
    atac_file <- file.path(directory, "atac_counts.rds")
    object_file <- file.path(directory, "metacell_object.rds")
    if (file.exists(rna_file)) rnas[[directory]] <- readRDS(rna_file)
    if (file.exists(atac_file)) atacs[[directory]] <- readRDS(atac_file)
    if (file.exists(object_file)) objects <- c(objects, object_file)
    manifest_file <- file.path(
      directory, "fragments", "fragment_manifest.tsv.gz"
    )
    if (file.exists(manifest_file)) {
      manifest <- read_tsv(manifest_file)
      manifest$stratum_dir <- directory
      fragment_manifests[[directory]] <- manifest
    } else {
      discovered <- Sys.glob(file.path(directory, "fragments", "*.tsv.gz"))
      discovered <- setdiff(discovered, manifest_file)
      if (length(discovered)) fragments <- c(fragments, discovered)
    }
  }
  if (!length(metas)) {
    stop("No metacell metadata files were found.", call. = FALSE)
  }
  metacell_meta <- do.call(rbind, metas)
  rownames(metacell_meta) <- NULL
  metacell_meta$metacell_id <- as.character(metacell_meta$metacell_id)
  if (anyDuplicated(metacell_meta$metacell_id)) {
    duplicated_ids <- unique(
      metacell_meta$metacell_id[duplicated(metacell_meta$metacell_id)]
    )
    stop("Duplicated metacell IDs across strata: ",
         paste(utils::head(duplicated_ids, 10L), collapse = ", "),
         call. = FALSE)
  }
  membership <- if (length(memberships)) {
    do.call(rbind, memberships)
  } else {
    data.frame()
  }
  if (!length(rnas)) {
    stop("No rna_counts.rds files were found in metacell directories.",
         call. = FALSE)
  }
  rna_counts <- .rc_cbind_sparse_feature_union(rnas)
  atac_counts <- .rc_cbind_sparse_feature_union(atacs)
  colnames(rna_counts) <- as.character(colnames(rna_counts))
  if (!is.null(atac_counts)) {
    colnames(atac_counts) <- as.character(colnames(atac_counts))
    if (!setequal(colnames(rna_counts), colnames(atac_counts))) {
      stop("RNA and ATAC metacell IDs differ after import.", call. = FALSE)
    }
    atac_counts <- atac_counts[, colnames(rna_counts), drop = FALSE]
  }
  metacell_meta <- metacell_meta[
    match(colnames(rna_counts), metacell_meta$metacell_id),
    , drop = FALSE
  ]
  if (anyNA(metacell_meta$metacell_id)) {
    stop("Metacell metadata are incomplete.", call. = FALSE)
  }
  rc_validate_metacell_inputs(
    rna_counts,
    metacell_meta,
    atac_metacell_counts = atac_counts,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  fragment_manifest <- if (length(fragment_manifests)) {
    do.call(rbind, fragment_manifests)
  } else {
    data.frame()
  }
  if (nrow(fragment_manifest)) {
    fragments <- unique(as.character(fragment_manifest$fragment_file))
  } else if (length(fragments)) {
    warning(
      "Legacy fragment discovery by glob was used because no fragment ",
      "manifest was found.",
      call. = FALSE
    )
  }
  if (isTRUE(require_fragments)) {
    missing_idx <- fragments[
      !file.exists(paste0(fragments, ".tbi")) &
        !file.exists(paste0(fragments, ".csi"))
    ]
    if (!length(fragments) || length(missing_idx)) {
      stop("Metacell fragment files or indexes are missing.", call. = FALSE)
    }
  }
  list(
    schema_version = "regcompass_metacell_v1.1_peak_union",
    metacell_meta = metacell_meta,
    membership = membership,
    rna_counts = rna_counts,
    atac_counts = atac_counts,
    metacell_objects = objects,
    fragment_manifest = fragment_manifest,
    fragment_files = fragments,
    diagnostics = data.frame(
      n_metacells = ncol(rna_counts),
      n_membership_rows = nrow(membership),
      n_atac_peaks = if (is.null(atac_counts)) 0L else nrow(atac_counts)
    )
  )
}
