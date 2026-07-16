#' Validate metacell-level RegCompass inputs
rc_validate_metacell_inputs <- function(rna_metacell_counts,
                                        metacell_meta,
                                        atac_metacell_counts = NULL,
                                        metacell_id_col = "metacell_id",
                                        sample_col = "sample_id",
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

  strict_columns <- intersect(
    c("sample_id", "condition", "cell_type"),
    colnames(x)
  )
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

.rc_safe_path_component <- function(x) gsub("[^A-Za-z0-9_.=-]+", "_", as.character(x))

.rc_write_tsv_gz <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(file, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(file)
}

.rc_get_assay_counts_safe <- function(object, assay) {
  .rc_get_assay_counts(object, assay)
}

.rc_extract_supercell_membership <- function(mc_object, original_cells, metacell_ids) {
  make_df <- function(cell_id, metacell_id) {
    out <- data.frame(cell_id = as.character(cell_id), metacell_id = as.character(metacell_id), stringsAsFactors = FALSE)
    out <- out[!is.na(out$cell_id) & nzchar(out$cell_id) & !is.na(out$metacell_id) & nzchar(out$metacell_id), , drop = FALSE]
    out <- out[!duplicated(out$cell_id), , drop = FALSE]
    rownames(out) <- NULL
    out
  }
  normalize_ids <- function(x) {
    x <- as.character(x)
    if (length(x) > 0L && all(x %in% metacell_ids)) return(x)
    suppressWarnings(ix <- as.integer(x))
    if (length(ix) > 0L && all(is.finite(ix)) && all(ix >= 1L) && all(ix <= length(metacell_ids))) return(metacell_ids[ix])
    x
  }
  misc_table <- tryCatch(mc_object@misc$membership_table, error = function(e) NULL)
  if (is.data.frame(misc_table) && all(c("cell_id", "metacell_id") %in% colnames(misc_table))) {
    out <- misc_table[, c("cell_id", "metacell_id"), drop = FALSE]
    out$cell_id <- as.character(out$cell_id)
    out$metacell_id <- as.character(out$metacell_id)
    if (anyDuplicated(out$cell_id)) stop("Duplicated cell IDs in SuperCell membership_table.", call. = FALSE)
    if (!setequal(unique(out$metacell_id), as.character(metacell_ids))) {
      stop("SuperCell membership_table metacell IDs do not match metacell object colnames.", call. = FALSE)
    }
    out <- out[out$cell_id %in% original_cells, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  meta <- mc_object@meta.data
  candidates <- c("cell_membership", "membership", "cells", "cell_ids", "single_cell_ids", "SC")
  for (nm in intersect(candidates, colnames(meta))) {
    vals <- meta[[nm]]
    if (is.list(vals)) {
      return(do.call(rbind, lapply(seq_along(vals), function(i) make_df(vals[[i]], metacell_ids[[i]]))))
    }
  }
  misc_mem <- tryCatch(mc_object@misc$membership, error = function(e) NULL)
  if (!is.null(misc_mem)) {
    if (is.data.frame(misc_mem)) {
      if (!"cell_id" %in% colnames(misc_mem) && "cell" %in% colnames(misc_mem)) misc_mem$cell_id <- misc_mem$cell
      if (!"metacell_id" %in% colnames(misc_mem) && "metacell" %in% colnames(misc_mem)) misc_mem$metacell_id <- misc_mem$metacell
      if (all(c("cell_id", "metacell_id") %in% colnames(misc_mem))) return(make_df(misc_mem$cell_id, normalize_ids(misc_mem$metacell_id)))
    }
    if (is.atomic(misc_mem) && length(misc_mem) == length(original_cells)) {
      cell_ids <- names(misc_mem)
      if (is.null(cell_ids) || any(!nzchar(cell_ids))) cell_ids <- original_cells
      return(make_df(cell_ids, normalize_ids(misc_mem)))
    }
  }
  wt <- tryCatch(mc_object@misc$walktrap_clusters, error = function(e) NULL)
  if (!is.null(wt) && is.atomic(wt) && length(wt) == length(original_cells)) {
    cell_ids <- names(wt)
    if (is.null(cell_ids) || any(!nzchar(cell_ids))) cell_ids <- original_cells
    return(make_df(cell_ids, normalize_ids(wt)))
  }
  hierarchy <- tryCatch(mc_object@misc$metacells_hierarchy, error = function(e) NULL)
  if (!is.null(hierarchy) && requireNamespace("igraph", quietly = TRUE)) {
    mem <- tryCatch(igraph::membership(hierarchy), error = function(e) NULL)
    if (!is.null(mem) && length(mem) == length(original_cells)) {
      cell_ids <- names(mem)
      if (is.null(cell_ids) || any(!nzchar(cell_ids))) cell_ids <- original_cells
      return(make_df(cell_ids, normalize_ids(mem)))
    }
  }
  attr_map <- attr(mc_object, "membership")
  if (!is.null(attr_map) && is.data.frame(attr_map) && all(c("cell_id", "metacell_id") %in% colnames(attr_map))) return(make_df(attr_map$cell_id, normalize_ids(attr_map$metacell_id)))
  data.frame(cell_id = character(0), metacell_id = character(0), stringsAsFactors = FALSE)
}

#' Aggregate ATAC fragments by metacell membership
rc_aggregate_fragments_by_membership <- function(fragment_files, membership, outdir, tmp_root = tempdir(), bgzip_path = "bgzip", tabix_path = "tabix", nb_cl = 1L) {
  if (!is.data.frame(membership) || !all(c("cell_id", "metacell_id") %in% colnames(membership))) stop("`membership` must be a data.frame containing cell_id and metacell_id.", call. = FALSE)
  if (anyDuplicated(membership$cell_id)) stop("Duplicated cell IDs in membership.", call. = FALSE)
  if (!requireNamespace("SuperCell", quietly = TRUE)) stop("Package 'SuperCell' is required for fragment aggregation.", call. = FALSE)
  agg <- getExportedValue("SuperCell", "AggregateFragmentFile")
  required_formals <- c("input_file", "membership", "output_name", "output_path")
  missing_formals <- setdiff(required_formals, names(formals(agg)))
  if (length(missing_formals)) {
    stop("Installed SuperCell::AggregateFragmentFile() has an incompatible API. Missing: ", paste(missing_formals, collapse = ", "), call. = FALSE)
  }
  files <- unique(as.character(unlist(fragment_files, use.names = FALSE)))
  files <- files[!is.na(files) & nzchar(files)]
  if (!length(files)) stop("No fragment files supplied.", call. = FALSE)
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) stop("Fragment files not found: ", paste(utils::head(missing_files, 10L), collapse = ", "), call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
  map <- stats::setNames(as.character(membership$metacell_id), as.character(membership$cell_id))
  rows <- lapply(seq_along(files), function(i) {
    input_file <- files[[i]]
    output_file <- agg(
      input_file = input_file,
      membership = map,
      output_name = paste0("MC_", sprintf("%03d", i), "_", basename(input_file)),
      output_path = outdir,
      tmp_path = file.path(tmp_root, paste0("fragment_", sprintf("%03d", i))),
      nb_cl = max(1L, as.integer(nb_cl)),
      bgzip_path = bgzip_path,
      tabix_path = tabix_path,
      returnOutputFileName = TRUE
    )
    data.frame(input_file = normalizePath(input_file, mustWork = FALSE),
               fragment_file = normalizePath(output_file, mustWork = FALSE),
               index_file = normalizePath(paste0(output_file, ".tbi"), mustWork = FALSE),
               status = "ok", stringsAsFactors = FALSE)
  })
  manifest <- do.call(rbind, rows)
  if (any(!file.exists(manifest$fragment_file)) || any(!file.exists(manifest$index_file))) {
    stop("Fragment aggregation did not produce all required fragment files and indexes.", call. = FALSE)
  }
  manifest
}

.rc_aggregate_fragments_by_stratum <- function(fragment_files,
                                               membership,
                                               outdir,
                                               stratum_cols,
                                               sample_col = "sample_id",
                                               atac_assay = "ATAC",
                                               tmp_root = tempdir(),
                                               bgzip_path = "bgzip",
                                               tabix_path = "tabix",
                                               nb_cl = 1L) {
  if (!is.data.frame(membership) || !all(c("cell_id", "metacell_id") %in% colnames(membership))) {
    stop("`membership` must be a data.frame containing cell_id and metacell_id.", call. = FALSE)
  }
  stratum_cols <- stratum_cols[!is.null(stratum_cols) & !is.na(stratum_cols) & nzchar(stratum_cols)]
  missing <- setdiff(c(stratum_cols, sample_col), colnames(membership))
  if (length(missing)) stop("Membership is missing stratum columns required for fragment aggregation: ", paste(missing, collapse = ", "), call. = FALSE)
  membership$.rc_fragment_stratum <- rc_make_stratum_id(membership, stratum_cols)
  input_manifest <- .rc_normalize_fragment_manifest(fragment_files, sample_ids = membership[[sample_col]], atac_assay = atac_assay)
  if (!nrow(input_manifest)) stop("No fragment files supplied for stratum-wise fragment aggregation.", call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pieces <- lapply(split(membership, membership$.rc_fragment_stratum), function(mem_i) {
    st <- mem_i$.rc_fragment_stratum[[1L]]
    vals <- mem_i[1, stratum_cols, drop = FALSE]
    sample_value <- as.character(mem_i[[sample_col]][[1L]])
    files_i <- input_manifest$fragment_file[input_manifest$sample_id == sample_value & input_manifest$assay == atac_assay]
    files_i <- unique(as.character(files_i[!is.na(files_i) & nzchar(files_i)]))
    if (!length(files_i)) stop("No fragment file was mapped to stratum `", st, "` sample `", sample_value, "`.", call. = FALSE)
    stratum_dir <- file.path(outdir, .rc_safe_path_component(st))
    manifest_i <- rc_aggregate_fragments_by_membership(
      fragment_files = files_i,
      membership = mem_i[, c("cell_id", "metacell_id"), drop = FALSE],
      outdir = stratum_dir,
      tmp_root = file.path(tmp_root, .rc_safe_path_component(st)),
      bgzip_path = bgzip_path,
      tabix_path = tabix_path,
      nb_cl = nb_cl
    )
    ids_i <- unique(as.character(mem_i$metacell_id))
    manifest_i <- do.call(rbind, lapply(seq_len(nrow(manifest_i)), function(i) {
      row <- manifest_i[i, , drop = FALSE]
      map <- data.frame(object_cell = ids_i, fragment_barcode = ids_i, stringsAsFactors = FALSE)
      cbind(row, vals[rep(1L, nrow(map)), , drop = FALSE], link_stratum = st, map)
    }))
    rownames(manifest_i) <- NULL
    manifest_i
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

.rc_normalize_fragment_manifest <- function(fragment_files, sample_ids, atac_assay = "ATAC") {
  if (is.null(fragment_files)) return(data.frame(sample_id = character(), assay = character(), fragment_file = character(), stringsAsFactors = FALSE))
  sample_ids <- unique(as.character(sample_ids))
  if (is.data.frame(fragment_files)) {
    required <- c("sample_id", "assay", "fragment_file")
    missing <- setdiff(required, colnames(fragment_files))
    if (length(missing)) stop("Fragment manifest missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    out <- fragment_files[, required, drop = FALSE]
    out$sample_id <- as.character(out$sample_id)
    out$assay <- as.character(out$assay)
    out$fragment_file <- as.character(out$fragment_file)
  } else {
    files <- as.character(unlist(fragment_files, use.names = FALSE))
    files <- files[!is.na(files) & nzchar(files)]
    if (length(sample_ids) != 1L && length(files) != 1L) stop("Multi-sample input requires a fragment manifest containing sample_id, assay and fragment_file, unless one shared fragment file is supplied for all samples.", call. = FALSE)
    out <- data.frame(sample_id = if (length(sample_ids) == 1L) sample_ids[[1L]] else sample_ids, assay = atac_assay, fragment_file = if (length(files) == 1L) files[[1L]] else files, stringsAsFactors = FALSE)
  }
  out <- out[!is.na(out$fragment_file) & nzchar(out$fragment_file), , drop = FALSE]
  if (any(!file.exists(out$fragment_file))) stop("One or more fragment files in the manifest do not exist.", call. = FALSE)
  rownames(out) <- NULL
  out
}


.rc_fragment_path_from_object <- function(fragment) {
  path <- NULL
  if (is.list(fragment) && !is.null(fragment$path)) path <- fragment$path
  if (is.null(path)) path <- attr(fragment, "path", exact = TRUE)
  if (is.null(path) && methods::is(fragment, "Fragment") && methods::isS4(fragment) && "path" %in% methods::slotNames(fragment)) {
    path <- methods::slot(fragment, "path")
  }
  path <- as.character(path %||% character(0))
  path[!is.na(path) & nzchar(path)]
}

.rc_fragment_files_from_atac <- function(object, atac_assay = "ATAC") {
  if (!requireNamespace("Signac", quietly = TRUE)) return(NULL)
  assay <- tryCatch(object[[atac_assay]], error = function(e) NULL)
  if (is.null(assay)) return(NULL)
  fragments <- tryCatch(Signac::Fragments(assay), error = function(e) list())
  paths <- unlist(lapply(fragments, .rc_fragment_path_from_object), use.names = FALSE)
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  if (length(paths) == 0L) return(NULL)
  stats::setNames(list(paths), atac_assay)
}

.rc_normalize_fragment_files <- function(fragment_files, atac_assay = "ATAC") {
  if (is.null(fragment_files)) return(NULL)
  if (is.character(fragment_files)) {
    if (length(fragment_files) == 1L && (is.null(names(fragment_files)) || !nzchar(names(fragment_files)[1L]))) {
      return(stats::setNames(list(fragment_files[[1L]]), atac_assay))
    }
    if (is.null(names(fragment_files)) || any(!nzchar(names(fragment_files)))) {
      stop("`fragment_files` must be a named character vector/list when more than one file is supplied.", call. = FALSE)
    }
    return(as.list(fragment_files))
  }
  if (is.list(fragment_files)) {
    if (length(fragment_files) == 1L && (is.null(names(fragment_files)) || !nzchar(names(fragment_files)[1L]))) {
      names(fragment_files) <- atac_assay
    }
    if (is.null(names(fragment_files)) || any(!nzchar(names(fragment_files)))) {
      stop("`fragment_files` list must be named by chromatin assay, e.g. list(ATAC = path).", call. = FALSE)
    }
    return(fragment_files)
  }
  stop("`fragment_files` must be NULL, a character path/vector, or a named list.", call. = FALSE)
}


.rc_require_supercell2 <- function() {
  if (!requireNamespace("SuperCell", quietly = TRUE)) {
    stop("Package 'SuperCell' is required for rc_make_supercell2_metacells(). Install the SuperCell2 branch with `remotes::install_github(\"1667857557/SuperCell_Seurat_V4@supercell-2.0\")` (or the upstream mirror `GfellerLab/SuperCell@supercell-2.0`).", call. = FALSE)
  }
  if (!exists("SCimplify_for_Seurat", envir = asNamespace("SuperCell"), inherits = FALSE)) {
    stop("Installed package 'SuperCell' does not export SCimplify_for_Seurat(); install the SuperCell2 branch before running metacells.", call. = FALSE)
  }
  version <- tryCatch(utils::packageVersion("SuperCell"), error = function(e) NULL)
  if (!is.null(version) && version < "2.0") {
    stop("RegCompass requires SuperCell2 (>= 2.0) for multimodal Seurat metacells. Reinstall with `remotes::install_github(\"1667857557/SuperCell_Seurat_V4@supercell-2.0\")`.", call. = FALSE)
  }
  invisible(TRUE)
}

.rc_validate_supercell2_inputs <- function(object, assays, reductions) {
  missing_assays <- setdiff(assays, names(object@assays))
  if (length(missing_assays) > 0L) stop("Seurat object is missing assay(s) required by SuperCell2: ", paste(missing_assays, collapse = ", "), call. = FALSE)
  reductions <- unlist(reductions, use.names = FALSE)
  reductions <- reductions[!is.na(reductions) & nzchar(reductions)]
  missing_reductions <- setdiff(reductions, names(object@reductions))
  if (length(missing_reductions) > 0L) {
    stop("Seurat object is missing reduction(s) required by SuperCell2: ", paste(missing_reductions, collapse = ", "), ". Run the corresponding Seurat/Signac dimensional reduction first or pass existing `rna_reduction`/`atac_reduction` names.", call. = FALSE)
  }
  invisible(TRUE)
}

.rc_supercell2_scimplify_for_seurat <- function(args) {
  .rc_require_supercell2()
  .rc_with_seurat4_filterobjects(do.call(getExportedValue("SuperCell", "SCimplify_for_Seurat"), args))
}

.rc_assert_shell_safe_paths <- function(...) {
  paths <- unlist(list(...), use.names = FALSE)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  bad <- paths[grepl("\\s", paths)]
  if (length(bad)) {
    warning("Some paths contain whitespace. RegCompassR sanitizes stratum directories, but external fragment tools may still fail if shell paths are not quoted: ", paste(unique(bad), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Build sample-aware SuperCell2.0 RNA+ATAC metacells and save outputs
rc_make_supercell2_metacells <- function(object,
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
                                         adaptive_gamma = FALSE,
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
  .rc_require_supercell2()
  if (!inherits(object, "Seurat")) stop("`object` must inherit from class 'Seurat'.", call. = FALSE)
  if (is.null(fragment_files)) fragment_files <- .rc_fragment_files_from_atac(object, atac_assay = atac_assay)
  if (isTRUE(require_fragment_aggregation)) {
    if (!isTRUE(save_fragments)) stop("Formal multiome workflow requires `save_fragments = TRUE`.", call. = FALSE)
    if (is.null(fragment_files)) stop("Formal multiome workflow requires `fragment_files` for metacell fragment aggregation, or a fragment file registered on the ATAC assay.", call. = FALSE)
  }
  .rc_validate_supercell2_inputs(object, assays = c(rna_assay, atac_assay), reductions = c(rna_reduction, atac_reduction))
  meta <- object@meta.data
  required <- c(sample_col, condition_col, celltype_col, state_col, label_col)
  required <- required[!is.null(required) & !is.na(required) & nzchar(required)]
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0L) stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  group_cols <- .rc_strict_stratum_cols(sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col)
  meta <- .rc_add_stratum_id(meta, group_cols)
  fragment_manifest <- .rc_normalize_fragment_manifest(fragment_files, sample_ids = meta[[sample_col]], atac_assay = atac_assay)
  .rc_assert_shell_safe_paths(outdir, fragment_manifest$fragment_file)
  meta$cell_id <- rownames(meta)
  keys <- meta$.rc_stratum_id
  groups <- split(meta$cell_id, keys)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  run_one <- function(key) {
    cells <- groups[[key]]
    one_meta <- meta[match(cells, meta$cell_id), , drop = FALSE]
    vals <- one_meta[1, group_cols, drop = FALSE]
    dir_name <- paste(paste0(group_cols, "=", vapply(vals, .rc_safe_path_component, character(1))), collapse = "__")
    stratum_dir <- file.path(outdir, dir_name)
    if (dir.exists(stratum_dir) && !overwrite && file.exists(file.path(stratum_dir, "metacell_metadata.tsv.gz"))) return(stratum_dir)
    dir.create(file.path(stratum_dir, "fragments"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(stratum_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
    min_required_cells <- as.integer(min_cells_per_stratum)
    gamma_i <- as.integer(gamma)
    if (length(cells) < min_required_cells) {
      reason <- "stratum_below_min_cells_per_stratum"
      diag <- data.frame(group_id = key, n_cells = length(cells), skipped = TRUE, skip_reason = reason, gamma = gamma, min_required_cells = min_required_cells, stringsAsFactors = FALSE)
      .rc_write_tsv_gz(diag, file.path(stratum_dir, "qc", "metacell_qc.tsv.gz"))
      return(stratum_dir)
    }
    if (isTRUE(adaptive_gamma)) {
      gamma_i <- min(as.integer(gamma), floor(length(cells) / as.integer(min_metacells_per_stratum)))
      gamma_i <- max(gamma_i, as.integer(min_metacell_size))
    }
    prefix <- paste(vapply(vals, .rc_safe_path_component, character(1)), collapse = "_")
    prefix_mc <- paste0(prefix, "_MC_")
    sample_value <- as.character(vals[[sample_col]][[1L]])
    fragment_files_i <- fragment_manifest$fragment_file[fragment_manifest$sample_id == sample_value & fragment_manifest$assay == atac_assay]
    if (isTRUE(require_fragment_aggregation) && save_fragments && !length(fragment_files_i)) {
      stop("No fragment file was mapped to sample: ", sample_value, call. = FALSE)
    }
    seu_sub <- subset(object, cells = cells)
    seed_i <- as.integer(seed) + match(key, names(groups)) - 1L
    args <- list(seurat = seu_sub, assay = c(rna_assay, atac_assay), reduction = list(rna_reduction, atac_reduction), dims = list(rna_dims, atac_dims), gamma = gamma_i, return.seurat = TRUE, prefixMC = prefix_mc, seed = seed_i)
    if (!is.null(label_col)) args$label <- label_col
    if (identical(fragment_aggregation_backend, "supercell") && save_fragments && length(fragment_files_i)) {
      args$fragmentFiles <- stats::setNames(list(fragment_files_i), atac_assay)
      args$outputDirMcFragment <- file.path(stratum_dir, "fragments")
      args$bgzip_path <- bgzip_path
      args$tabix_path <- tabix_path
      args$nb_cl <- max(1L, as.integer(fragment_nb_cl))
    }
    mc <- tryCatch(
      .rc_supercell2_scimplify_for_seurat(args),
      error = function(e) {
        if (identical(fragment_aggregation_backend, "supercell")) {
          if (isTRUE(require_fragment_aggregation)) stop("SuperCell fragment aggregation failed before metacell output was returned: ", conditionMessage(e), call. = FALSE)
          warning("SuperCell fragment aggregation failed; retrying metacell construction without fragment arguments because `require_fragment_aggregation = FALSE`.", call. = FALSE)
          args2 <- args[setdiff(names(args), c("fragmentFiles", "outputDirMcFragment", "bgzip_path", "tabix_path", "nb_cl"))]
          return(.rc_supercell2_scimplify_for_seurat(args2))
        }
        stop("SuperCell2 metacell construction failed: ", conditionMessage(e), call. = FALSE)
      }
    )
    mc_ids <- colnames(.rc_get_assay_counts_safe(mc, rna_assay))
    mc_ids <- as.character(mc_ids)
    if (anyDuplicated(mc_ids)) stop("Duplicated metacell IDs within stratum.", call. = FALSE)
    display_ids <- paste0(prefix, "_MC", sprintf(paste0("%0", max(3, nchar(length(mc_ids))), "d"), seq_along(mc_ids)))
    membership <- .rc_extract_supercell_membership(mc, cells, mc_ids)
    if (nrow(membership) == 0L) stop("Could not infer single-cell membership from SuperCell output for stratum ", key, call. = FALSE)
    fragment_manifest_i <- NULL
    if (identical(fragment_aggregation_backend, "supercell") && save_fragments) {
      fragment_manifest_i <- tryCatch(mc@misc$fragment_manifest, error = function(e) NULL)
      if (is.null(fragment_manifest_i) || !is.data.frame(fragment_manifest_i) || !nrow(fragment_manifest_i)) {
        produced <- Sys.glob(file.path(stratum_dir, "fragments", "*.tsv.gz"))
        produced <- produced[basename(produced) != "fragment_manifest.tsv.gz"]
        if (length(produced)) {
          fragment_manifest_i <- data.frame(
            input_file = NA_character_,
            fragment_file = normalizePath(produced, mustWork = FALSE),
            index_file = normalizePath(paste0(produced, ".tbi"), mustWork = FALSE),
            status = "ok",
            stringsAsFactors = FALSE
          )
        }
      }
      if (!is.null(fragment_manifest_i) && is.data.frame(fragment_manifest_i) && nrow(fragment_manifest_i)) {
        if (!"input_file" %in% colnames(fragment_manifest_i)) fragment_manifest_i$input_file <- NA_character_
        if (!"fragment_file" %in% colnames(fragment_manifest_i) && "output_file" %in% colnames(fragment_manifest_i)) fragment_manifest_i$fragment_file <- fragment_manifest_i$output_file
        if (!"index_file" %in% colnames(fragment_manifest_i)) fragment_manifest_i$index_file <- paste0(fragment_manifest_i$fragment_file, ".tbi")
        if (!"status" %in% colnames(fragment_manifest_i)) fragment_manifest_i$status <- "ok"
        fragment_manifest_i$stratum_id <- key
        fragment_manifest_i$sample_id <- sample_value
        fragment_manifest_i$assay <- atac_assay
        fragment_manifest_i$metacell_prefix <- prefix_mc
        .rc_write_tsv_gz(fragment_manifest_i, file.path(stratum_dir, "fragments", "fragment_manifest.tsv.gz"))
      }
    }
    if (identical(fragment_aggregation_backend, "regcompass") && save_fragments && length(fragment_files_i)) {
      fragment_manifest_i <- tryCatch(
        rc_aggregate_fragments_by_membership(fragment_files = fragment_files_i, membership = membership, outdir = file.path(stratum_dir, "fragments"), bgzip_path = bgzip_path, tabix_path = tabix_path, nb_cl = max(1L, as.integer(fragment_nb_cl))),
        error = function(e) {
          if (isTRUE(require_fragment_aggregation)) stop("Metacell fragment aggregation failed: ", conditionMessage(e), call. = FALSE)
          warning("Fragment aggregation failed; continuing only because `require_fragment_aggregation = FALSE`: ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (!is.null(fragment_manifest_i)) {
        fragment_manifest_i$stratum_id <- key
        fragment_manifest_i$sample_id <- sample_value
        fragment_manifest_i$assay <- atac_assay
        fragment_manifest_i$metacell_prefix <- prefix_mc
        .rc_write_tsv_gz(fragment_manifest_i, file.path(stratum_dir, "fragments", "fragment_manifest.tsv.gz"))
      }
    }
    for (col in group_cols) membership[[col]] <- vals[[col]][[1]]
    mc_meta <- rc_build_metacell_metadata(membership)
    if (nrow(mc_meta) == 0L) mc_meta <- data.frame(metacell_id = mc_ids, n_cells = NA_integer_, stringsAsFactors = FALSE)
    for (col in group_cols) if (!col %in% colnames(mc_meta)) mc_meta[[col]] <- vals[[col]][[1]]
    mc_meta$stratum_id <- key
    mc_meta$metacell_display_id <- display_ids[match(as.character(mc_meta$metacell_id), mc_ids)]
    mc_meta$low_power_metacell <- !is.na(mc_meta$n_cells) & mc_meta$n_cells < min_metacell_size
    mc_meta$effective_gamma <- gamma_i
    rna_counts <- .rc_as_sparse(.rc_get_assay_counts_safe(mc, rna_assay))
    atac_counts <- .rc_as_sparse(.rc_get_assay_counts_safe(mc, atac_assay))
    if (save_metacell_object) saveRDS(mc, file.path(stratum_dir, "metacell_object.rds"))
    .rc_write_tsv_gz(membership, file.path(stratum_dir, "membership.tsv.gz"))
    .rc_write_tsv_gz(mc_meta, file.path(stratum_dir, "metacell_metadata.tsv.gz"))
    if (save_counts) {
      saveRDS(rna_counts, file.path(stratum_dir, "rna_counts.rds"))
      saveRDS(atac_counts, file.path(stratum_dir, "atac_counts.rds"))
    }
    diag <- data.frame(group_id = key, n_cells = length(cells), n_metacells = length(mc_ids), gamma = gamma_i, requested_gamma = gamma, min_metacell_size = min_metacell_size, min_required_cells = min_required_cells, skipped = FALSE, stringsAsFactors = FALSE)
    .rc_write_tsv_gz(diag, file.path(stratum_dir, "qc", "metacell_qc.tsv.gz"))
    if (requireNamespace("yaml", quietly = TRUE)) {
      run_params <- list(
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
        adaptive_gamma = adaptive_gamma,
        label_col = label_col,
        fragment_files = fragment_manifest,
        bgzip_path = bgzip_path,
        tabix_path = tabix_path,
        fragment_nb_cl = max(1L, as.integer(fragment_nb_cl)),
        save_metacell_object = save_metacell_object,
        save_counts = save_counts,
        save_fragments = save_fragments,
        require_fragment_aggregation = require_fragment_aggregation,
        fragment_aggregation_backend = fragment_aggregation_backend
      )
      yaml::write_yaml(run_params, file.path(stratum_dir, "qc", "run_params.yaml"))
    }
    stratum_dir
  }
  stratum_status_row <- function(key, status, output_dir = NA_character_, error = NULL) {
    cells <- groups[[key]]
    one_meta <- meta[match(cells, meta$cell_id), , drop = FALSE]
    vals <- one_meta[1, group_cols, drop = FALSE]
    target_metacells <- suppressWarnings(floor(length(cells) / as.integer(gamma)))
    if (!is.finite(target_metacells)) target_metacells <- NA_integer_
    actual_metacells <- NA_integer_
    if (!is.na(output_dir)) {
      mm <- file.path(output_dir, "metacell_metadata.tsv.gz")
      if (file.exists(mm)) actual_metacells <- tryCatch(nrow(utils::read.delim(gzfile(mm), stringsAsFactors = FALSE)), error = function(e) NA_integer_)
    }
    data.frame(
      stratum_id = key,
      vals,
      n_input_cells = length(cells),
      gamma = as.integer(gamma),
      target_metacells = as.integer(target_metacells),
      actual_metacells = as.integer(actual_metacells),
      status = status,
      output_dir = output_dir,
      error_class = if (is.null(error)) NA_character_ else class(error)[[1L]],
      error_message = if (is.null(error)) NA_character_ else conditionMessage(error),
      intermediate_files = if (!is.na(output_dir) && dir.exists(output_dir)) length(list.files(output_dir, recursive = TRUE, all.files = FALSE, no.. = TRUE)) else 0L,
      resumable = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  run_one_safe <- function(key) {
    tryCatch(
      {
        out <- run_one(key)
        list(status = stratum_status_row(key, "ok", output_dir = out), output_dir = out)
      },
      error = function(e) {
        if (identical(on_stratum_error, "stop")) stop(e)
        list(status = stratum_status_row(key, "failed", error = e), output_dir = NA_character_)
      }
    )
  }
  results <- rc_parallel_lapply(names(groups), run_one_safe, BPPARAM = BPPARAM)
  status <- do.call(rbind, lapply(results, `[[`, "status"))
  .rc_write_tsv_gz(status, file.path(outdir, "metacell_stratum_status.tsv.gz"))
  dirs <- vapply(results, `[[`, character(1), "output_dir")
  dirs <- dirs[!is.na(dirs) & nzchar(dirs)]
  if (!length(dirs)) stop("All metacell strata failed. See metacell_stratum_status.tsv.gz for details.", call. = FALSE)
  out <- rc_import_supercell2_metacells(dirs, rna_assay = rna_assay, atac_assay = atac_assay, sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col, require_fragments = require_fragment_aggregation)
  out$stratum_status <- status
  out
}

#' Import saved SuperCell2.0 metacell outputs
rc_import_supercell2_metacells <- function(metacell_dirs,
                                           rna_assay = "RNA",
                                           atac_assay = "ATAC",
                                           sample_col = "sample_id",
                                           condition_col = "condition",
                                           celltype_col = "cell_type",
                                           require_fragments = FALSE) {
  metacell_dirs <- metacell_dirs[dir.exists(metacell_dirs)]
  if (length(metacell_dirs) == 0L) stop("No valid metacell directories supplied.", call. = FALSE)
  read_tsv <- function(path) utils::read.delim(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
  metas <- memberships <- fragment_manifests <- list(); rnas <- atacs <- list(); objects <- fragments <- character(0)
  for (d in metacell_dirs) {
    meta_file <- file.path(d, "metacell_metadata.tsv.gz")
    if (!file.exists(meta_file)) next
    mm <- read_tsv(meta_file)
    metas[[d]] <- mm
    mem_file <- file.path(d, "membership.tsv.gz")
    if (file.exists(mem_file)) memberships[[d]] <- read_tsv(mem_file)
    rna_file <- file.path(d, "rna_counts.rds")
    atac_file <- file.path(d, "atac_counts.rds")
    obj_file <- file.path(d, "metacell_object.rds")
    if (file.exists(rna_file)) rnas[[d]] <- readRDS(rna_file)
    if (file.exists(atac_file)) atacs[[d]] <- readRDS(atac_file)
    if (file.exists(obj_file)) objects <- c(objects, obj_file)
    manifest_file <- file.path(d, "fragments", "fragment_manifest.tsv.gz")
    if (file.exists(manifest_file)) {
      fm <- read_tsv(manifest_file)
      fm$stratum_dir <- d
      fragment_manifests[[d]] <- fm
    } else {
      frag <- Sys.glob(file.path(d, "fragments", "*.tsv.gz"))
      frag <- setdiff(frag, manifest_file)
      if (length(frag)) fragments <- c(fragments, frag)
    }
  }
  metacell_meta <- do.call(rbind, metas); rownames(metacell_meta) <- NULL
  metacell_meta$metacell_id <- as.character(metacell_meta$metacell_id)
  if (anyDuplicated(metacell_meta$metacell_id)) {
    duplicated_ids <- unique(metacell_meta$metacell_id[duplicated(metacell_meta$metacell_id)])
    stop("Duplicated metacell IDs across strata: ", paste(utils::head(duplicated_ids, 10L), collapse = ", "), call. = FALSE)
  }
  membership <- if (length(memberships)) do.call(rbind, memberships) else data.frame()
  if (length(rnas) == 0L) stop("No rna_counts.rds files were found in metacell directories.", call. = FALSE)
  rna_counts <- do.call(cbind, lapply(rnas, .rc_as_sparse))
  atac_counts <- if (length(atacs)) do.call(cbind, lapply(atacs, .rc_as_sparse)) else NULL
  colnames(rna_counts) <- as.character(colnames(rna_counts))
  if (!is.null(atac_counts)) {
    colnames(atac_counts) <- as.character(colnames(atac_counts))
    if (!setequal(colnames(rna_counts), colnames(atac_counts))) stop("RNA and ATAC metacell IDs differ after import.", call. = FALSE)
    atac_counts <- atac_counts[, colnames(rna_counts), drop = FALSE]
  }
  metacell_meta <- metacell_meta[match(colnames(rna_counts), metacell_meta$metacell_id), , drop = FALSE]
  if (anyNA(metacell_meta$metacell_id)) stop("Metacell metadata are incomplete.", call. = FALSE)
  rc_validate_metacell_inputs(rna_counts, metacell_meta, atac_metacell_counts = atac_counts, sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col)
  fragment_manifest <- if (length(fragment_manifests)) do.call(rbind, fragment_manifests) else data.frame()
  if (nrow(fragment_manifest)) {
    fragments <- unique(as.character(fragment_manifest$fragment_file))
  } else if (length(fragments)) {
    warning("Legacy fragment discovery by glob was used because no fragment manifest was found.", call. = FALSE)
  }
  if (require_fragments) {
    missing_idx <- fragments[!file.exists(paste0(fragments, ".tbi"))]
    if (length(fragments) == 0L || length(missing_idx) > 0L) stop("Metacell fragment files or tabix indexes are missing.", call. = FALSE)
  }
  list(schema_version = "regcompass_metacell_v1", metacell_meta = metacell_meta, membership = membership, rna_counts = rna_counts, atac_counts = atac_counts, metacell_objects = objects, fragment_manifest = fragment_manifest, fragment_files = fragments, diagnostics = data.frame(n_metacells = ncol(rna_counts), n_membership_rows = nrow(membership)))
}

.rc_apply_used_metacell_ids <- function(out, used_ids) {
  used_ids <- unique(as.character(used_ids))
  if (!length(used_ids)) stop("No metacells remain after filtering.", call. = FALSE)
  missing_meta <- setdiff(used_ids, as.character(out$metacell_meta$metacell_id))
  if (length(missing_meta)) stop("Metacell metadata are missing used metacells: ", paste(utils::head(missing_meta, 10L), collapse = ", "), call. = FALSE)
  missing_rna <- setdiff(used_ids, colnames(out$rna_counts))
  if (length(missing_rna)) stop("RNA counts are missing used metacells: ", paste(utils::head(missing_rna, 10L), collapse = ", "), call. = FALSE)
  out$used_metacell_ids <- used_ids
  out$metacell_meta_used <- out$metacell_meta[match(used_ids, as.character(out$metacell_meta$metacell_id)), , drop = FALSE]
  out$rna_counts <- out$rna_counts[, used_ids, drop = FALSE]
  if (!is.null(out$atac_counts)) {
    missing_atac <- setdiff(used_ids, colnames(out$atac_counts))
    if (length(missing_atac)) stop("ATAC counts are missing used metacells: ", paste(utils::head(missing_atac, 10L), collapse = ", "), call. = FALSE)
    out$atac_counts <- out$atac_counts[, used_ids, drop = FALSE]
  }
  if (is.data.frame(out$membership) && "metacell_id" %in% colnames(out$membership)) {
    out$membership_used <- out$membership[as.character(out$membership$metacell_id) %in% used_ids, , drop = FALSE]
  } else {
    out$membership_used <- out$membership
  }
  if (is.data.frame(out$fragment_manifest) && nrow(out$fragment_manifest) && "object_cell" %in% colnames(out$fragment_manifest)) {
    out$fragment_manifest_used <- out$fragment_manifest[as.character(out$fragment_manifest$object_cell) %in% used_ids, , drop = FALSE]
  } else {
    out$fragment_manifest_used <- out$fragment_manifest
  }
  out
}



.rc_seurat4_filterobjects <- function(object, classes.keep = c("Assay", "Assay5", "ChromatinAssay")) {
  assays <- names(object@assays)
  assays[vapply(assays, function(a) {
    any(vapply(classes.keep, function(cl) inherits(object@assays[[a]], cl), logical(1)))
  }, logical(1))]
}

.rc_install_seurat4_filterobjects <- function(envir = .GlobalEnv) {
  if (!exists(".FilterObjects", envir = envir, inherits = FALSE)) {
    assign(".FilterObjects", .rc_seurat4_filterobjects, envir = envir)
  }
  invisible(TRUE)
}

.rc_with_seurat4_filterobjects <- function(expr) {
  envir <- .GlobalEnv
  had_old <- exists(".FilterObjects", envir = envir, inherits = FALSE)
  old <- if (had_old) get(".FilterObjects", envir = envir, inherits = FALSE) else NULL
  .rc_install_seurat4_filterobjects(envir)
  on.exit({
    if (had_old) assign(".FilterObjects", old, envir = envir) else if (exists(".FilterObjects", envir = envir, inherits = FALSE)) rm(".FilterObjects", envir = envir)
  }, add = TRUE)
  force(expr)
}
#' Load and merge saved metacell Seurat objects
rc_load_or_merge_metacell_objects <- function(metacell_objects, fragment_manifest = NULL, metacell_meta = NULL, fragment_files = NULL, rna_assay = "RNA", atac_assay = "ATAC", require_complete_fragments = TRUE) {
  if (is.null(metacell_objects) || length(metacell_objects) == 0L) stop("No metacell Seurat objects supplied.", call. = FALSE)
  objs <- lapply(metacell_objects, function(x) if (inherits(x, "Seurat")) x else readRDS(x))

  object_cells_by_input <- lapply(objs, colnames)
  input_cells <- unlist(object_cells_by_input, use.names = FALSE)
  duplicated_before_merge <- unique(input_cells[duplicated(input_cells)])
  if (length(duplicated_before_merge)) {
    stop("Metacell IDs are not globally unique before merge: ", paste(utils::head(duplicated_before_merge, 10L), collapse = ", "), call. = FALSE)
  }

  objs <- lapply(objs, .rc_clear_signac_fragments, atac_assay = atac_assay)
  obj <- if (length(objs) == 1L) objs[[1L]] else Reduce(function(a, b) merge(x = a, y = b, merge.data = FALSE), objs)
  if (anyDuplicated(colnames(obj))) stop("Merged metacell object contains duplicated cell names.", call. = FALSE)

  if (!is.null(metacell_meta)) {
    metacell_meta$metacell_id <- as.character(metacell_meta$metacell_id)
    expected <- metacell_meta$metacell_id
    observed <- colnames(obj)
    missing_in_object <- setdiff(expected, observed)
    if (length(missing_in_object)) stop("Merged metacell object is missing expected IDs: ", paste(utils::head(missing_in_object, 10L), collapse = ", "), call. = FALSE)
    extra_in_object <- setdiff(observed, expected)
    obj <- subset(obj, cells = expected)
    if (!identical(colnames(obj), expected)) stop("Merged object could not be subset and reordered to expected metacell IDs.", call. = FALSE)
    attr(obj, "removed_extra_metacell_ids") <- extra_in_object
  }

  if (!is.null(fragment_manifest) && is.data.frame(fragment_manifest) && nrow(fragment_manifest)) {
    registration <- .rc_fragment_registration_from_manifest(
      fragment_manifest = fragment_manifest,
      metacell_meta = metacell_meta,
      object_cells = colnames(obj)
    )
    fragment_files <- registration$fragment_files
    cells_by_fragment <- registration$cell_maps
  } else {
    cells_by_fragment <- NULL
  }

  .rc_register_signac_fragments(
    obj,
    fragment_files = fragment_files,
    cells_by_fragment = cells_by_fragment,
    atac_assay = atac_assay,
    replace_existing = TRUE,
    require_complete = require_complete_fragments
  )
}

.rc_clear_signac_fragments <- function(object, atac_assay = "ATAC") {
  if (!inherits(object, "Seurat")) return(object)
  if (!requireNamespace("Signac", quietly = TRUE)) return(object)
  if (!atac_assay %in% names(object@assays)) return(object)
  if (!inherits(object[[atac_assay]], "ChromatinAssay")) return(object)
  fragment_setter <- get("Fragments<-", envir = asNamespace("Signac"))
  object[[atac_assay]] <- fragment_setter(object[[atac_assay]], value = list())
  object
}

.rc_normalize_fragment_cell_map <- function(cell_map, object_cells, fragment_file = NULL) {
  if (is.data.frame(cell_map)) {
    required <- c("object_cell", "fragment_barcode")
    missing <- setdiff(required, colnames(cell_map))
    if (length(missing)) stop("`cell_map` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    cell_map <- stats::setNames(as.character(cell_map$fragment_barcode), as.character(cell_map$object_cell))
  } else {
    cell_map <- as.character(cell_map)
    if (is.null(names(cell_map))) names(cell_map) <- cell_map
  }
  if (!length(cell_map)) stop("Fragment cell map is empty", if (!is.null(fragment_file)) paste0(": ", fragment_file) else ".", call. = FALSE)
  if (anyNA(cell_map) || any(!nzchar(cell_map)) || anyNA(names(cell_map)) || any(!nzchar(names(cell_map)))) stop("Fragment cell map contains missing or empty identifiers.", call. = FALSE)
  if (anyDuplicated(names(cell_map))) {
    duplicated_cells <- unique(names(cell_map)[duplicated(names(cell_map))])
    stop("Duplicated object cells within one fragment mapping: ", paste(utils::head(duplicated_cells, 10L), collapse = ", "), call. = FALSE)
  }
  unknown <- setdiff(names(cell_map), object_cells)
  if (length(unknown)) stop("Fragment mapping contains cells absent from the merged object: ", paste(utils::head(unknown, 10L), collapse = ", "), call. = FALSE)
  cell_map
}

.rc_validate_fragment_registration_plan <- function(fragment_files, cell_maps, object_cells, require_complete = TRUE) {
  fragment_files <- as.character(fragment_files)
  if (length(fragment_files) != length(cell_maps)) stop("`fragment_files` and `cell_maps` must have the same length.", call. = FALSE)
  if (!length(fragment_files)) stop("No fragment files were supplied.", call. = FALSE)
  missing_files <- fragment_files[!file.exists(fragment_files)]
  if (length(missing_files)) stop("Metacell fragment files are missing: ", paste(utils::head(missing_files, 10L), collapse = ", "), call. = FALSE)
  missing_indexes <- vapply(fragment_files, function(path) !file.exists(paste0(path, ".tbi")) && !file.exists(paste0(path, ".csi")), logical(1))
  if (any(missing_indexes)) stop("Metacell fragment tabix indexes are missing: ", paste(utils::head(paste0(fragment_files[missing_indexes], ".tbi"), 10L), collapse = ", "), call. = FALSE)
  registered <- unlist(lapply(cell_maps, names), use.names = FALSE)
  if (anyDuplicated(registered)) {
    duplicated_cells <- unique(registered[duplicated(registered)])
    stop("Object cells are assigned to multiple fragment files: ", paste(utils::head(duplicated_cells, 10L), collapse = ", "), call. = FALSE)
  }
  if (isTRUE(require_complete)) {
    missing_cells <- setdiff(object_cells, registered)
    extra_cells <- setdiff(registered, object_cells)
    if (length(missing_cells) || length(extra_cells)) stop("Fragment registration does not exactly cover the merged object. Missing cells: ", paste(utils::head(missing_cells, 10L), collapse = ", "), "; extra cells: ", paste(utils::head(extra_cells, 10L), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.rc_fragment_registration_from_manifest <- function(fragment_manifest, metacell_meta = NULL, object_cells) {
  if (!"fragment_file" %in% colnames(fragment_manifest)) stop("`fragment_manifest` must contain `fragment_file`.", call. = FALSE)
  manifest <- fragment_manifest
  manifest$fragment_file <- as.character(manifest$fragment_file)
  if (all(c("object_cell", "fragment_barcode") %in% colnames(manifest))) {
    manifest$object_cell <- as.character(manifest$object_cell)
    manifest$fragment_barcode <- as.character(manifest$fragment_barcode)
    manifest <- manifest[manifest$object_cell %in% object_cells, , drop = FALSE]
    missing_maps <- setdiff(object_cells, unique(manifest$object_cell))
    if (length(missing_maps)) stop("Fragment manifest is missing mappings for metacells: ", paste(utils::head(missing_maps, 10L), collapse = ", "), call. = FALSE)
    manifest <- unique(manifest[, c("fragment_file", "object_cell", "fragment_barcode"), drop = FALSE])
    cell_path <- paste(manifest$fragment_file, manifest$object_cell, sep = "\001")
    barcode_by_cell_path <- tapply(manifest$fragment_barcode, cell_path, function(x) length(unique(x)))
    conflicts <- names(barcode_by_cell_path)[barcode_by_cell_path > 1L]
    if (length(conflicts)) {
      conflict_cells <- sub("^.*\\001", "", conflicts)
      stop("Fragment manifest assigns one object cell to multiple barcodes: ", paste(utils::head(conflict_cells, 10L), collapse = ", "), call. = FALSE)
    }
    files <- unique(manifest$fragment_file)
    maps <- lapply(files, function(path) {
      x <- manifest[manifest$fragment_file == path, , drop = FALSE]
      .rc_normalize_fragment_cell_map(x[, c("object_cell", "fragment_barcode"), drop = FALSE], object_cells = object_cells, fragment_file = path)
    })
    return(list(fragment_files = files, cell_maps = maps))
  }
  stop(
    paste(
      "Fragment manifest entries must contain explicit",
      "`object_cell` and `fragment_barcode` columns."
    ),
    call. = FALSE
  )
}

.rc_register_signac_fragments <- function(object, fragment_files = NULL, cells_by_fragment = NULL, atac_assay = "ATAC", replace_existing = TRUE, require_complete = TRUE, validate_fragments = TRUE) {
  if (is.null(fragment_files) || length(fragment_files) == 0L) return(object)
  fragment_files <- as.character(fragment_files)
  if (is.null(cells_by_fragment)) {
    if (length(fragment_files) != 1L) stop("`cells_by_fragment` is required when registering multiple fragment files.", call. = FALSE)
    ids <- colnames(object)
    cells_by_fragment <- list(stats::setNames(ids, ids))
  }
  if (length(cells_by_fragment) != length(fragment_files)) stop("`cells_by_fragment` must have one cell vector per fragment file.", call. = FALSE)
  cell_maps <- Map(function(cell_map, path) .rc_normalize_fragment_cell_map(cell_map, object_cells = colnames(object), fragment_file = path), cells_by_fragment, fragment_files)
  .rc_validate_fragment_registration_plan(fragment_files, cell_maps, object_cells = colnames(object), require_complete = require_complete)
  if (!requireNamespace("Signac", quietly = TRUE)) stop("Package 'Signac' is required to register metacell fragment files.", call. = FALSE)
  if (!inherits(object, "Seurat") || !atac_assay %in% names(object@assays)) stop("Metacell object is missing ATAC assay `", atac_assay, "`.", call. = FALSE)
  fragment_files <- normalizePath(fragment_files, mustWork = TRUE)
  frag_setter <- get("Fragments<-", envir = asNamespace("Signac"))
  if (isTRUE(replace_existing)) object[[atac_assay]] <- frag_setter(object[[atac_assay]], value = list())
  fragments <- Map(function(path, cell_map) {
    tryCatch(
      Signac::CreateFragmentObject(path = path, cells = cell_map, validate.fragments = validate_fragments),
      error = function(e) stop("Failed to register metacell fragment file `", path, "`: ", conditionMessage(e), call. = FALSE)
    )
  }, fragment_files, cell_maps)
  object[[atac_assay]] <- frag_setter(object[[atac_assay]], value = fragments)
  registered <- unlist(lapply(fragments, SeuratObject::Cells), use.names = FALSE)
  if (anyDuplicated(registered)) stop("Post-registration validation detected cells in multiple Fragment objects.", call. = FALSE)
  object
}
