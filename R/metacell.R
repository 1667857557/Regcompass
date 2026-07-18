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
  .rc_expand_fragment_manifest(
    manifest,
    unique(as.character(membership$metacell_id))
  )
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
  if (is.null(path) && methods::is(fragment, "Fragment") && isS4(fragment) && "path" %in% methods::slotNames(fragment)) {
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

.rc_fragments_requested <- function(fragment_files) {
  !identical(fragment_files, FALSE)
}

.rc_fragment_processing_enabled <- function(fragment_files, save_fragments,
                                            backend,
                                            require_fragments = TRUE) {
  .rc_fragments_requested(fragment_files) &&
    (!is.null(fragment_files) || isTRUE(require_fragments)) &&
    isTRUE(save_fragments) &&
    !identical(backend, "none")
}

#' Build sample-aware SuperCell2.0 RNA+ATAC metacells and save outputs
.rc_build_supercell2_strata <- function(object,
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
  controls <- c(
    gamma = gamma,
    min_cells_per_stratum = min_cells_per_stratum,
    min_metacell_size = min_metacell_size,
    min_metacells_per_stratum = min_metacells_per_stratum,
    fragment_nb_cl = fragment_nb_cl
  )
  if (any(!is.finite(controls)) || any(controls < 1) ||
      any(abs(controls - round(controls)) > sqrt(.Machine$double.eps))) {
    stop(
      "Metacell size, count, gamma, and worker controls must be positive integers.",
      call. = FALSE
    )
  }
  .rc_require_supercell2()
  if (!inherits(object, "Seurat")) stop("`object` must inherit from class 'Seurat'.", call. = FALSE)
  if (identical(fragment_files, FALSE)) {
    fragment_files <- NULL
    save_fragments <- FALSE
    require_fragment_aggregation <- FALSE
    fragment_aggregation_backend <- "none"
  }
  if (is.null(fragment_files) && isTRUE(require_fragment_aggregation)) {
    fragment_files <- .rc_fragment_files_from_atac(object, atac_assay = atac_assay)
  }
  if (isTRUE(require_fragment_aggregation)) {
    if (!isTRUE(save_fragments)) stop("Formal multiome workflow requires `save_fragments = TRUE`.", call. = FALSE)
    if (is.null(fragment_files)) stop("Formal multiome workflow requires `fragment_files` for metacell fragment aggregation, or a fragment file registered on the ATAC assay. Use `fragment_files = FALSE` to skip fragment aggregation and use ATAC peak raw counts from the object.", call. = FALSE)
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
    checkpoint_files <- c(
      file.path(stratum_dir, "metacell_metadata.tsv.gz"),
      file.path(stratum_dir, "rna_counts.rds"),
      file.path(stratum_dir, "atac_counts.rds")
    )
    if (isTRUE(save_metacell_object)) {
      checkpoint_files <- c(
        checkpoint_files,
        file.path(stratum_dir, "metacell_object.rds")
      )
    }
    if (isTRUE(require_fragment_aggregation)) {
      checkpoint_files <- c(
        checkpoint_files,
        file.path(stratum_dir, "fragments", "fragment_manifest.tsv.gz")
      )
    }
    if (dir.exists(stratum_dir) && !overwrite &&
        all(file.exists(checkpoint_files))) {
      return(stratum_dir)
    }
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
    if (length(mc_ids) < as.integer(min_metacells_per_stratum)) {
      diag <- data.frame(
        group_id = key,
        n_cells = length(cells),
        n_metacells = length(mc_ids),
        gamma = gamma_i,
        requested_gamma = gamma,
        min_metacells_per_stratum = as.integer(min_metacells_per_stratum),
        skipped = TRUE,
        skip_reason = "stratum_below_min_metacells_per_stratum",
        stringsAsFactors = FALSE
      )
      .rc_write_tsv_gz(
        diag,
        file.path(stratum_dir, "qc", "metacell_qc.tsv.gz")
      )
      return(stratum_dir)
    }
    display_ids <- paste0(prefix, "_MC", sprintf(paste0("%0", max(3, nchar(length(mc_ids))), "d"), seq_along(mc_ids)))
    membership <- .rc_extract_supercell_membership(mc, cells, mc_ids)
    if (nrow(membership) == 0L) stop("Could not infer single-cell membership from SuperCell output for stratum ", key, call. = FALSE)
    fragment_manifest_i <- NULL
    if (identical(fragment_aggregation_backend, "supercell") && save_fragments) {
      fragment_manifest_i <- tryCatch(mc@misc$fragment_manifest, error = function(e) NULL)
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
        completed <- file.exists(file.path(out, "metacell_metadata.tsv.gz")) &&
          file.exists(file.path(out, "rna_counts.rds"))
        state <- if (completed) "ok" else "skipped"
        list(
          status = stratum_status_row(key, state, output_dir = out),
          output_dir = if (completed) out else NA_character_
        )
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

#' Build sample-aware SuperCell2.0 RNA+ATAC metacells and save outputs
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
    on_stratum_error = c("record", "stop"),
    call_peaks_from_fragments = TRUE,
    macs2_path = NULL,
    peak_calling_effective_genome_size = NULL,
    peak_calling_args = list()) {
  fragment_aggregation_backend <- match.arg(fragment_aggregation_backend)
  on_stratum_error <- match.arg(on_stratum_error)
  if (!is.logical(call_peaks_from_fragments) ||
      length(call_peaks_from_fragments) != 1L ||
      is.na(call_peaks_from_fragments)) {
    stop("call_peaks_from_fragments must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(peak_calling_args)) {
    stop("peak_calling_args must be a list.", call. = FALSE)
  }
  if (!isTRUE(save_counts)) {
    stop(
      "save_counts must be TRUE because the returned bundle is imported from ",
      "the per-stratum count artifacts.",
      call. = FALSE
    )
  }
  fragment_enabled <- .rc_fragment_processing_enabled(
    fragment_files,
    save_fragments,
    fragment_aggregation_backend,
    require_fragments = require_fragment_aggregation
  )
  if (fragment_enabled && !isTRUE(save_metacell_object)) {
    stop(
      "save_metacell_object must be TRUE when fragment recounting is enabled.",
      call. = FALSE
    )
  }

  built <- .rc_build_supercell2_strata(
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

  if (!fragment_enabled) {
    built$atac_count_source <- "aggregated_object_peak_counts"
    built$atac_peak_source <- "existing_object_peak_ranges"
    return(built)
  }
  if (!is.data.frame(built$fragment_manifest) ||
      !nrow(built$fragment_manifest)) {
    if (isTRUE(require_fragment_aggregation)) {
      stop(
        "Fragment aggregation completed without an explicit usable manifest.",
        call. = FALSE
      )
    }
    warning(
      "No fragment manifest was produced; retaining aggregated ATAC peak ",
      "counts because fragment aggregation was not required.",
      call. = FALSE
    )
    built$atac_count_source <- "aggregated_object_peak_counts"
    built$atac_peak_source <- "existing_object_peak_ranges"
    return(built)
  }

  object_files <- as.character(built$metacell_objects)
  if (!length(object_files)) {
    stop("Fragment recounting requires saved metacell objects.", call. = FALSE)
  }
  for (object_file in object_files) {
    mc <- readRDS(object_file)
    stratum_dir <- dirname(object_file)
    manifest_i <- built$fragment_manifest
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
    already_recounted <- identical(
      tryCatch(mc@misc$atac_count_source, error = function(e) NULL),
      "recomputed_from_metacell_fragments"
    )
    expected_peak_source <- if (isTRUE(call_peaks_from_fragments)) {
      "de_novo_macs2_from_metacell_fragments"
    } else {
      "existing_object_peak_ranges"
    }
    already_recounted <- already_recounted && identical(
      tryCatch(mc@misc$atac_peak_source, error = function(e) NULL),
      expected_peak_source
    )
    if (already_recounted && !isTRUE(overwrite)) next
    mc <- .rc_recount_atac_from_fragment_manifest(
      object = mc,
      fragment_manifest = manifest_i,
      atac_assay = atac_assay,
      require_complete = TRUE,
      call_peaks = call_peaks_from_fragments,
      macs2_path = macs2_path,
      effective_genome_size = peak_calling_effective_genome_size,
      peak_calling_args = peak_calling_args,
      peak_calling_outdir = file.path(stratum_dir, "peaks", "macs2")
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
    if (identical(
      tryCatch(mc@misc$atac_peak_source, error = function(e) NULL),
      "de_novo_macs2_from_metacell_fragments"
    )) {
      peak_dir <- file.path(stratum_dir, "peaks")
      dir.create(peak_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(
        methods::slot(mc[[atac_assay]], "ranges"),
        file.path(peak_dir, "called_peaks.rds")
      )
    }
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
  refreshed$stratum_status <- built$stratum_status
  refreshed$atac_count_source <- "recomputed_from_metacell_fragments"
  refreshed$atac_peak_source <- if (isTRUE(call_peaks_from_fragments)) {
    "de_novo_macs2_from_metacell_fragments"
  } else {
    "existing_object_peak_ranges"
  }
  refreshed
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
