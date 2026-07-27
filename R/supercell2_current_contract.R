# Current SuperCell2 contract: explicit RegCompass strata plus the SuperCell2
# `label` argument. No biological-sample column or artificial condition-pool
# metadata field is created or interpreted.

rc_validate_metacell_inputs <- function(
    rna_metacell_counts, metacell_meta, atac_metacell_counts = NULL,
    metacell_id_col = "metacell_id", condition_col = "condition",
    celltype_col = "cell_type") {
  if (is.null(dim(rna_metacell_counts)) || length(dim(rna_metacell_counts)) != 2L) {
    stop("`rna_metacell_counts` must be a feature-by-metacell matrix.", call. = FALSE)
  }
  ids <- colnames(rna_metacell_counts)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("RNA metacell IDs must be unique and non-empty.", call. = FALSE)
  }
  if (!is.data.frame(metacell_meta)) {
    stop("`metacell_meta` must be a data.frame.", call. = FALSE)
  }
  required <- unique(c(metacell_id_col, condition_col, celltype_col))
  required <- required[!is.na(required) & nzchar(required)]
  missing <- setdiff(required, colnames(metacell_meta))
  if (length(missing)) {
    stop("`metacell_meta` is missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  meta_ids <- trimws(as.character(metacell_meta[[metacell_id_col]]))
  if (anyNA(meta_ids) || any(!nzchar(meta_ids)) || anyDuplicated(meta_ids)) {
    stop("Metacell metadata IDs must be unique and non-empty.", call. = FALSE)
  }
  if (!setequal(ids, meta_ids)) {
    stop("RNA counts and metacell metadata contain different metacell IDs.",
         call. = FALSE)
  }
  for (column in setdiff(required, metacell_id_col)) {
    value <- trimws(as.character(metacell_meta[[column]]))
    if (anyNA(value) || any(!nzchar(value))) {
      stop("Metacell metadata column `", column, "` is incomplete.", call. = FALSE)
    }
  }
  if (!is.null(atac_metacell_counts)) {
    atac_ids <- colnames(atac_metacell_counts)
    if (is.null(dim(atac_metacell_counts)) || length(dim(atac_metacell_counts)) != 2L ||
        is.null(atac_ids) || !identical(ids, atac_ids)) {
      stop("RNA and ATAC metacell matrices must align exactly.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

rc_build_metacell_metadata <- function(
    membership, metacell_id_col = "metacell_id", cell_id_col = "cell_id",
    strict_cols = NULL) {
  if (!is.data.frame(membership)) {
    stop("`membership` must be a data.frame.", call. = FALSE)
  }
  required <- c(metacell_id_col, cell_id_col)
  missing <- setdiff(required, colnames(membership))
  if (length(missing)) {
    stop("`membership` is missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  x <- membership
  x[[metacell_id_col]] <- trimws(as.character(x[[metacell_id_col]]))
  x[[cell_id_col]] <- trimws(as.character(x[[cell_id_col]]))
  keep <- !is.na(x[[metacell_id_col]]) & nzchar(x[[metacell_id_col]]) &
    !is.na(x[[cell_id_col]]) & nzchar(x[[cell_id_col]])
  x <- x[keep, , drop = FALSE]
  if (anyDuplicated(x[[cell_id_col]])) {
    stop("A single cell cannot belong to more than one metacell.", call. = FALSE)
  }
  strict_cols <- unique(as.character(strict_cols %||% character()))
  strict_cols <- strict_cols[!is.na(strict_cols) & nzchar(strict_cols)]
  missing_strict <- setdiff(strict_cols, colnames(x))
  if (length(missing_strict)) {
    stop("Membership lacks strict provenance columns: ",
         paste(missing_strict, collapse = ", "), call. = FALSE)
  }
  split_rows <- split(seq_len(nrow(x)), x[[metacell_id_col]])
  for (id in names(split_rows)) {
    rows <- split_rows[[id]]
    for (column in strict_cols) {
      value <- trimws(as.character(x[[column]][rows]))
      unique_value <- unique(value[!is.na(value) & nzchar(value)])
      if (length(unique_value) != 1L || anyNA(value) || any(!nzchar(value))) {
        stop("Metacell `", id, "` mixes or lacks `", column, "` labels.",
             call. = FALSE)
      }
    }
  }
  columns <- unique(c(metacell_id_col, strict_cols))
  out <- x[!duplicated(x[[metacell_id_col]]), columns, drop = FALSE]
  out$n_cells <- as.integer(vapply(
    as.character(out[[metacell_id_col]]),
    function(id) length(split_rows[[id]]), integer(1)
  ))
  rownames(out) <- NULL
  out
}

rc_import_supercell2_metacells <- function(
    metacell_dirs, rna_assay = "RNA", atac_assay = "ATAC",
    condition_col = "condition", celltype_col = "cell_type",
    require_fragments = FALSE) {
  metacell_dirs <- unique(as.character(metacell_dirs))
  metacell_dirs <- metacell_dirs[dir.exists(metacell_dirs)]
  if (!length(metacell_dirs)) stop("No valid metacell directories supplied.", call. = FALSE)
  read_tsv <- function(path) {
    out <- utils::read.delim(gzfile(path), stringsAsFactors = FALSE,
                             check.names = FALSE, colClasses = "character")
    out[] <- lapply(out, utils::type.convert, as.is = TRUE, tryLogical = FALSE)
    out
  }
  metas <- memberships <- fragment_manifests <- rnas <- atacs <- list()
  objects <- character()
  for (directory in metacell_dirs) {
    meta_file <- file.path(directory, "metacell_metadata.tsv.gz")
    rna_file <- file.path(directory, "rna_counts.rds")
    if (!file.exists(meta_file) || !file.exists(rna_file)) next
    metas[[directory]] <- read_tsv(meta_file)
    rnas[[directory]] <- readRDS(rna_file)
    membership_file <- file.path(directory, "membership.tsv.gz")
    atac_file <- file.path(directory, "atac_counts.rds")
    object_file <- file.path(directory, "metacell_object.rds")
    manifest_file <- file.path(directory, "fragments", "fragment_manifest.tsv.gz")
    if (file.exists(membership_file)) memberships[[directory]] <- read_tsv(membership_file)
    if (file.exists(atac_file)) atacs[[directory]] <- readRDS(atac_file)
    if (file.exists(object_file)) objects <- c(objects, object_file)
    if (file.exists(manifest_file)) {
      one <- read_tsv(manifest_file)
      one$stratum_dir <- directory
      fragment_manifests[[directory]] <- one
    }
  }
  if (!length(metas) || !length(rnas)) {
    stop("No complete metacell strata were found.", call. = FALSE)
  }
  metacell_meta <- do.call(rbind, metas)
  rownames(metacell_meta) <- NULL
  membership <- if (length(memberships)) do.call(rbind, memberships) else data.frame()
  rna_counts <- .rc_cbind_sparse_feature_union(rnas)
  atac_counts <- .rc_cbind_sparse_feature_union(atacs)
  if (!is.null(atac_counts)) atac_counts <- atac_counts[, colnames(rna_counts), drop = FALSE]
  metacell_meta <- metacell_meta[
    match(colnames(rna_counts), as.character(metacell_meta$metacell_id)),
    , drop = FALSE
  ]
  rc_validate_metacell_inputs(
    rna_counts, metacell_meta, atac_counts,
    condition_col = condition_col, celltype_col = celltype_col
  )
  fragment_manifest <- if (length(fragment_manifests)) {
    do.call(rbind, fragment_manifests)
  } else data.frame()
  fragment_files <- if (nrow(fragment_manifest) &&
                        "fragment_file" %in% colnames(fragment_manifest)) {
    unique(as.character(fragment_manifest$fragment_file))
  } else character()
  if (isTRUE(require_fragments) && (!length(fragment_files) ||
      any(!file.exists(fragment_files)))) {
    stop("Required metacell fragment files are unavailable.", call. = FALSE)
  }
  list(
    schema_version = "regcompass_metacell_v2_supercell_current",
    metacell_meta = metacell_meta,
    membership = membership,
    rna_counts = rna_counts,
    atac_counts = atac_counts,
    metacell_objects = objects,
    fragment_manifest = fragment_manifest,
    fragment_files = fragment_files,
    diagnostics = data.frame(
      n_metacells = ncol(rna_counts),
      n_membership_rows = nrow(membership),
      n_atac_peaks = if (is.null(atac_counts)) 0L else nrow(atac_counts),
      stringsAsFactors = FALSE
    )
  )
}

.rc_supercell2_fragment_files_for_stratum <- function(fragment_files, key) {
  if (identical(fragment_files, FALSE) || is.null(fragment_files)) return(character())
  if (is.data.frame(fragment_files)) {
    required <- c("stratum_id", "fragment_file")
    missing <- setdiff(required, colnames(fragment_files))
    if (length(missing)) {
      stop("Fragment manifest must contain stratum_id and fragment_file.", call. = FALSE)
    }
    return(as.character(fragment_files$fragment_file[
      as.character(fragment_files$stratum_id) == key
    ]))
  }
  if (is.list(fragment_files)) {
    if (!is.null(names(fragment_files)) && key %in% names(fragment_files)) {
      return(as.character(unlist(fragment_files[[key]], use.names = FALSE)))
    }
    value <- as.character(unlist(fragment_files, use.names = FALSE))
    if (length(unique(value)) == 1L) return(unique(value))
    stop("Named fragment files must use exact RegCompass stratum IDs.", call. = FALSE)
  }
  value <- as.character(fragment_files)
  if (length(value) == 1L) return(value)
  stop("Multiple fragment files require a stratum_id manifest or named list.",
       call. = FALSE)
}

#' Build SuperCell2 metacells from explicit RegCompass strata
#'
#' `strata_cols` are used by RegCompass to split the input object before each
#' SuperCell2 call. `label_col` is passed through the current
#' `SuperCell::SCimplify_for_Seurat(label = ...)` API. No sample-column argument
#' is sent to SuperCell2 or required in metadata.
#'
#' @export
rc_make_supercell2_metacells <- function(
    object, outdir, strata_cols = "condition", label_col = NULL,
    rna_assay = "RNA", atac_assay = "ATAC",
    rna_reduction = "pca", atac_reduction = "lsi",
    rna_dims = 1:30, atac_dims = 2:30,
    gamma = 30, seed = 12345L,
    min_cells_per_stratum = 100L, min_metacell_size = 20L,
    min_metacells_per_stratum = 2L,
    fragment_files = FALSE, bgzip_path = "bgzip", tabix_path = "tabix",
    fragment_nb_cl = 1L, save_metacell_object = TRUE, save_counts = TRUE,
    overwrite = FALSE, BPPARAM = NULL,
    on_stratum_error = c("record", "stop"), ...) {
  on_stratum_error <- match.arg(on_stratum_error)
  .rc_require_supercell2()
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  if (!isTRUE(save_counts)) stop("`save_counts = TRUE` is required.", call. = FALSE)
  strata_cols <- unique(trimws(as.character(strata_cols)))
  strata_cols <- strata_cols[!is.na(strata_cols) & nzchar(strata_cols)]
  if (!length(strata_cols)) stop("`strata_cols` must contain at least one column.", call. = FALSE)
  required <- unique(c(strata_cols, label_col))
  missing <- setdiff(required, colnames(object@meta.data))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  .rc_validate_supercell2_inputs(
    object, assays = c(rna_assay, atac_assay),
    reductions = c(rna_reduction, atac_reduction)
  )
  meta <- object@meta.data
  invalid <- vapply(meta[, required, drop = FALSE], function(value) {
    anyNA(value) || any(!nzchar(trimws(as.character(value))))
  }, logical(1))
  if (any(invalid)) stop("SuperCell2 strata and labels must be complete.", call. = FALSE)
  meta$cell_id <- rownames(meta)
  meta$.rc_stratum_id <- rc_make_stratum_id(meta, strata_cols)
  groups <- split(meta$cell_id, meta$.rc_stratum_id)
  dots <- list(...)
  supercell_formals <- names(formals(getExportedValue("SuperCell", "SCimplify_for_Seurat")))
  unknown <- setdiff(names(dots), supercell_formals)
  if (length(unknown)) {
    stop("Unknown SuperCell2 arguments: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  run_one <- function(key) {
    cells <- groups[[key]]
    one_meta <- meta[match(cells, meta$cell_id), , drop = FALSE]
    values <- one_meta[1L, strata_cols, drop = FALSE]
    directory <- file.path(outdir, paste0("stratum=", .rc_safe_path_component(key)))
    checkpoints <- file.path(directory, c(
      "metacell_metadata.tsv.gz", "membership.tsv.gz",
      "rna_counts.rds", "atac_counts.rds"
    ))
    if (dir.exists(directory) && !isTRUE(overwrite) && all(file.exists(checkpoints))) {
      return(directory)
    }
    dir.create(file.path(directory, "qc"), recursive = TRUE, showWarnings = FALSE)
    if (length(cells) < as.integer(min_cells_per_stratum)) {
      .rc_write_tsv_gz(data.frame(
        stratum_id = key, n_cells = length(cells), status = "skipped",
        reason = "below_min_cells_per_stratum", stringsAsFactors = FALSE
      ), file.path(directory, "qc", "metacell_qc.tsv.gz"))
      return(NA_character_)
    }
    prefix <- paste0(.rc_safe_path_component(key), "_MC_")
    args <- c(list(
      seurat = subset(object, cells = cells),
      assay = c(rna_assay, atac_assay),
      reduction = list(rna_reduction, atac_reduction),
      dims = list(rna_dims, atac_dims),
      gamma = gamma,
      prefixMC = prefix,
      seed = as.integer(seed) + match(key, names(groups)) - 1L,
      return.seurat = TRUE,
      return_membership_table = TRUE
    ), dots)
    if (!is.null(label_col)) args$label <- label_col
    fragments <- .rc_supercell2_fragment_files_for_stratum(fragment_files, key)
    fragments <- unique(fragments[!is.na(fragments) & nzchar(fragments)])
    if (length(fragments)) {
      args$fragmentFiles <- stats::setNames(list(fragments), atac_assay)
      args$outputDirMcFragment <- file.path(directory, "fragments")
      args$bgzip_path <- bgzip_path
      args$tabix_path <- tabix_path
      args$nb_cl <- as.integer(fragment_nb_cl)
    }
    mc <- .rc_supercell2_scimplify_for_seurat(args)
    mc_ids <- as.character(colnames(.rc_get_assay_counts_safe(mc, rna_assay)))
    if (length(mc_ids) < as.integer(min_metacells_per_stratum)) {
      .rc_write_tsv_gz(data.frame(
        stratum_id = key, n_cells = length(cells), n_metacells = length(mc_ids),
        status = "skipped", reason = "below_min_metacells_per_stratum",
        stringsAsFactors = FALSE
      ), file.path(directory, "qc", "metacell_qc.tsv.gz"))
      return(NA_character_)
    }
    membership <- .rc_extract_supercell_membership(mc, cells, mc_ids)
    if (!nrow(membership)) stop("SuperCell2 returned no membership table for `", key, "`.", call. = FALSE)
    original_index <- match(membership$cell_id, rownames(meta))
    for (column in required) membership[[column]] <- meta[[column]][original_index]
    mc_meta <- rc_build_metacell_metadata(
      membership, strict_cols = required
    )
    mc_meta$stratum_id <- key
    mc_meta$low_power_metacell <- mc_meta$n_cells < as.integer(min_metacell_size)
    mc_meta$effective_gamma <- as.integer(gamma)
    rna_counts <- .rc_as_sparse(.rc_get_assay_counts_safe(mc, rna_assay))
    atac_counts <- .rc_as_sparse(.rc_get_assay_counts_safe(mc, atac_assay))
    if (isTRUE(save_metacell_object)) saveRDS(mc, file.path(directory, "metacell_object.rds"))
    .rc_write_tsv_gz(membership, file.path(directory, "membership.tsv.gz"))
    .rc_write_tsv_gz(mc_meta, file.path(directory, "metacell_metadata.tsv.gz"))
    saveRDS(rna_counts, file.path(directory, "rna_counts.rds"))
    saveRDS(atac_counts, file.path(directory, "atac_counts.rds"))
    manifest <- tryCatch(mc@misc$fragment_manifest, error = function(e) NULL)
    if (is.data.frame(manifest) && nrow(manifest)) {
      manifest$stratum_id <- key
      dir.create(file.path(directory, "fragments"), recursive = TRUE, showWarnings = FALSE)
      .rc_write_tsv_gz(manifest, file.path(directory, "fragments", "fragment_manifest.tsv.gz"))
    }
    .rc_write_tsv_gz(data.frame(
      stratum_id = key, n_cells = length(cells), n_metacells = length(mc_ids),
      gamma = as.integer(gamma), status = "ok", stringsAsFactors = FALSE
    ), file.path(directory, "qc", "metacell_qc.tsv.gz"))
    directory
  }

  run_safe <- function(key) tryCatch(
    list(directory = run_one(key), error = NULL),
    error = function(error) {
      if (identical(on_stratum_error, "stop")) stop(error)
      list(directory = NA_character_, error = conditionMessage(error))
    }
  )
  answers <- rc_parallel_lapply(names(groups), run_safe, BPPARAM = BPPARAM)
  status <- do.call(rbind, lapply(seq_along(answers), function(i) {
    data.frame(
      stratum_id = names(groups)[[i]],
      status = if (!is.na(answers[[i]]$directory)) "ok" else "failed",
      output_dir = answers[[i]]$directory,
      error_message = answers[[i]]$error %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  .rc_write_tsv_gz(status, file.path(outdir, "metacell_stratum_status.tsv.gz"))
  directories <- vapply(answers, function(x) x$directory, character(1))
  directories <- directories[!is.na(directories) & nzchar(directories)]
  if (!length(directories)) stop("All SuperCell2 strata failed.", call. = FALSE)
  out <- rc_import_supercell2_metacells(
    directories, rna_assay = rna_assay, atac_assay = atac_assay,
    condition_col = strata_cols[[1L]], celltype_col = label_col %||% "cell_type",
    require_fragments = !identical(fragment_files, FALSE) && !is.null(fragment_files)
  )
  out$stratum_status <- status
  out$strata_cols <- strata_cols
  out$label_col <- label_col
  out$construction_contract <- "regcompass_split_strata_then_supercell2_label"
  out
}

.rc_make_condition_pooled_metacells <- function(
    object, outdir, condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC", fragment_files = FALSE,
    metacell_args = list()) {
  if (!is.list(metacell_args)) stop("`metacell_args` must be a list.", call. = FALSE)
  required <- c(condition_col, celltype_col)
  missing <- setdiff(required, colnames(object@meta.data))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (is.null(metacell_args$gamma)) metacell_args$gamma <- 30L
  cache_contract <- .rc_condition_metacell_cache_contract(
    object, condition_col, celltype_col, rna_assay, atac_assay, metacell_args
  )
  .rc_validate_condition_metacell_cache(
    outdir, cache_contract, overwrite = isTRUE(metacell_args$overwrite %||% FALSE)
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache_contract, file.path(outdir, "condition_metacell_cache_contract.rds"))
  reserved <- intersect(names(metacell_args), c(
    "object", "outdir", "strata_cols", "label_col", "rna_assay",
    "atac_assay", "fragment_files"
  ))
  if (length(reserved)) {
    stop("`metacell_args` cannot override workflow fields: ",
         paste(reserved, collapse = ", "), call. = FALSE)
  }
  defaults <- list(
    object = object, outdir = outdir, strata_cols = condition_col,
    label_col = celltype_col, rna_assay = rna_assay, atac_assay = atac_assay,
    fragment_files = fragment_files, on_stratum_error = "stop"
  )
  defaults[names(metacell_args)] <- NULL
  pooled <- do.call(rc_make_supercell2_metacells, c(defaults, metacell_args))
  pooled <- .rc_assign_metacell_dominant_celltype(pooled, object, celltype_col)
  if (any(pooled$celltype_composition_summary$n_celltypes != 1L) ||
      any(pooled$celltype_composition_summary$mixed_celltype_metacell)) {
    stop("SuperCell2 `label` contract was violated: a metacell mixes cell types.",
         call. = FALSE)
  }
  pooled$metacell_meta$pooling_scope <- "condition_only"
  pooled$metacell_meta$celltype_role <- "supercell2_label_exact"
  pooled$condition_col <- condition_col
  pooled$celltype_col <- celltype_col
  pooled$pooling_scope <- "condition_only"
  pooled$cache_contract <- cache_contract
  pooled$input_design <- list(
    metacell_grouping = condition_col,
    supercell_label_col = celltype_col,
    celltype_assignment = "exact_label_preserving_membership",
    gamma = metacell_args$gamma,
    cache_contract_schema = cache_contract$schema_version,
    inference_policy = paste(
      "RegCompass splits cells by condition before each SuperCell2 call;",
      "cell type is supplied through label and must remain pure within metacells"
    )
  )
  pooled
}
