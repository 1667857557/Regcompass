.rc_condition_only_sample_col <- function(meta) {
  candidate <- ".rc_condition_pool_id"
  while (candidate %in% colnames(meta)) candidate <- paste0(candidate, "_")
  candidate
}

.rc_condition_only_celltype_col <- function(meta) {
  candidate <- ".rc_condition_only_celltype"
  while (candidate %in% colnames(meta)) candidate <- paste0(candidate, "_")
  candidate
}

.rc_assign_metacell_dominant_celltype <- function(
    pooled, object, celltype_col) {
  membership <- pooled$membership
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership))) {
    stop("Condition-only metacells require cell-to-metacell membership.",
         call. = FALSE)
  }
  cell_index <- match(as.character(membership$cell_id), rownames(object@meta.data))
  if (anyNA(cell_index)) {
    stop("Metacell membership contains cells absent from the input object.",
         call. = FALSE)
  }
  cell_type <- trimws(as.character(object@meta.data[[celltype_col]][cell_index]))
  if (anyNA(cell_type) || any(!nzchar(cell_type))) {
    stop("Cell-type metadata are incomplete in metacell membership.",
         call. = FALSE)
  }
  membership[[celltype_col]] <- cell_type
  composition <- stats::aggregate(
    rep.int(1L, nrow(membership)),
    by = list(
      metacell_id = as.character(membership$metacell_id),
      cell_type = cell_type
    ),
    FUN = sum
  )
  colnames(composition) <- c("metacell_id", celltype_col, "n_cells")
  split_rows <- split(seq_len(nrow(composition)), composition$metacell_id)
  summary <- do.call(rbind, lapply(names(split_rows), function(id) {
    z <- composition[split_rows[[id]], , drop = FALSE]
    ord <- order(-z$n_cells, as.character(z[[celltype_col]]))
    z <- z[ord, , drop = FALSE]
    total <- sum(z$n_cells)
    top_count <- max(z$n_cells)
    tied <- sum(z$n_cells == top_count) > 1L
    data.frame(
      metacell_id = id,
      dominant_celltype = if (tied) NA_character_ else
        as.character(z[[celltype_col]][[1L]]),
      dominant_celltype_fraction = top_count / total,
      n_celltypes = nrow(z),
      mixed_celltype_metacell = nrow(z) > 1L,
      dominant_celltype_tied = tied,
      stringsAsFactors = FALSE
    )
  }))
  tied_ids <- as.character(summary$metacell_id[
    summary$dominant_celltype_tied %in% TRUE
  ])
  if (length(tied_ids)) {
    stop(
      "Condition-only SuperCell2 produced metacells with tied dominant cell types: ",
      paste(utils::head(tied_ids, 10L), collapse = ", "),
      ". A condition-by-cell-type GRN cannot be assigned unambiguously.",
      call. = FALSE
    )
  }
  meta <- pooled$metacell_meta
  index <- match(as.character(meta$metacell_id), summary$metacell_id)
  if (anyNA(index)) {
    stop("Dominant cell-type summaries do not align with metacell metadata.",
         call. = FALSE)
  }
  meta[[celltype_col]] <- summary$dominant_celltype[index]
  meta$dominant_celltype_fraction <-
    summary$dominant_celltype_fraction[index]
  meta$n_celltypes <- summary$n_celltypes[index]
  meta$mixed_celltype_metacell <- summary$mixed_celltype_metacell[index]
  meta$dominant_celltype_tied <- summary$dominant_celltype_tied[index]
  pooled$membership <- membership
  pooled$metacell_meta <- meta
  pooled$celltype_composition <- composition
  pooled$celltype_composition_summary <- summary
  pooled
}

.rc_make_condition_pooled_metacells <- function(
    object, outdir,
    sample_col = NULL,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list()) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
  }
  required <- c(condition_col, celltype_col)
  missing <- setdiff(required, colnames(object@meta.data))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invalid <- vapply(
    object@meta.data[, required, drop = FALSE],
    function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))),
    logical(1)
  )
  if (any(invalid)) {
    stop("Condition and cell-type metadata must be complete.", call. = FALSE)
  }
  if (!identical(fragment_files, FALSE) && !is.null(fragment_files)) {
    stop(
      paste(
        "The canonical condition-only path requires `fragment_files = FALSE`",
        "and aggregates the existing ATAC peak-count assay."
      ),
      call. = FALSE
    )
  }
  unsupported <- intersect(
    names(metacell_args), c("sample_balance", "sample_balance_seed")
  )
  if (length(unsupported)) {
    stop(
      "Sample balancing is not part of the canonical workflow: ",
      paste(unsupported, collapse = ", "), call. = FALSE
    )
  }
  if (is.null(metacell_args$gamma)) metacell_args$gamma <- 75L
  internal_sample_col <- .rc_condition_only_sample_col(object@meta.data)
  object@meta.data[[internal_sample_col]] <- paste0(
    as.character(object@meta.data[[condition_col]]), "__condition_pool"
  )
  internal_celltype_col <- .rc_condition_only_celltype_col(object@meta.data)
  object@meta.data[[internal_celltype_col]] <- "all_celltypes"
  reserved <- intersect(names(metacell_args), c(
    "object", "outdir", "sample_col", "condition_col", "celltype_col",
    "rna_assay", "atac_assay", "fragment_files", "save_metacell_object",
    "save_counts", "save_fragments", "require_fragment_aggregation",
    "fragment_aggregation_backend", "on_stratum_error"
  ))
  if (length(reserved)) {
    stop(
      "`metacell_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "), call. = FALSE
    )
  }
  defaults <- list(
    object = object,
    outdir = outdir,
    sample_col = internal_sample_col,
    condition_col = condition_col,
    celltype_col = internal_celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = FALSE,
    save_metacell_object = TRUE,
    save_counts = TRUE,
    save_fragments = FALSE,
    require_fragment_aggregation = FALSE,
    fragment_aggregation_backend = "none",
    on_stratum_error = "stop"
  )
  defaults[names(metacell_args)] <- NULL
  pooled <- do.call(
    rc_make_supercell2_metacells,
    c(defaults, metacell_args)
  )
  pooled <- .rc_assign_metacell_dominant_celltype(
    pooled = pooled,
    object = object,
    celltype_col = celltype_col
  )
  meta <- pooled$metacell_meta
  if (!is.data.frame(meta) || !nrow(meta)) {
    stop("Condition-only SuperCell2 produced no metacells.", call. = FALSE)
  }
  meta$pooled_sample_id <- meta[[internal_sample_col]]
  meta$pooling_scope <- "condition_only"
  meta$sample_weighting <- "none"
  meta$sample_col_role <- "internal_condition_pool_id"
  meta$celltype_role <- "posthoc_dominant_membership_label"
  pooled$metacell_meta <- meta
  pooled$input_sample_col <- sample_col
  pooled$analysis_sample_col <- internal_sample_col
  pooled$condition_col <- condition_col
  pooled$celltype_col <- celltype_col
  pooled$internal_celltype_col <- internal_celltype_col
  pooled$pooling_scope <- "condition_only"
  pooled$sample_weighting <- "none"
  pooled$input_design <- list(
    metacell_grouping = condition_col,
    condition_only_stratification = TRUE,
    celltype_assignment = "dominant membership after condition-only SuperCell2",
    ambiguous_celltype_policy = "error_on_tied_dominant_membership",
    gamma = metacell_args$gamma,
    inference_policy = paste(
      "cells are stratified only by condition; sample and cell-type metadata",
      "are not used for selection, weighting or metacell grouping"
    )
  )
  pooled
}
