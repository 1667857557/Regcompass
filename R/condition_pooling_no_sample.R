# Canonical condition-only wrapper loaded after condition_pooling.R.

.rc_condition_only_pool_col <- function(meta) {
  candidate <- ".rc_condition_pool_id"
  while (candidate %in% colnames(meta)) candidate <- paste0(candidate, "_")
  candidate
}

.rc_make_condition_pooled_metacells <- function(
    object, outdir,
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
  required <- unique(c(condition_col, celltype_col))
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
      "Biological-sample balancing is not part of the condition-only workflow: ",
      paste(unsupported, collapse = ", "), call. = FALSE
    )
  }
  reserved <- intersect(names(metacell_args), c(
    "object", "outdir", "condition_col", "celltype_col", "label_col",
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
  if (is.null(metacell_args$gamma)) metacell_args$gamma <- 30L
  cache_contract <- .rc_condition_metacell_cache_contract(
    object = object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    metacell_args = metacell_args
  )
  .rc_validate_condition_metacell_cache(
    outdir = outdir,
    contract = cache_contract,
    overwrite = isTRUE(metacell_args$overwrite %||% FALSE)
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    cache_contract,
    file.path(outdir, "condition_metacell_cache_contract.rds")
  )

  # The generic SuperCell2 adapter predates the condition-only workflow and
  # requires one technical stratum identifier. This private field is generated
  # exclusively from condition and is never interpreted as biological sample
  # metadata or exposed in the public result contract.
  internal_pool_col <- .rc_condition_only_pool_col(object@meta.data)
  object@meta.data[[internal_pool_col]] <- paste0(
    as.character(object@meta.data[[condition_col]]), "__condition_pool"
  )
  internal_celltype_col <- .rc_condition_only_celltype_col(object@meta.data)
  object@meta.data[[internal_celltype_col]] <- "all_celltypes"
  defaults <- list(
    object = object,
    outdir = outdir,
    sample_col = internal_pool_col,
    condition_col = condition_col,
    celltype_col = internal_celltype_col,
    label_col = celltype_col,
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
  meta$condition_pool_id <- meta[[internal_pool_col]]
  meta[[internal_pool_col]] <- NULL
  meta$pooling_scope <- "condition_only"
  meta$celltype_role <- "label_guided_posthoc_dominant_membership"
  pooled$metacell_meta <- meta
  pooled$condition_col <- condition_col
  pooled$celltype_col <- celltype_col
  pooled$internal_celltype_col <- internal_celltype_col
  pooled$analysis_pool_col <- internal_pool_col
  pooled$pooling_scope <- "condition_only"
  pooled$cache_contract <- cache_contract
  pooled$input_design <- list(
    metacell_grouping = condition_col,
    condition_only_stratification = TRUE,
    supercell_label_col = celltype_col,
    celltype_assignment = paste0(
      "SuperCell2 label-guided construction using `", celltype_col,
      "`, followed by dominant membership assignment"
    ),
    ambiguous_celltype_policy = "error_on_tied_dominant_membership",
    gamma = metacell_args$gamma,
    cache_contract_schema = cache_contract$schema_version,
    inference_policy = paste(
      "cells are stratified only by condition; the supplied cell-type label",
      "is passed to SuperCell2 before aggregation to discourage label mixing"
    )
  )
  pooled
}
