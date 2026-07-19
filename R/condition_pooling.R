.rc_condition_pool_col <- function(meta) {
  candidate <- ".rc_condition_pool_v170"
  while (candidate %in% colnames(meta)) candidate <- paste0(candidate, "_")
  candidate
}

.rc_condition_group_cols <- function(condition_col, celltype_col) {
  cols <- unique(c(condition_col, celltype_col))
  cols <- cols[!is.na(cols) & nzchar(cols)]
  if (length(cols) != 2L) {
    stop("Condition and cell-type columns must be distinct and non-empty.", call. = FALSE)
  }
  cols
}

.rc_make_condition_pooled_metacells <- function(
    object, outdir,
    sample_col = "sample_id",
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
  required <- c(sample_col, condition_col, celltype_col)
  missing <- setdiff(required, colnames(object@meta.data))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!identical(fragment_files, FALSE) && !is.null(fragment_files)) {
    stop(
      paste(
        "RegCompassR 1.7.0 condition-pooled metacells require",
        "`fragment_files = FALSE` and aggregate the existing ATAC peak-count",
        "assay. Pooling fragment files from multiple samples requires an",
        "explicit per-file barcode map and is not part of the canonical path."
      ),
      call. = FALSE
    )
  }

  object_pool <- object
  pool_col <- .rc_condition_pool_col(object_pool@meta.data)
  object_pool@meta.data[[pool_col]] <- paste0(
    as.character(object_pool@meta.data[[condition_col]]),
    "__condition_pool"
  )

  reserved <- intersect(
    names(metacell_args),
    c(
      "object", "outdir", "sample_col", "condition_col", "celltype_col",
      "rna_assay", "atac_assay", "fragment_files", "save_metacell_object",
      "save_counts", "save_fragments", "require_fragment_aggregation",
      "fragment_aggregation_backend", "on_stratum_error"
    )
  )
  if (length(reserved)) {
    stop(
      "`metacell_args` cannot override condition-pooled workflow fields: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }

  defaults <- list(
    object = object_pool,
    outdir = outdir,
    sample_col = pool_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
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
  pooled <- do.call(rc_make_supercell2_metacells, c(defaults, metacell_args))

  meta <- pooled$metacell_meta
  if (!is.data.frame(meta) || !nrow(meta)) {
    stop("Condition-pooled SuperCell2 produced no metacells.", call. = FALSE)
  }
  if (!all(c(condition_col, celltype_col) %in% colnames(meta))) {
    stop("Condition-pooled metacell metadata are incomplete.", call. = FALSE)
  }
  meta[[sample_col]] <- paste0(as.character(meta[[condition_col]]), "__pooled")
  meta$pooling_scope <- "condition_x_celltype"
  meta$samples_mixed_within_condition <- TRUE
  pooled$metacell_meta <- meta
  pooled$pooling_scope <- "condition_x_celltype"
  pooled$pooled_sample_column <- pool_col
  pooled
}

.rc_normalize_condition_metacell_object <- function(
    pooled, rna_assay = "RNA", atac_assay = "ATAC") {
  object <- rc_load_or_merge_metacell_objects(
    pooled$metacell_objects,
    fragment_manifest = NULL,
    metacell_meta = pooled$metacell_meta,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    require_complete_fragments = FALSE
  )
  object <- Seurat::NormalizeData(object, assay = rna_assay, verbose = FALSE)
  object <- Signac::RunTFIDF(object, assay = atac_assay)
  object
}

.rc_run_condition_pando_modules <- function(
    metacell_object, gem, outdir, pfm, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    pando_args = list(), layer1_args = list()) {
  if (!is.list(pando_args) || !is.list(layer1_args)) {
    stop("`pando_args` and `layer1_args` must be lists.", call. = FALSE)
  }
  group_cols <- .rc_condition_group_cols(condition_col, celltype_col)
  pando_args$group_cols <- NULL
  pando_args$sample_col <- NULL
  pando_args$condition_col <- NULL
  pando_args$celltype_col <- NULL
  pando_args$metacell_object <- NULL
  pando_args$gem <- NULL
  pando_args$outdir <- NULL
  pando_args$pfm <- NULL
  pando_args$genome <- NULL
  pando_args$rna_assay <- NULL
  pando_args$atac_assay <- NULL
  pando_args$BPPARAM <- NULL

  defaults <- list(
    metacell_object = metacell_object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    sample_col = condition_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    group_cols = group_cols,
    single_cell_genes = rownames(.rc_get_assay_counts(metacell_object, rna_assay)),
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    save_sample_metacell_objects = TRUE,
    BPPARAM = FALSE,
    on_sample_error = "stop"
  )
  defaults[names(pando_args)] <- NULL
  meta_modules <- do.call(rc_run_pando_meta_modules, c(defaults, pando_args))
  status <- meta_modules$sample_status
  if (!is.data.frame(status) || !nrow(status) || !"status" %in% colnames(status)) {
    stop("Pando did not return a complete condition-level status table.", call. = FALSE)
  }
  failed <- is.na(status$status) | status$status != "ok"
  if (any(failed)) {
    group_label <- if ("group_id" %in% colnames(status)) {
      as.character(status$group_id[failed])
    } else {
      which(failed)
    }
    stop(
      "Every condition-by-cell-type Pando GRN must complete successfully. Failed: ",
      paste(group_label, collapse = ";"),
      call. = FALSE
    )
  }

  local_fastcore_args <- layer1_args$local_fastcore_args %||% list()
  local_fastcore_args$enabled <- layer1_args$local_fastcore %||%
    local_fastcore_args$enabled %||% TRUE
  completion <- .rc_complete_stratum_meta_modules(
    meta_modules,
    gem,
    outdir = file.path(outdir, "local_fastcore"),
    local_fastcore_args = local_fastcore_args
  )
  meta_modules$biological_reaction_membership <- meta_modules$reaction_membership
  meta_modules$local_completed_reaction_membership <-
    completion$completed_reaction_membership
  meta_modules$local_fastcore_summary <- completion$summary
  meta_modules$local_fastcore_diagnostics <- completion$diagnostics
  meta_modules$local_fastcore_completion_iterations <-
    completion$completion_iterations
  meta_modules$local_fastcore_parent_scope <- completion$parent_scope
  meta_modules
}

.rc_edge_activity_deviation <- function(edge_activity, min_scale = 0.05) {
  edge_activity <- as.matrix(edge_activity)
  centers <- matrixStats::rowMedians(edge_activity, na.rm = TRUE)
  mad_scale <- matrixStats::rowMads(
    edge_activity,
    constant = 1.4826,
    na.rm = TRUE
  )
  iqr_scale <- matrixStats::rowIQRs(edge_activity, na.rm = TRUE) / 1.349
  scale <- pmax(mad_scale, iqr_scale, min_scale, na.rm = TRUE)
  standardized <- sweep(edge_activity, 1L, centers, "-")
  standardized <- sweep(standardized, 1L, scale, "/")
  tanh(standardized)
}

.rc_condition_gene_regulatory_modifier <- function(
    significant_edges, object, unit_meta,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    target_genes = NULL,
    min_scale = 0.05) {
  if (!is.data.frame(significant_edges)) {
    stop("`significant_edges` must be a data.frame.", call. = FALSE)
  }
  required_edges <- c(
    "target", "region", "tf", "estimate", condition_col, celltype_col
  )
  missing_edges <- setdiff(required_edges, colnames(significant_edges))
  if (length(missing_edges)) {
    stop(
      "Pando edge table is missing columns: ",
      paste(missing_edges, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.data.frame(unit_meta) ||
      !all(c("pool_id", condition_col, celltype_col) %in% colnames(unit_meta))) {
    stop("`unit_meta` is incomplete for condition-pooled regulatory scoring.", call. = FALSE)
  }

  units <- colnames(object)
  unit_meta <- unit_meta[match(units, as.character(unit_meta$pool_id)), , drop = FALSE]
  if (anyNA(unit_meta$pool_id)) {
    stop("Metacell metadata do not align to the Pando object.", call. = FALSE)
  }
  genes <- unique(tolower(trimws(as.character(target_genes))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) {
    genes <- unique(tolower(trimws(as.character(significant_edges$target))))
  }
  modifier <- matrix(
    0,
    nrow = length(genes),
    ncol = length(units),
    dimnames = list(genes, units)
  )
  if (!nrow(significant_edges) || !length(genes)) return(modifier)

  edges <- significant_edges
  edges$target <- toupper(trimws(as.character(edges$target)))
  edges$tf <- toupper(trimws(as.character(edges$tf)))
  edges$region <- trimws(as.character(edges$region))
  edges$estimate <- suppressWarnings(as.numeric(edges$estimate))
  edges <- edges[
    !is.na(edges$target) & nzchar(edges$target) &
      !is.na(edges$tf) & nzchar(edges$tf) &
      !is.na(edges$region) & nzchar(edges$region) &
      is.finite(edges$estimate) & edges$estimate != 0,
    , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  atac <- .rc_pando_assay_data(object, atac_assay)
  rna <- .rc_pando_assay_data(object, rna_assay)
  peak_keys <- toupper(.rc_pando_region_key(rownames(atac)))
  peak_keep <- !is.na(peak_keys) & nzchar(peak_keys) & !duplicated(peak_keys)
  peak_lookup <- stats::setNames(rownames(atac)[peak_keep], peak_keys[peak_keep])
  tf_lookup <- .rc_case_insensitive_lookup(rownames(rna))
  edges$.peak_id <- unname(peak_lookup[toupper(.rc_pando_region_key(edges$region))])
  edges$.tf_id <- unname(tf_lookup[edges$tf])
  edges <- edges[
    !is.na(edges$.peak_id) & nzchar(edges$.peak_id) &
      !is.na(edges$.tf_id) & nzchar(edges$.tf_id),
    , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  group_key_edges <- paste(
    as.character(edges[[condition_col]]),
    as.character(edges[[celltype_col]]),
    sep = "\001"
  )
  group_key_units <- paste(
    as.character(unit_meta[[condition_col]]),
    as.character(unit_meta[[celltype_col]]),
    sep = "\001"
  )
  for (group_key in unique(group_key_edges)) {
    group_edges <- edges[group_key_edges == group_key, , drop = FALSE]
    group_units <- units[group_key_units == group_key]
    if (!nrow(group_edges) || !length(group_units)) next

    for (target in unique(group_edges$target)) {
      selected <- group_edges[group_edges$target == target, , drop = FALSE]
      gene_id <- tolower(target)
      if (!gene_id %in% rownames(modifier) || !nrow(selected)) next

      peak_score <- rc_gene_score(
        as.matrix(atac[selected$.peak_id, units, drop = FALSE]),
        mode = "absolute",
        half_saturation = getOption("RegCompassR.atac_half_saturation", 1)
      )
      tf_score <- rc_gene_score(
        as.matrix(rna[selected$.tf_id, units, drop = FALSE]),
        mode = "absolute",
        half_saturation = getOption("RegCompassR.tf_half_saturation", 1)
      )
      edge_deviation <- .rc_edge_activity_deviation(
        peak_score * tf_score,
        min_scale = min_scale
      )

      weight <- abs(selected$estimate)
      weight[!is.finite(weight)] <- 0
      if (!any(weight > 0)) next
      weight <- weight / sum(weight)
      model_rsq <- if ("rsq" %in% colnames(selected)) {
        suppressWarnings(as.numeric(selected$rsq))
      } else {
        numeric()
      }
      model_rsq <- model_rsq[is.finite(model_rsq)]
      reliability <- if (length(model_rsq)) {
        sqrt(min(max(stats::median(model_rsq), 0), 1))
      } else {
        1
      }
      signed_weight <- weight * sign(selected$estimate)
      value <- reliability * as.numeric(crossprod(
        signed_weight,
        edge_deviation[, group_units, drop = FALSE]
      ))
      modifier[gene_id, group_units] <- pmax(pmin(value, 1), -1)
    }
  }
  attr(modifier, "score_semantics") <- paste(
    "condition-specific Pando coefficient sign and magnitude applied to",
    "pooled-reference TF-by-ATAC activity deviations"
  )
  modifier
}
