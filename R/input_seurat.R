#' Validate an annotated Seurat v4 RNA+ATAC counts object
#'
#' RegCompassR Layer 1 starts from raw counts. The object must contain paired RNA
#' and ATAC assays with identical cell barcodes plus sample and cell-type labels.
#' Optional condition labels define biological strata, but clustering/WNN is not
#' rerun inside RegCompassR.
#' @export
rc_validate_seurat <- function(object,
                               rna_assay = "RNA",
                               atac_assay = "ATAC",
                               sample_col = "sample_id",
                               celltype_col = "cell_type",
                               condition_col = NULL,
                               batch_col = NULL,
                               state_col = NULL,
                               embedding = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class 'Seurat'.", call. = FALSE)
  }

  meta <- object@meta.data
  required_meta <- c(sample_col, celltype_col, condition_col, batch_col, state_col)
  required_meta <- required_meta[!is.na(required_meta) & nzchar(required_meta)]
  missing_meta <- setdiff(required_meta, colnames(meta))
  if (length(missing_meta) > 0L) stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "), call. = FALSE)

  assay_names <- names(object@assays)
  if (!rna_assay %in% assay_names) stop("RNA assay not found: ", rna_assay, call. = FALSE)
  if (!atac_assay %in% assay_names) stop("ATAC assay not found: ", atac_assay, call. = FALSE)

  canonical_cells <- colnames(object)
  rna_cells <- colnames(.rc_get_assay_counts(object, rna_assay))
  atac_cells <- colnames(.rc_get_assay_counts(object, atac_assay))
  if (!setequal(canonical_cells, rna_cells)) stop("RNA assay cell set differs from the Seurat object.", call. = FALSE)
  if (!setequal(canonical_cells, atac_cells)) stop("ATAC assay cell set differs from the Seurat object.", call. = FALSE)

  if (!is.null(embedding) && !embedding %in% names(object@reductions)) {
    stop("Embedding/reduction not found: ", embedding, call. = FALSE)
  }

  invisible(TRUE)
}

#' Extract RNA/ATAC counts and metadata from an annotated Seurat v4 object
#' @export
rc_extract_inputs <- function(object,
                              rna_assay = "RNA",
                              atac_assay = "ATAC",
                              sample_col = "sample_id",
                              celltype_col = "cell_type",
                              condition_col = NULL,
                              batch_col = NULL,
                              state_col = NULL,
                              embedding = NULL) {
  rc_validate_seurat(
    object = object,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    sample_col = sample_col,
    celltype_col = celltype_col,
    condition_col = condition_col,
    batch_col = batch_col,
    state_col = state_col,
    embedding = embedding
  )

  rna <- .rc_get_assay_counts(object, rna_assay)
  atac <- .rc_get_assay_counts(object, atac_assay)
  meta <- object@meta.data

  emb <- NULL
  if (!is.null(embedding)) {
    emb <- object@reductions[[embedding]]@cell.embeddings
  }

  list(rna_counts = rna, atac_counts = atac, rna = rna, atac = atac, meta = meta, embedding = emb)
}

#' Check cell metadata distributions for input diagnostics
#' @export
rc_check_metadata <- function(meta, sample_col = "sample_id", celltype_col = "cell_type", condition_col = NULL, batch_col = NULL, state_col = NULL, min_cells = 30, state_source = NA_character_, state_resolution = NA_character_) {
  required <- c(sample_col, celltype_col, condition_col, batch_col, state_col)
  required <- required[!is.null(required) & !is.na(required)]
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0L) stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  na_counts <- data.frame(column = required, n_na = vapply(required, function(x) sum(is.na(meta[[x]])), integer(1)), stringsAsFactors = FALSE)
  group_cols <- c(sample_col, condition_col, celltype_col, state_col)
  group_cols <- group_cols[!is.null(group_cols) & !is.na(group_cols)]
  cell_counts <- as.data.frame(table(meta[, group_cols, drop = FALSE]), stringsAsFactors = FALSE)
  names(cell_counts)[ncol(cell_counts)] <- "n_cells"
  cell_counts$low_cell_count <- cell_counts$n_cells < min_cells
  condition_batch <- NULL
  if (!is.null(condition_col) && !is.null(batch_col)) condition_batch <- as.data.frame.matrix(table(meta[[condition_col]], meta[[batch_col]]))
  list(cell_counts = cell_counts, na_counts = na_counts, condition_batch = condition_batch,
       state_record = data.frame(state_col = ifelse(is.null(state_col), NA_character_, state_col), state_source = state_source, state_resolution = state_resolution, stringsAsFactors = FALSE))
}

#' Write input metadata summary tables
#' @export
rc_write_input_summary <- function(summary, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(summary$cell_counts, file.path(out_dir, "input_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(summary$na_counts, file.path(out_dir, "input_na_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(summary$condition_batch)) utils::write.table(summary$condition_batch, file.path(out_dir, "condition_batch_table.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  utils::write.table(summary$state_record, file.path(out_dir, "state_source.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(out_dir)
}
#' Recompute metabolic peak-gene links with Signac
#'
#' Runs `Signac::LinkPeaks()` inside RegCompassR for the metabolic genes in a
#' Human-GEM/RegCompass GPR table, then converts the resulting Signac links into
#' the `peak_id`, `gene`, `weight` table consumed by `rc_run_layer1_from_counts()`.
#' Signac is an optional dependency and must be installed by the caller.
rc_recompute_signac_peak_gene_links <- function(object,
                                                gpr_table = NULL,
                                                metabolic_genes = NULL,
                                                peak_assay = "ATAC",
                                                expression_assay = "RNA",
                                                genes_use = NULL,
                                                score_col = c("score", "zscore"),
                                                keep_negative = FALSE,
                                                ...) {
  if (!requireNamespace("Signac", quietly = TRUE)) {
    stop("Package 'Signac' is required to recompute peak-gene links. Install Signac or provide `peak_gene_links` directly.", call. = FALSE)
  }
  if (!requireNamespace("qlcMatrix", quietly = TRUE)) {
    stop("Package 'qlcMatrix' is required by Signac::LinkPeaks(). Install qlcMatrix or provide `peak_gene_links` directly.", call. = FALSE)
  }
  fragments <- Signac::Fragments(object[[peak_assay]])
  if (length(fragments) == 0L) {
    stop("No fragment file is registered for the peak assay; Signac::LinkPeaks() cannot be recomputed. Provide `peak_gene_links` directly or attach fragments to the Seurat object.", call. = FALSE)
  }
  if (is.null(metabolic_genes)) {
    if (is.null(gpr_table)) stop("Provide either `gpr_table` or `metabolic_genes`.", call. = FALSE)
    metabolic_genes <- rc_metabolic_gpr_genes(gpr_table)
  }
  if (!is.null(genes_use)) metabolic_genes <- genes_use
  metabolic_genes <- rc_match_linkpeaks_genes(metabolic_genes, object = object, expression_assay = expression_assay)
  if (length(metabolic_genes) == 0L) stop("No metabolic genes were available for Signac::LinkPeaks().", call. = FALSE)

  object <- Signac::LinkPeaks(
    object = object,
    peak.assay = peak_assay,
    expression.assay = expression_assay,
    genes.use = metabolic_genes,
    ...
  )
  links <- Signac::Links(object[[peak_assay]])
  link_table <- rc_signac_links_to_peak_gene_table(links, score_col = score_col, keep_negative = keep_negative)
  link_table <- rc_filter_peak_gene_links_to_gpr(link_table, tolower(metabolic_genes))
  if (nrow(link_table) == 0L) {
    attr(link_table, "diagnostics") <- data.frame(status = "zero_links", reason = "Signac LinkPeaks returned no metabolic peak-gene links after filtering", n_metacells = ncol(object), stringsAsFactors = FALSE)
  }
  attr(link_table, "seurat_object") <- object
  link_table
}


rc_match_linkpeaks_genes <- function(genes, object, expression_assay = "RNA") {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) return(character(0))
  features <- tryCatch(rownames(object[[expression_assay]]), error = function(e) NULL)
  if (is.null(features) || length(features) == 0L) return(genes)
  feature_key <- stats::setNames(as.character(features), tolower(as.character(features)))
  matched <- unname(feature_key[tolower(genes)])
  matched <- matched[!is.na(matched) & nzchar(matched)]
  unique(matched)
}

#' Convert Signac Links output to a RegCompass peak-gene table
#' @export
rc_signac_links_to_peak_gene_table <- function(links, peak_assay = "ATAC", score_col = c("score", "zscore", "correlation", "cor"), keep_negative = FALSE) {
  score_col <- match.arg(score_col)
  if (inherits(links, "GRanges") || inherits(links, "GenomicRanges")) {
    df <- as.data.frame(links)
  } else if (is.data.frame(links)) {
    df <- links
  } else {
    stop("Unsupported Signac Links object type: ", paste(class(links), collapse = "/"), call. = FALSE)
  }
  if (nrow(df) == 0L) {
    out <- data.frame(peak_id = character(), gene = character(), weight = numeric(), pvalue = numeric(), stringsAsFactors = FALSE)
    attr(out, "peak_assay") <- peak_assay
    return(out)
  }
  gene_col <- intersect(c("gene", "gene_name", "gene.name", "symbol"), colnames(df))[1]
  if (is.na(gene_col)) stop("Signac links must contain a gene column.", call. = FALSE)
  score_candidates <- unique(c(score_col, "score", "zscore", "correlation", "cor"))
  weight_col <- intersect(score_candidates, colnames(df))[1]
  if (is.na(weight_col)) stop("Signac links must contain a score/correlation column for link weights.", call. = FALSE)
  p_col <- intersect(c("pvalue", "p.value", "p_val", "p_val_adj"), colnames(df))[1]
  if ("peak_id" %in% colnames(df)) {
    peak_id <- as.character(df$peak_id)
  } else if ("peak" %in% colnames(df)) {
    peak_id <- as.character(df$peak)
  } else if (any(c("peak.id", "peak_name") %in% colnames(df))) {
    peak_id <- as.character(df[[intersect(c("peak.id", "peak_name"), colnames(df))[1]]])
  } else if (all(c("seqnames", "start", "end") %in% colnames(df))) {
    peak_id <- paste0(df$seqnames, "-", df$start, "-", df$end)
  } else {
    peak_id <- rownames(df)
  }
  out <- data.frame(
    peak_id = peak_id,
    gene = toupper(as.character(df[[gene_col]])),
    weight = as.numeric(df[[weight_col]]),
    pvalue = if (!is.na(p_col)) as.numeric(df[[p_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$peak_id) & nzchar(out$peak_id) & !is.na(out$gene) & nzchar(out$gene) & is.finite(out$weight), , drop = FALSE]
  if (keep_negative) out$weight <- abs(out$weight) else out$weight <- pmax(0, out$weight)
  out <- out[is.finite(out$weight) & out$weight > 0, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "peak_assay") <- peak_assay
  out
}

#' Run Layer 1 from Seurat with internal Signac metabolic relinking by default
#'
#' Extracts RNA/ATAC counts from an annotated Seurat object and, by default,
#' recomputes metabolic peak-gene links internally with `Signac::LinkPeaks()`
#' before calling `rc_run_layer1_from_counts()`.
rc_run_layer1_from_seurat <- function(gpr_table,
                                      object,
                                      unit_map,
                                      unit_meta = NULL,
                                      rna_assay = "RNA",
                                      atac_assay = "ATAC",
                                      sample_col = "sample_id",
                                      celltype_col = "cell_type",
                                      condition_col = NULL,
                                      batch_col = NULL,
                                      state_col = NULL,
                                      embedding = NULL,
                                      peak_gene_links = NULL,
                                      recompute_peak_gene_links = TRUE,
                                      signac_args = list(),
                                      stratum_col = "cell_type",
                                      promiscuity_mode = "sqrt",
                                      and_method = "boltzmann",
                                      tau = 0.20,
                                      or_method = c("sum_sqrtK", "max", "prob_or", "sum"),
                                      bootstrap = FALSE,
                                      low_confidence_threshold = 0.25,
                                      low_confidence_quantile = NULL,
                                      reaction_confidence_method = c("gpr_aware"),
                                      B = 500,
                                      BPPARAM = NULL) {
  or_method <- match.arg(or_method)
  reaction_confidence_method <- match.arg(reaction_confidence_method)
  inputs <- rc_extract_inputs(
    object = object,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    sample_col = sample_col,
    celltype_col = celltype_col,
    condition_col = condition_col,
    batch_col = batch_col,
    state_col = state_col,
    embedding = embedding
  )
  link_source <- "supplied_peak_gene_links"
  signac_object <- NULL
  if (isTRUE(recompute_peak_gene_links)) {
    args <- c(list(
      object = object,
      gpr_table = gpr_table,
      peak_assay = atac_assay,
      expression_assay = rna_assay
    ), signac_args)
    peak_gene_links <- do.call(rc_recompute_signac_peak_gene_links, args)
    signac_object <- attr(peak_gene_links, "seurat_object")
    link_source <- "signac_recomputed_metabolic_links"
  } else if (is.null(peak_gene_links)) {
    link_source <- "none"
  }

  out <- rc_run_layer1_from_counts(
    gpr_table = gpr_table,
    rna_counts = inputs$rna_counts,
    unit_map = unit_map,
    unit_meta = unit_meta,
    atac_counts = inputs$atac_counts,
    peak_gene_links = peak_gene_links,
    stratum_col = stratum_col,
    promiscuity_mode = promiscuity_mode,
    and_method = and_method,
    tau = tau,
    or_method = or_method,
    bootstrap = bootstrap,
    low_confidence_threshold = low_confidence_threshold,
    low_confidence_quantile = low_confidence_quantile,
    reaction_confidence_method = reaction_confidence_method,
    B = B,
    BPPARAM = BPPARAM
  )
  out$peak_gene_link_source <- link_source
  if (!is.null(signac_object)) out$signac_object <- signac_object
  out
}

#' Validate an annotated Seurat/Signac multiome raw-count object
#' @export
rc_validate_multiome_input <- function(object, rna_assay = "RNA", atac_assay = "ATAC", sample_col = "sample_id", condition_col = "condition", celltype_col = "cell_type", state_col = NULL, require_atac = TRUE, require_raw_counts = TRUE) {
  if (!inherits(object, "Seurat")) stop("`object` must inherit from class 'Seurat'.", call. = FALSE)
  meta <- object@meta.data
  required <- c(sample_col, condition_col, celltype_col, state_col); required <- required[!is.null(required) & !is.na(required) & nzchar(required)]
  miss <- setdiff(required, colnames(meta)); if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse = ", "), call. = FALSE)
  assays <- names(object@assays); if (!rna_assay %in% assays) stop("RNA assay not found: ", rna_assay, call. = FALSE); if (require_atac && !atac_assay %in% assays) stop("ATAC assay not found: ", atac_assay, call. = FALSE)
  rna <- .rc_get_assay_counts(object, rna_assay)
  if (require_raw_counts && is.null(rna)) stop("RNA raw counts could not be extracted.", call. = FALSE)
  canonical_cells <- colnames(object)
  if (!setequal(canonical_cells, colnames(rna))) stop("RNA assay cell set differs from the Seurat object.", call. = FALSE)
  rna <- rna[, canonical_cells, drop = FALSE]
  if (require_atac) {
    atac <- .rc_get_assay_counts(object, atac_assay)
    if (!setequal(canonical_cells, colnames(atac))) stop("ATAC assay cell set differs from the Seurat object.", call. = FALSE)
    atac <- atac[, canonical_cells, drop = FALSE]
    stopifnot(identical(colnames(rna), colnames(atac)))
  }
  invisible(TRUE)
}
.rc_get_assay_counts <- function(object, assay) {
  get_assay_data <- SeuratObject::GetAssayData
  args <- list(object = object, assay = assay)
  if ("layer" %in% names(formals(get_assay_data))) {
    args$layer <- "counts"
  } else {
    args[["slot"]] <- "counts"
  }
  do.call(get_assay_data, args)
}
#' @export
rc_get_assay_counts <- function(object, assay) .rc_get_assay_counts(object, assay)
