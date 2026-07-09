#' Validate metacell-level RegCompass inputs
#' @export
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
    if (!identical(colnames(rna_metacell_counts), colnames(atac_metacell_counts))) stop("RNA and ATAC metacell count matrices must have identical colnames in the same order.", call. = FALSE)
  }
  invisible(TRUE)
}

.rc_as_sparse <- function(x) {
  if (inherits(x, "sparseMatrix")) return(x)
  methods::as(x, "dgCMatrix")
}

.rc_metacell_meta_for_pool_apis <- function(metacell_meta, metacell_id_col = "metacell_id") {
  out <- metacell_meta
  if (!"pool_id" %in% colnames(out)) out$pool_id <- as.character(out[[metacell_id_col]])
  if (!"unit_id" %in% colnames(out)) out$unit_id <- out$pool_id
  out
}

#' Aggregate raw counts by a cell-to-metacell membership table
.rc_aggregate_counts_by_membership <- function(counts, membership, fun = c("sum", "mean"), metacell_id_col = "metacell_id", cell_id_col = "cell_id", BPPARAM = NULL) {
  fun <- match.arg(fun)
  if (!is.data.frame(membership)) stop("`membership` must be a data.frame.", call. = FALSE)
  missing <- setdiff(c(metacell_id_col, cell_id_col), colnames(membership))
  if (length(missing) > 0L) stop("`membership` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  map <- data.frame(pool_id = as.character(membership[[metacell_id_col]]), cell_id = as.character(membership[[cell_id_col]]), stringsAsFactors = FALSE)
  rc_pseudobulk_counts(counts, map, fun = fun, BPPARAM = BPPARAM)
}

#' Filter empty metacells before normalization
#' @export
rc_filter_empty_metacells <- function(metacell_counts, metacell_meta, metacell_id_col = "metacell_id") {
  lib <- Matrix::colSums(metacell_counts)
  keep <- lib > 0
  if (any(!keep)) warning(sum(!keep), " empty metacells removed before normalization", call. = FALSE)
  list(counts = metacell_counts[, keep, drop = FALSE], metacell_meta = metacell_meta[match(colnames(metacell_counts)[keep], metacell_meta[[metacell_id_col]]), , drop = FALSE])
}

#' Build one-row-per-metacell metadata from membership
#' @export
rc_build_metacell_metadata <- function(membership, metacell_id_col = "metacell_id", cell_id_col = "cell_id") {
  if (!is.data.frame(membership)) stop("`membership` must be a data.frame.", call. = FALSE)
  missing <- setdiff(c(metacell_id_col, cell_id_col), colnames(membership))
  if (length(missing) > 0L) stop("`membership` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  x <- membership[!is.na(membership[[metacell_id_col]]), , drop = FALSE]
  cols <- setdiff(colnames(x), cell_id_col)
  out <- x[!duplicated(x[[metacell_id_col]]), cols, drop = FALSE]
  out$n_cells <- as.integer(tabulate(match(x[[metacell_id_col]], out[[metacell_id_col]]), nbins = nrow(out)))
  rownames(out) <- NULL
  out
}

#' Compute metacell-level detection from raw metacell counts
#' @export
rc_metacell_detection <- function(metacell_counts) {
  methods::as(metacell_counts > 0, "dgCMatrix")
}

#' Filter ATAC peaks detected in metacells and compute metacell logCPM
#' @export
rc_atac_metacell_logcpm <- function(atac_metacell_counts, min_metacells = 3) {
  detected <- Matrix::rowSums(atac_metacell_counts > 0) >= min_metacells
  rc_logcpm(atac_metacell_counts[detected, , drop = FALSE])
}

#' Run Layer 1 from SuperCell2.0 metacell raw counts
#' @export
rc_run_layer1_from_metacells <- function(gpr_table,
                                         rna_metacell_counts,
                                         metacell_meta,
                                         atac_metacell_counts,
                                         metacell_seurat,
                                         peak_gene_links = NULL,
                                         allow_supplied_links = FALSE,
                                         force_metacell_relink = TRUE,
                                         link_stratum_cols = "cell_type",
                                         min_metacells_for_linkpeaks = 10,
                                         metabolic_genes = NULL,
                                         linkpeaks_args = list(),
                                         stratum_col = "cell_type",
                                         promiscuity_mode = "sqrt",
                                         and_method = "boltzmann",
                                         tau = 0.20,
                                         reaction_confidence_method = "gpr_aware",
                                         bootstrap = TRUE,
                                         B = 500,
                                         BPPARAM = NULL) {
  if (!is.null(peak_gene_links) && !isTRUE(allow_supplied_links)) {
    stop("`peak_gene_links` supplied by users are not accepted in formal multiome mode. RegCompassR recomputes metabolic peak-gene links on metacells.", call. = FALSE)
  }
  if (is.null(atac_metacell_counts)) stop("`atac_metacell_counts` is required in formal multiome mode.", call. = FALSE)
  if (isTRUE(force_metacell_relink)) {
    if (is.null(metacell_seurat)) stop("`metacell_seurat` is required for metacell-level LinkPeaks.", call. = FALSE)
    peak_gene_links <- do.call(rc_recompute_metacell_peak_gene_links_by_stratum, c(list(
      metacell_object = metacell_seurat,
      metacell_meta = metacell_meta,
      gpr_table = gpr_table,
      metabolic_genes = metabolic_genes,
      link_stratum_cols = link_stratum_cols,
      min_metacells_for_linkpeaks = min_metacells_for_linkpeaks
    ), linkpeaks_args))
  }
  if (is.null(peak_gene_links) || nrow(peak_gene_links) == 0L) stop("Metacell-level metabolic peak-gene relinking produced no usable links.", call. = FALSE)
  rc_validate_metacell_inputs(rna_metacell_counts, metacell_meta, atac_metacell_counts = atac_metacell_counts)
  rna_metacell_counts <- .rc_as_sparse(rna_metacell_counts)
  filtered <- rc_filter_empty_metacells(rna_metacell_counts, metacell_meta)
  rna_metacell_counts <- filtered$counts
  metacell_meta <- filtered$metacell_meta
  rna_logcpm <- rc_logcpm(rna_metacell_counts)
  rna_detection <- rc_metacell_detection(rna_metacell_counts)
  atac_metacell_counts <- atac_metacell_counts[, colnames(rna_logcpm), drop = FALSE]
  atac_peak <- rc_atac_metacell_logcpm(atac_metacell_counts, min_metacells = 3)
  pool_meta <- .rc_metacell_meta_for_pool_apis(metacell_meta)
  p_rna <- rc_percentile_by_stratum(rna_logcpm, pool_meta = pool_meta, stratum_col = stratum_col)
  p_atac_peak <- rc_percentile_by_stratum(atac_peak, pool_meta = pool_meta, stratum_col = stratum_col)
  link_conf <- rc_link_confidence_by_stratum(p_atac_peak = p_atac_peak, peak_gene_links = peak_gene_links, pool_meta = pool_meta, link_stratum_cols = link_stratum_cols)
  genes <- intersect(rownames(p_rna), rownames(link_conf))
  if (length(genes) == 0L) stop("No overlap between metacell RNA genes and linked metabolic genes.", call. = FALSE)
  gene_conf <- rc_concordance_null_correct(p_rna[genes, , drop = FALSE], link_conf[genes, , drop = FALSE], pool_meta = pool_meta, stratum_col = stratum_col)
  out <- rc_run_layer1_capacity(gpr_table = gpr_table, pool_expression = rna_logcpm, pool_detection = rna_detection, pool_meta = pool_meta, stratum_col = stratum_col, gene_confidence = gene_conf, promiscuity_mode = promiscuity_mode, and_method = and_method, tau = tau, reaction_confidence_method = reaction_confidence_method, bootstrap = bootstrap, B = B, BPPARAM = BPPARAM)
  out$metacell_meta <- metacell_meta
  out$pool_meta <- .rc_metacell_meta_for_pool_apis(out$metacell_meta)
  out$rna_metacell_logcpm <- rna_logcpm
  out$rna_metacell_detection <- rna_detection
  out$metacell_peak_gene_links <- peak_gene_links
  out$peak_gene_link_source <- "cell_type_stratified_metacell_recomputed_metabolic_links"
  out$layer1_unit <- "metacell"
  out
}


rc_prepare_metacell_linkpeaks_object <- function(object,
                                                peak_assay = "ATAC",
                                                expression_assay = "RNA",
                                                normalize_expression = TRUE,
                                                run_region_stats = TRUE,
                                                genome = NULL,
                                                genome_package = "BSgenome.Hsapiens.UCSC.hg38",
                                                annotation_package = "EnsDb.Hsapiens.v86") {
  rc_load_linkpeaks_reference_packages(genome_package = genome_package, annotation_package = annotation_package)
  old_assay <- tryCatch(SeuratObject::DefaultAssay(object), error = function(e) NULL)

  if (isTRUE(normalize_expression)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required to run NormalizeData on metacell RNA before LinkPeaks.", call. = FALSE)
    object <- rc_set_default_assay(object, expression_assay)
    object <- Seurat::NormalizeData(object = object, assay = expression_assay, verbose = FALSE)
  }

  if (isTRUE(run_region_stats) && !rc_has_region_stats(object, peak_assay = peak_assay)) {
    genome <- rc_resolve_linkpeaks_genome(genome = genome, genome_package = genome_package)
    object <- rc_set_default_assay(object, peak_assay)
    object <- Signac::RegionStats(object = object, assay = peak_assay, genome = genome, verbose = FALSE)
  }
  if (!is.null(old_assay)) object <- rc_set_default_assay(object, old_assay)
  object
}

rc_load_linkpeaks_reference_packages <- function(genome_package = "BSgenome.Hsapiens.UCSC.hg38",
                                                annotation_package = "EnsDb.Hsapiens.v86") {
  for (pkg in c(annotation_package, genome_package)) {
    if (!is.null(pkg) && nzchar(pkg) && !requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for default metacell LinkPeaks preprocessing. Install it, or pass a precomputed `genome`/precomputed RegionStats and set the package argument to NULL.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

rc_resolve_linkpeaks_genome <- function(genome = NULL, genome_package = "BSgenome.Hsapiens.UCSC.hg38") {
  if (!is.null(genome)) return(genome)
  if (is.null(genome_package) || !nzchar(genome_package)) {
    stop("Metacell LinkPeaks requires peak-level region statistics. Provide `genome` (for example BSgenome.Hsapiens.UCSC.hg38) so Signac::RegionStats() can run, or precompute RegionStats on the metacell ATAC assay.", call. = FALSE)
  }
  if (!requireNamespace(genome_package, quietly = TRUE)) {
    stop("Package '", genome_package, "' is required to compute default Signac::RegionStats() for metacell LinkPeaks.", call. = FALSE)
  }
  getExportedValue(genome_package, genome_package)
}

rc_set_default_assay <- function(object, assay) {
  setter <- get("DefaultAssay<-", envir = asNamespace("SeuratObject"))
  setter(object, value = assay)
}

rc_has_region_stats <- function(object, peak_assay = "ATAC") {
  meta_features <- tryCatch(object[[peak_assay]]@meta.features, error = function(e) NULL)
  if (is.null(meta_features)) return(FALSE)
  all(c("GC.percent", "sequence.length") %in% colnames(meta_features))
}

#' Recompute metabolic peak-gene links on a metacell Signac object
#' @export
rc_recompute_metacell_peak_gene_links <- function(metacell_object,
                                                  metabolic_genes = NULL,
                                                  peak_assay = "ATAC",
                                                  expression_assay = "RNA",
                                                  distance = 5e5,
                                                  min_cells = 5,
                                                  out_file = NULL,
                                                  gpr_table = NULL,
                                                  require_fragments = TRUE,
                                                  normalize_expression = TRUE,
                                                  run_region_stats = TRUE,
                                                  genome = NULL,
                                                  genome_package = "BSgenome.Hsapiens.UCSC.hg38",
                                                  annotation_package = "EnsDb.Hsapiens.v86",
                                                  ...) {
  if (!inherits(metacell_object, "Seurat")) stop("`metacell_object` must be a metacell-level Seurat/Signac object.", call. = FALSE)
  if (!requireNamespace("Signac", quietly = TRUE)) stop("Package 'Signac' is required for metacell-level LinkPeaks.", call. = FALSE)
  frags <- Signac::Fragments(metacell_object[[peak_assay]])
  if (require_fragments && length(frags) == 0L) stop("Metacell-level LinkPeaks requires fragment files registered on the metacell ATAC assay. Run metacell fragment aggregation successfully before Layer 1 multiome analysis.", call. = FALSE)
  metacell_object <- rc_prepare_metacell_linkpeaks_object(metacell_object, peak_assay = peak_assay, expression_assay = expression_assay, normalize_expression = normalize_expression, run_region_stats = run_region_stats, genome = genome, genome_package = genome_package, annotation_package = annotation_package)
  links <- rc_recompute_signac_peak_gene_links(object = metacell_object, gpr_table = gpr_table, metabolic_genes = metabolic_genes, peak_assay = peak_assay, expression_assay = expression_assay, distance = distance, min.cells = min_cells, ...)
  if (nrow(links) == 0L) stop("Metacell-level metabolic peak-gene relinking returned 0 links.", call. = FALSE)
  if (!is.null(out_file)) .rc_write_tsv_gz(links, out_file)
  links
}

#' Recompute metacell peak-gene links separately within metadata strata
#' @export
rc_recompute_metacell_peak_gene_links_by_stratum <- function(metacell_object,
                                                             metacell_meta,
                                                             gpr_table,
                                                             metabolic_genes = NULL,
                                                             link_stratum_cols = "cell_type",
                                                             min_metacells_for_linkpeaks = 10,
                                                             on_too_few_metacells = c("skip", "pool_by_cell_type", "global", "stop"),
                                                             diagnostics_file = NULL,
                                                             peak_assay = "ATAC",
                                                             expression_assay = "RNA",
                                                             ...) {
  on_too_few_metacells <- match.arg(on_too_few_metacells)
  if (!inherits(metacell_object, "Seurat")) stop("`metacell_object` must be a metacell-level Seurat object.", call. = FALSE)
  if (!"metacell_id" %in% colnames(metacell_meta)) stop("`metacell_meta` must contain `metacell_id`.", call. = FALSE)
  missing_cols <- setdiff(link_stratum_cols, colnames(metacell_meta))
  if (length(missing_cols) > 0L) stop("Missing link stratum columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  if (is.null(metabolic_genes)) metabolic_genes <- rc_metabolic_gpr_genes(gpr_table)
  metabolic_genes <- unique(as.character(metabolic_genes))
  metabolic_genes <- metabolic_genes[!is.na(metabolic_genes) & nzchar(metabolic_genes)]
  if (length(metabolic_genes) == 0L) stop("No metabolic/GPR genes available for metacell LinkPeaks.", call. = FALSE)
  stratum_id <- interaction(metacell_meta[, link_stratum_cols, drop = FALSE], sep = "|", drop = TRUE)
  metacell_meta$link_stratum <- as.character(stratum_id)
  strata <- split(as.character(metacell_meta$metacell_id), metacell_meta$link_stratum)
  diagnostics <- list()
  links <- lapply(names(strata), function(st) {
    cells <- intersect(strata[[st]], colnames(metacell_object))
    if (length(cells) < min_metacells_for_linkpeaks) {
      reason <- paste0("too_few_metacells: ", length(cells), " < ", min_metacells_for_linkpeaks)
      if (identical(on_too_few_metacells, "stop")) stop("Too few metacells for LinkPeaks in stratum `", st, "`: ", length(cells), " < ", min_metacells_for_linkpeaks, call. = FALSE)
      diagnostics[[st]] <<- data.frame(link_stratum = st, n_metacells = length(cells), n_links = 0L, status = "skipped", reason = reason, stringsAsFactors = FALSE)
      return(NULL)
    }
    obj_st <- subset(metacell_object, cells = cells)
    x <- tryCatch(
      rc_recompute_metacell_peak_gene_links(metacell_object = obj_st, gpr_table = gpr_table, metabolic_genes = metabolic_genes, peak_assay = peak_assay, expression_assay = expression_assay, ...),
      error = function(e) {
        if (identical(on_too_few_metacells, "stop")) stop(e)
        diagnostics[[st]] <<- data.frame(link_stratum = st, n_metacells = length(cells), n_links = 0L, status = "failed", reason = conditionMessage(e), stringsAsFactors = FALSE)
        NULL
      }
    )
    if (is.null(x)) return(NULL)
    x$link_stratum <- st
    diagnostics[[st]] <<- data.frame(link_stratum = st, n_metacells = length(cells), n_links = nrow(x), status = "ok", reason = NA_character_, stringsAsFactors = FALSE)
    x
  })
  links <- links[!vapply(links, is.null, logical(1))]
  diag <- if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame(link_stratum = character(), n_metacells = integer(), n_links = integer(), status = character(), reason = character())
  rownames(diag) <- NULL
  if (!is.null(diagnostics_file)) .rc_write_tsv_gz(diag, diagnostics_file)
  if (!length(links)) {
    out <- data.frame()
    attr(out, "diagnostics") <- diag
    return(out)
  }
  out <- do.call(rbind, links)
  rownames(out) <- NULL
  attr(out, "diagnostics") <- diag
  out
}

#' Sample-level long summary with metacell diagnostics
#' @export
rc_metacell_sample_summary <- function(score_mat, metacell_meta, sample_col = "sample_id", celltype_col = "cell_type", condition_col = NULL) {
  score_mat <- as.matrix(score_mat)
  required <- c("metacell_id", sample_col, celltype_col)
  missing <- setdiff(required, colnames(metacell_meta))
  if (length(missing) > 0L) stop("`metacell_meta` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  metacell_meta <- metacell_meta[metacell_meta$metacell_id %in% colnames(score_mat), , drop = FALSE]
  group_cols <- c(sample_col, condition_col, celltype_col)
  group_cols <- group_cols[!is.null(group_cols) & !is.na(group_cols)]
  keys <- interaction(metacell_meta[, group_cols, drop = FALSE], drop = TRUE, sep = "|")
  pieces <- lapply(split(metacell_meta, keys), function(mm) {
    vals <- score_mat[, mm$metacell_id, drop = FALSE]
    data.frame(group_id = unique(keys[match(mm$metacell_id, metacell_meta$metacell_id)])[1], reaction_id = rownames(score_mat), value = matrixStats::rowMedians(vals, na.rm = TRUE), n_metacells_used = ncol(vals), n_cells_used = if ("n_cells" %in% colnames(mm)) sum(mm$n_cells, na.rm = TRUE) else NA_real_, low_power_group_flag = if ("low_power_metacell" %in% colnames(mm)) any(mm$low_power_metacell, na.rm = TRUE) else NA, single_metacell_group_flag = ncol(vals) == 1L, stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}

#' Minimal metacell diagnostics for the Layer 1 workflow
#' @export
rc_metacell_diagnostics <- function(metacell_meta, rna_metacell_counts = NULL, atac_metacell_counts = NULL, gpr_genes = NULL) {
  if (!is.data.frame(metacell_meta) || !"metacell_id" %in% colnames(metacell_meta)) stop("`metacell_meta` must contain `metacell_id`.", call. = FALSE)
  ids <- as.character(metacell_meta$metacell_id)
  out <- metacell_meta
  out$RNA_depth <- if (is.null(rna_metacell_counts)) NA_real_ else as.numeric(Matrix::colSums(rna_metacell_counts[, ids, drop = FALSE]))
  out$RNA_detected_genes <- if (is.null(rna_metacell_counts)) NA_real_ else as.numeric(Matrix::colSums(rna_metacell_counts[, ids, drop = FALSE] > 0))
  out$ATAC_detected_peaks <- if (is.null(atac_metacell_counts)) NA_real_ else as.numeric(Matrix::colSums(atac_metacell_counts[, ids, drop = FALSE] > 0))
  gpr_genes <- rc_match_matrix_features(gpr_genes, rna_metacell_counts)
  out$GPR_gene_detection_rate <- if (is.null(rna_metacell_counts) || length(gpr_genes) == 0L) NA_real_ else as.numeric(Matrix::colMeans(rna_metacell_counts[gpr_genes, ids, drop = FALSE] > 0))
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
  rc_get_assay_counts(object, assay)
}

.rc_extract_supercell_membership <- function(mc_object, original_cells, metacell_ids) {
  meta <- mc_object@meta.data
  candidates <- c("cell_membership", "membership", "cells", "cell_ids", "single_cell_ids", "SC")
  for (nm in intersect(candidates, colnames(meta))) {
    vals <- meta[[nm]]
    if (is.list(vals)) {
      return(do.call(rbind, lapply(seq_along(vals), function(i) data.frame(cell_id = as.character(vals[[i]]), metacell_id = metacell_ids[[i]], stringsAsFactors = FALSE))))
    }
  }
  misc_mem <- tryCatch(mc_object@misc$membership, error = function(e) NULL)
  if (!is.null(misc_mem)) {
    if (is.data.frame(misc_mem)) return(misc_mem)
    if (is.numeric(misc_mem) || is.integer(misc_mem)) {
      cell_ids <- names(misc_mem)
      if (is.null(cell_ids) || any(!nzchar(cell_ids))) cell_ids <- original_cells
      idx <- as.integer(misc_mem)
      mc <- metacell_ids[idx]
      ok <- !is.na(mc) & nzchar(mc)
      return(data.frame(cell_id = as.character(cell_ids[ok]), metacell_id = as.character(mc[ok]), stringsAsFactors = FALSE))
    }
  }
  attr_map <- attr(mc_object, "membership")
  if (!is.null(attr_map) && is.data.frame(attr_map)) return(attr_map)
  data.frame(cell_id = character(0), metacell_id = character(0), stringsAsFactors = FALSE)
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
#' @export
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
                                         overwrite = FALSE,
                                         BPPARAM = NULL) {
  if (!requireNamespace("SuperCell", quietly = TRUE)) stop("Package 'SuperCell' is required for rc_make_supercell2_metacells(). By default, install SuperCell from https://github.com/1667857557/SuperCell-Seurat-V4/tree/supercell-2.0, or import existing metacells with rc_import_supercell2_metacells().", call. = FALSE)
  if (!inherits(object, "Seurat")) stop("`object` must inherit from class 'Seurat'.", call. = FALSE)
  if (is.null(fragment_files)) fragment_files <- .rc_fragment_files_from_atac(object, atac_assay = atac_assay)
  if (isTRUE(require_fragment_aggregation)) {
    if (!isTRUE(save_fragments)) stop("Formal multiome workflow requires `save_fragments = TRUE`.", call. = FALSE)
    if (is.null(fragment_files)) stop("Formal multiome workflow requires `fragment_files` for metacell fragment aggregation, or a fragment file registered on the ATAC assay.", call. = FALSE)
  }
  fragment_files <- .rc_normalize_fragment_files(fragment_files, atac_assay = atac_assay)
  .rc_assert_shell_safe_paths(outdir, unlist(fragment_files, use.names = FALSE))
  meta <- object@meta.data
  required <- c(sample_col, condition_col, celltype_col, state_col, label_col)
  required <- required[!is.null(required) & !is.na(required) & nzchar(required)]
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0L) stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  group_cols <- c(sample_col, condition_col, celltype_col, state_col)
  group_cols <- group_cols[!is.null(group_cols) & !is.na(group_cols) & nzchar(group_cols)]
  meta <- rc_drop_na_grouping(meta, group_cols)
  meta$cell_id <- rownames(meta)
  keys <- interaction(meta[, group_cols, drop = FALSE], drop = TRUE, sep = "|")
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
    min_required_cells <- max(as.integer(min_cells_per_stratum), as.integer(min_metacells_per_stratum) * as.integer(gamma))
    gamma_i <- as.integer(gamma)
    if (length(cells) < min_required_cells) {
      if (isTRUE(adaptive_gamma) && length(cells) >= as.integer(min_metacells_per_stratum) * as.integer(min_metacell_size)) {
        gamma_i <- max(1L, floor(length(cells) / as.integer(min_metacells_per_stratum)))
      } else {
        reason <- if (length(cells) < min_cells_per_stratum) "stratum_below_min_cells_per_stratum" else "stratum_too_small_for_min_metacells_at_current_gamma"
        diag <- data.frame(group_id = key, n_cells = length(cells), skipped = TRUE, skip_reason = reason, gamma = gamma, min_required_cells = min_required_cells, stringsAsFactors = FALSE)
        .rc_write_tsv_gz(diag, file.path(stratum_dir, "qc", "metacell_qc.tsv.gz"))
        return(stratum_dir)
      }
    }
    prefix <- paste(vapply(vals, .rc_safe_path_component, character(1)), collapse = "_")
    seu_sub <- subset(object, cells = cells)
    args <- list(seurat = seu_sub, assay = c(rna_assay, atac_assay), reduction = list(rna_reduction, atac_reduction), dims = list(rna_dims, atac_dims), gamma = gamma_i, return.seurat = TRUE, prefixMC = "Metacell_")
    if (!is.null(label_col)) args$label <- label_col
    if (save_fragments && !is.null(fragment_files)) {
      args$fragmentFiles <- fragment_files
      args$outputDirMcFragment <- file.path(stratum_dir, "fragments")
      args$bgzip_path <- bgzip_path
      args$tabix_path <- tabix_path
      args$nb_cl <- max(1L, as.integer(fragment_nb_cl))
    }
    mc <- tryCatch(
      do.call(SuperCell::SCimplify_for_Seurat, args),
      error = function(e) {
        if (isTRUE(require_fragment_aggregation)) stop("Metacell fragment aggregation failed; metacell-level LinkPeaks cannot be recomputed: ", conditionMessage(e), call. = FALSE)
        warning("Fragment aggregation failed; continuing only because `require_fragment_aggregation = FALSE`.", call. = FALSE)
        args2 <- args[setdiff(names(args), c("fragmentFiles", "outputDirMcFragment", "bgzip_path", "tabix_path"))]
        do.call(SuperCell::SCimplify_for_Seurat, args2)
      }
    )
    mc_ids <- colnames(.rc_get_assay_counts_safe(mc, rna_assay))
    display_ids <- paste0(prefix, "_MC", sprintf(paste0("%0", max(3, nchar(length(mc_ids))), "d"), seq_along(mc_ids)))
    membership <- .rc_extract_supercell_membership(mc, cells, mc_ids)
    if (nrow(membership) == 0L) warning("Could not infer single-cell membership from SuperCell output for stratum ", key, "; membership.tsv.gz will be empty.", call. = FALSE)
    for (col in group_cols) membership[[col]] <- vals[[col]][[1]]
    mc_meta <- rc_build_metacell_metadata(membership)
    if (nrow(mc_meta) == 0L) mc_meta <- data.frame(metacell_id = mc_ids, n_cells = NA_integer_, stringsAsFactors = FALSE)
    for (col in group_cols) if (!col %in% colnames(mc_meta)) mc_meta[[col]] <- vals[[col]][[1]]
    mc_meta$metacell_display_id <- display_ids[match(as.character(mc_meta$metacell_id), mc_ids)]
    mc_meta$low_power_metacell <- !is.na(mc_meta$n_cells) & mc_meta$n_cells < min_metacell_size
    rna_counts <- .rc_as_sparse(.rc_get_assay_counts_safe(mc, rna_assay))
    atac_counts <- .rc_as_sparse(.rc_get_assay_counts_safe(mc, atac_assay))
    if (save_metacell_object) saveRDS(mc, file.path(stratum_dir, "metacell_object.rds"))
    .rc_write_tsv_gz(membership, file.path(stratum_dir, "membership.tsv.gz"))
    .rc_write_tsv_gz(mc_meta, file.path(stratum_dir, "metacell_metadata.tsv.gz"))
    if (save_counts) {
      saveRDS(rna_counts, file.path(stratum_dir, "rna_counts.rds"))
      saveRDS(atac_counts, file.path(stratum_dir, "atac_counts.rds"))
    }
    diag <- data.frame(group_id = key, n_cells = length(cells), n_metacells = length(mc_ids), gamma = gamma_i, requested_gamma = gamma, min_metacell_size = min_metacell_size, skipped = FALSE, stringsAsFactors = FALSE)
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
        min_cells_per_stratum = min_cells_per_stratum,
        min_metacell_size = min_metacell_size,
        min_metacells_per_stratum = min_metacells_per_stratum,
        adaptive_gamma = adaptive_gamma,
        label_col = label_col,
        fragment_files = fragment_files,
        bgzip_path = bgzip_path,
        tabix_path = tabix_path,
        fragment_nb_cl = max(1L, as.integer(fragment_nb_cl)),
        save_metacell_object = save_metacell_object,
        save_counts = save_counts,
        save_fragments = save_fragments,
        require_fragment_aggregation = require_fragment_aggregation
      )
      yaml::write_yaml(run_params, file.path(stratum_dir, "qc", "run_params.yaml"))
    }
    stratum_dir
  }
  dirs <- unlist(rc_parallel_lapply(names(groups), run_one, BPPARAM = BPPARAM), use.names = FALSE)
  rc_import_supercell2_metacells(dirs, rna_assay = rna_assay, atac_assay = atac_assay, sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col, require_fragments = require_fragment_aggregation)
}

#' Import saved SuperCell2.0 metacell outputs
#' @export
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
  metas <- memberships <- list(); rnas <- atacs <- list(); objects <- fragments <- character(0)
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
    frag <- Sys.glob(file.path(d, "fragments", "*.tsv.gz"))
    if (length(frag)) fragments <- c(fragments, frag)
  }
  metacell_meta <- do.call(rbind, metas); rownames(metacell_meta) <- NULL
  membership <- if (length(memberships)) do.call(rbind, memberships) else data.frame()
  if (length(rnas) == 0L) stop("No rna_counts.rds files were found in metacell directories.", call. = FALSE)
  rna_counts <- do.call(cbind, lapply(rnas, .rc_as_sparse))
  atac_counts <- if (length(atacs)) do.call(cbind, lapply(atacs, .rc_as_sparse)) else NULL
  rc_validate_metacell_inputs(rna_counts, metacell_meta, atac_metacell_counts = atac_counts, sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col)
  if (require_fragments) {
    missing_idx <- fragments[!file.exists(paste0(fragments, ".tbi"))]
    if (length(fragments) == 0L || length(missing_idx) > 0L) stop("Metacell fragment files or tabix indexes are missing.", call. = FALSE)
  }
  list(metacell_meta = metacell_meta, membership = membership, rna_counts = rna_counts, atac_counts = atac_counts, metacell_objects = objects, fragment_files = fragments, diagnostics = data.frame(n_metacells = ncol(rna_counts), n_membership_rows = nrow(membership)))
}

#' Write a metacell QC report
#' @export
rc_write_metacell_report <- function(metacell_meta, file, diagnostics = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  n_cells <- if ("n_cells" %in% colnames(metacell_meta)) metacell_meta$n_cells else NA_real_
  lines <- c("# RegCompassR metacell report", "", paste0("- Metacells: ", nrow(metacell_meta)), paste0("- Median cells per metacell: ", stats::median(n_cells, na.rm = TRUE)), paste0("- Low-power metacell fraction: ", if ("low_power_metacell" %in% colnames(metacell_meta)) mean(metacell_meta$low_power_metacell, na.rm = TRUE) else NA_real_))
  if (!is.null(diagnostics)) lines <- c(lines, "", "## Diagnostics", utils::capture.output(print(utils::head(diagnostics, 20L))))
  writeLines(lines, con = file)
  invisible(file)
}

#' Build sample-aware RNA+ATAC metacells
#' @export
rc_make_metacells <- function(..., allow_empty_membership = FALSE, filter_low_power_metacells = TRUE, min_metacell_size = 20, require_fragment_aggregation = TRUE) {
  out <- rc_make_supercell2_metacells(..., min_metacell_size = min_metacell_size, require_fragment_aggregation = require_fragment_aggregation)
  if (nrow(out$membership) == 0L && !allow_empty_membership) stop("SuperCell membership could not be extracted. Check SuperCell version/output schema.", call. = FALSE)
  if ("n_cells" %in% colnames(out$metacell_meta)) {
    out$metacell_meta_all <- out$metacell_meta
    out$metacell_meta$low_power_metacell <- out$metacell_meta$n_cells < min_metacell_size
    out$metacell_meta_used <- if (filter_low_power_metacells) out$metacell_meta[!out$metacell_meta$low_power_metacell, , drop=FALSE] else out$metacell_meta
    out$rna_counts_all <- out$rna_counts
    if (!is.null(out$atac_counts)) out$atac_counts_all <- out$atac_counts
    if (filter_low_power_metacells) {
      ids <- as.character(out$metacell_meta_used$metacell_id)
      out$rna_counts <- out$rna_counts[, ids, drop = FALSE]
      if (!is.null(out$atac_counts)) out$atac_counts <- out$atac_counts[, ids, drop = FALSE]
    }
    out$low_power_metacell_fraction <- mean(out$metacell_meta$low_power_metacell, na.rm=TRUE)
  }
  out
}
#' @export
rc_import_metacells <- function(metacell_dirs, ..., filter_low_power_metacells = TRUE, min_metacell_size = 20) {
  out <- rc_import_supercell2_metacells(metacell_dirs, ...)
  out$metacell_meta_all <- out$metacell_meta
  if ("n_cells" %in% colnames(out$metacell_meta)) out$metacell_meta$low_power_metacell <- out$metacell_meta$n_cells < min_metacell_size
  keep <- if (filter_low_power_metacells && "low_power_metacell" %in% colnames(out$metacell_meta)) !out$metacell_meta$low_power_metacell else rep(TRUE, nrow(out$metacell_meta))
  out$metacell_meta_used <- out$metacell_meta[keep,,drop=FALSE]
  ids <- out$metacell_meta_used$metacell_id; out$rna_counts <- out$rna_counts[, ids, drop=FALSE]; if (!is.null(out$atac_counts)) out$atac_counts <- out$atac_counts[, ids, drop=FALSE]
  out$low_power_metacell_fraction <- if ("low_power_metacell" %in% colnames(out$metacell_meta)) mean(out$metacell_meta$low_power_metacell, na.rm=TRUE) else NA_real_
  out
}


.rc_seurat4_filterobjects <- function(object, classes.keep = c("Assay", "Assay5", "ChromatinAssay")) {
  assays <- names(object@assays)
  assays[vapply(assays, function(a) {
    any(vapply(classes.keep, function(cl) inherits(object@assays[[a]], cl), logical(1)))
  }, logical(1))]
}

.rc_install_seurat4_filterobjects <- function(envir = .GlobalEnv) invisible(TRUE)

.rc_with_seurat4_filterobjects <- function(expr) force(expr)
#' Load and merge saved metacell Seurat objects
#' @export
rc_load_or_merge_metacell_objects <- function(metacell_objects, fragment_files = NULL, rna_assay = "RNA", atac_assay = "ATAC") {
  if (is.null(metacell_objects) || length(metacell_objects) == 0L) stop("No metacell Seurat objects supplied.", call. = FALSE)
  objs <- lapply(metacell_objects, function(x) if (inherits(x, "Seurat")) x else readRDS(x))
  cells_by_fragment <- lapply(objs, colnames)
  obj <- if (length(objs) == 1L) objs[[1L]] else Reduce(function(a, b) merge(a, y = b), objs)
  .rc_register_signac_fragments(obj, fragment_files = fragment_files, cells_by_fragment = cells_by_fragment, atac_assay = atac_assay, replace_existing = TRUE)
}

.rc_register_signac_fragments <- function(object, fragment_files = NULL, cells_by_fragment = NULL, atac_assay = "ATAC", replace_existing = TRUE) {
  if (is.null(fragment_files) || length(fragment_files) == 0L) return(object)
  fragment_files <- as.character(fragment_files)
  if (is.null(cells_by_fragment)) cells_by_fragment <- rep(list(colnames(object)), length(fragment_files))
  if (length(cells_by_fragment) != length(fragment_files)) stop("`cells_by_fragment` must have one cell vector per fragment file.", call. = FALSE)
  missing <- fragment_files[!file.exists(fragment_files)]
  missing_index <- fragment_files[!file.exists(paste0(fragment_files, ".tbi"))]
  if (length(missing) > 0L) stop("Metacell fragment files are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (length(missing_index) > 0L) stop("Metacell fragment tabix indexes are missing: ", paste(paste0(missing_index, ".tbi"), collapse = ", "), call. = FALSE)
  if (!requireNamespace("Signac", quietly = TRUE)) stop("Package 'Signac' is required to register metacell fragment files.", call. = FALSE)
  if (!atac_assay %in% names(object@assays)) stop("Metacell object is missing ATAC assay `", atac_assay, "`.", call. = FALSE)
  frag_setter <- get("Fragments<-", envir = asNamespace("Signac"))
  if (isTRUE(replace_existing)) object[[atac_assay]] <- frag_setter(object[[atac_assay]], value = list())
  fragments <- Map(function(path, cells) {
    tryCatch(
      Signac::CreateFragmentObject(path = path, cells = as.character(cells), validate.fragments = FALSE),
      error = function(e) stop("Failed to register metacell fragment file `", path, "`: ", conditionMessage(e), call. = FALSE)
    )
  }, fragment_files, cells_by_fragment)
  object[[atac_assay]] <- frag_setter(object[[atac_assay]], value = fragments)
  object
}

#' Run the formal sample-aware metacell multiome workflow
#' @export
rc_run_regcompass_multiome_metacell <- function(object, gpr_table, outdir, fragment_files = NULL, sample_col = "sample_id", condition_col = "condition", celltype_col = "cell_type", state_col = NULL, label_col = NULL, rna_assay = "RNA", atac_assay = "ATAC", rna_reduction = "pca", atac_reduction = "lsi", rna_dims = 1:30, atac_dims = 2:30, gamma = 100, min_cells_per_stratum = 100, min_metacell_size = 20, min_metacells_per_stratum = 2L, adaptive_gamma = FALSE, fragment_nb_cl = 1L, require_fragment_aggregation = TRUE, save_fragments = TRUE, save_metacell_object = TRUE, save_counts = TRUE, overwrite = FALSE, BPPARAM_metacell = FALSE, link_stratum_cols = "cell_type", min_metacells_for_linkpeaks = 10, linkpeaks_args = list(), layer1_args = list(), future_plan = c("sequential", "current"), future_globals_max_size = 8 * 1024^3) {
  future_plan <- match.arg(future_plan)
  if (future_plan == "sequential" && requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    old_max_size <- getOption("future.globals.maxSize")
    on.exit(future::plan(old_plan), add = TRUE)
    on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)
    future::plan(future::sequential)
    options(future.globals.maxSize = future_globals_max_size)
  }
  mc <- rc_make_metacells(object = object, outdir = file.path(outdir, "01_metacells"), sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col, state_col = state_col, label_col = label_col, rna_assay = rna_assay, atac_assay = atac_assay, rna_reduction = rna_reduction, atac_reduction = atac_reduction, rna_dims = rna_dims, atac_dims = atac_dims, gamma = gamma, min_cells_per_stratum = min_cells_per_stratum, min_metacell_size = min_metacell_size, min_metacells_per_stratum = min_metacells_per_stratum, adaptive_gamma = adaptive_gamma, fragment_files = fragment_files, fragment_nb_cl = fragment_nb_cl, save_fragments = save_fragments, save_metacell_object = save_metacell_object, save_counts = save_counts, require_fragment_aggregation = require_fragment_aggregation, overwrite = overwrite, BPPARAM = BPPARAM_metacell)
  metacell_seurat <- rc_load_or_merge_metacell_objects(mc$metacell_objects, fragment_files = mc$fragment_files, rna_assay = rna_assay, atac_assay = atac_assay)
  do.call(rc_run_layer1_from_metacells, c(list(gpr_table = gpr_table, rna_metacell_counts = mc$rna_counts, atac_metacell_counts = mc$atac_counts, metacell_meta = mc$metacell_meta, metacell_seurat = metacell_seurat, force_metacell_relink = TRUE, allow_supplied_links = FALSE, link_stratum_cols = link_stratum_cols, min_metacells_for_linkpeaks = min_metacells_for_linkpeaks, linkpeaks_args = linkpeaks_args), layer1_args))
}
