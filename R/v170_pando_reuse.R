.rc_require_normalized_assay <- function(object, assay, context) {
  value <- .rc_pando_assay_data(object, assay)
  if (!identical(colnames(value), colnames(object))) {
    stop(context, " normalized assay data are not aligned to the Pando units.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Canonical correction: condition x cell-type Pando models reuse normalized
# RNA and the cell-type-shared ATAC TF-IDF from Step 1. No group-specific
# NormalizeData/RunTFIDF pass is performed here. Peaks with zero counts in the
# current Pando group are excluded before motif and GRN inference.
.rc_run_pando_meta_modules_v170 <- function(metacell_object,
                                            gem,
                                            outdir,
                                            pfm,
                                            genome,
                                            sample_col = "sample_id",
                                            condition_col = "condition",
                                            celltype_col = "cell_type",
                                            group_cols = NULL,
                                            single_cell_genes = NULL,
                                            rna_assay = "RNA",
                                            atac_assay = "ATAC",
                                            min_metacells = 20L,
                                            pando_initiate_args = list(exclude_exons = TRUE),
                                            pando_motif_args = list(),
                                            pando_infer_args = list(method = "glm", tf_cor = 0.1, peak_cor = 0,
                                                                    adjust_method = "fdr", parallel = FALSE),
                                            padj_threshold = 0.05,
                                            min_abs_estimate = 0,
                                            min_model_rsq = 0.1,
                                            require_padj = TRUE,
                                            top_k_neighbors = 5L,
                                            min_shared_tfs = 1L,
                                            min_tf_jaccard = 0,
                                            max_targets_per_tf = 200L,
                                            subsystem_table = NULL,
                                            expansion_mode = c("ordered_once", "fixed_point"),
                                            save_sample_metacell_objects = TRUE,
                                            save_pando_objects = TRUE,
                                            BPPARAM = NULL,
                                            on_sample_error = c("record", "stop")) {
  expansion_mode <- match.arg(expansion_mode)
  on_sample_error <- match.arg(on_sample_error)
  if (!inherits(metacell_object, "Seurat")) stop("`metacell_object` must inherit from Seurat.", call. = FALSE)
  if (!requireNamespace("Pando", quietly = TRUE)) stop("Install the pinned Pando fork before running v1.2 meta-modules.", call. = FALSE)
  pando_install <- .rc_validate_pando_repository()
  installed <- pando_install$version
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Signac", quietly = TRUE)) {
    stop("Packages 'Seurat' and 'Signac' are required.", call. = FALSE)
  }
  if (is.null(group_cols)) {
    group_cols <- .rc_strict_stratum_cols(
      sample_col = sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col
    )
  }
  group_cols <- unique(as.character(group_cols))
  missing_group_cols <- setdiff(group_cols, colnames(metacell_object@meta.data))
  if (length(missing_group_cols)) {
    stop("Missing metadata column(s): ", paste(missing_group_cols, collapse = ", "), call. = FALSE)
  }
  if (!sample_col %in% group_cols) stop("`group_cols` must include `sample_col`.", call. = FALSE)
  display_cols <- unique(c("group_id", group_cols))
  module_cols <- unique(c("group_id", group_cols, "sample_id", "module_id"))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "pando_objects"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "sample_metacell_objects"), recursive = TRUE, showWarnings = FALSE)

  metabolic_genes <- gem$metabolic_genes %||% rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(.rc_get_assay_counts(metacell_object, rna_assay))
  if (is.null(single_cell_genes)) single_cell_genes <- rna_genes
  target_upper <- intersect(toupper(.rc_mm_trim_unique(single_cell_genes)),
                            toupper(.rc_mm_trim_unique(metabolic_genes)))
  target_genes <- rna_genes[toupper(rna_genes) %in% target_upper]
  target_genes <- .rc_mm_trim_unique(target_genes)
  if (!length(target_genes)) stop("No overlap between single-cell RNA genes and Human-GEM metabolic genes.", call. = FALSE)

  .rc_require_normalized_assay(metacell_object, rna_assay, "RNA")
  .rc_require_normalized_assay(metacell_object, atac_assay, "ATAC")
  normalization <- metacell_object@misc$regcompass_atac_normalization %||% list()
  if (!identical(normalization$scope, "cell_type_across_conditions")) {
    stop("Pando requires Step 1 ATAC TF-IDF normalized within cell type across conditions.", call. = FALSE)
  }

  meta <- metacell_object@meta.data
  meta$.rc_pando_group_id <- rc_make_stratum_id(meta, group_cols)
  group_ids <- unique(as.character(meta$.rc_pando_group_id))

  run_one_group <- function(group_id) {
    cells <- rownames(meta)[as.character(meta$.rc_pando_group_id) == group_id]
    vals <- meta[match(cells[[1L]], rownames(meta)), group_cols, drop = FALSE]
    sample <- as.character(vals[[sample_col]][[1L]])
    status <- data.frame(group_id = group_id, vals,
      n_metacells = length(cells), n_target_genes = length(target_genes),
      n_atac_peaks_input = NA_integer_, n_zero_count_peaks_excluded = NA_integer_,
      n_atac_peaks_used = NA_integer_,
      status = "pending", n_edges = 0L, n_significant_edges = 0L,
      error_class = NA_character_, error_message = NA_character_,
      stringsAsFactors = FALSE, check.names = FALSE)
    if (length(cells) < as.integer(min_metacells)) {
      status$status <- "skipped_too_few_metacells"
      return(list(status = status, all = data.frame(), significant = data.frame()))
    }
    one <- tryCatch({
      obj <- subset(metacell_object, cells = cells)
      .rc_require_normalized_assay(obj, rna_assay, "RNA")
      .rc_require_normalized_assay(obj, atac_assay, "ATAC")
      filtered <- .rc_drop_zero_count_atac_features(
        obj,
        atac_assay = atac_assay,
        context = paste0("Pando group ", group_id)
      )
      obj <- filtered$object
      .rc_require_normalized_assay(obj, atac_assay, "ATAC")
      group_file_id <- gsub("[^A-Za-z0-9_.-]+", "_", group_id)
      if (isTRUE(save_sample_metacell_objects)) saveRDS(obj, file.path(outdir, "sample_metacell_objects", paste0(group_file_id, ".rds")))
      init_defaults <- list(object = obj, peak_assay = atac_assay, rna_assay = rna_assay)
      init_defaults[names(pando_initiate_args)] <- NULL
      grn <- do.call(Pando::initiate_grn, c(init_defaults, pando_initiate_args))
      motif_defaults <- list(object = grn, pfm = pfm, genome = genome)
      motif_defaults[names(pando_motif_args)] <- NULL
      grn <- do.call(Pando::find_motifs, c(motif_defaults, pando_motif_args))
      infer_defaults <- list(object = grn, genes = target_genes)
      infer_defaults[names(pando_infer_args)] <- NULL
      grn <- do.call(Pando::infer_grn, c(infer_defaults, pando_infer_args))
      tab <- rc_extract_pando_tf_peak_gene(grn, sample_id = sample,
        padj_threshold = padj_threshold, min_abs_estimate = min_abs_estimate,
        min_model_rsq = min_model_rsq, require_padj = require_padj)
      add_group_meta <- function(x) {
        if (!nrow(x)) return(x)
        x$group_id <- group_id
        for (col in group_cols) x[[col]] <- vals[[col]][[1L]]
        x[, c(display_cols, setdiff(colnames(x), display_cols)), drop = FALSE]
      }
      tab$all <- add_group_meta(tab$all)
      tab$significant <- add_group_meta(tab$significant)
      tab$peak_diagnostics <- filtered$diagnostics
      if (isTRUE(save_pando_objects)) saveRDS(grn, file.path(outdir, "pando_objects", paste0(group_file_id, ".rds")))
      tab
    }, error = function(e) e)
    if (inherits(one, "error")) {
      status$status <- "failed"; status$error_class <- class(one)[[1L]]; status$error_message <- conditionMessage(one)
      if (identical(on_sample_error, "stop")) stop(one)
      return(list(status = status, all = data.frame(), significant = data.frame()))
    }
    status$n_atac_peaks_input <- one$peak_diagnostics$n_input_peaks
    status$n_zero_count_peaks_excluded <- one$peak_diagnostics$n_zero_count_peaks_excluded
    status$n_atac_peaks_used <- one$peak_diagnostics$n_retained_peaks
    status$status <- "ok"; status$n_edges <- nrow(one$all); status$n_significant_edges <- nrow(one$significant)
    list(status = status, all = one$all, significant = one$significant)
  }

  group_results <- rc_parallel_lapply(group_ids, run_one_group, BPPARAM = BPPARAM)
  status_table <- do.call(rbind, lapply(group_results, `[[`, "status"))
  all_table <- do.call(rbind, lapply(group_results, `[[`, "all"))
  sig_table <- do.call(rbind, lapply(group_results, `[[`, "significant"))
  if (is.null(all_table)) all_table <- data.frame()
  if (is.null(sig_table)) sig_table <- data.frame()
  .rc_mm_write_tsv_gz(status_table, file.path(outdir, "pando_sample_status.tsv.gz"))
  .rc_mm_write_tsv_gz(all_table, file.path(outdir, "pando_tf_peak_gene_all.tsv.gz"))
  .rc_mm_write_tsv_gz(sig_table, file.path(outdir, "pando_tf_peak_gene_significant.tsv.gz"))
  if (!nrow(sig_table)) stop("No significant Pando TF-peak-gene edges were available across samples.", call. = FALSE)

  sig_for_projection <- sig_table; sig_for_projection$sample_id <- sig_for_projection$group_id
  projection <- rc_project_metabolic_grn(sig_for_projection, metabolic_genes = metabolic_genes,
    top_k = top_k_neighbors, min_shared_tfs = min_shared_tfs,
    min_tf_jaccard = min_tf_jaccard, max_targets_per_tf = max_targets_per_tf,
    include_direct_metabolic_tf = TRUE)
  group_meta <- unique(status_table[, c("group_id", group_cols), drop = FALSE])
  projection$nodes <- .rc_remap_projection_metadata(projection$nodes, group_meta, sample_col, display_cols)
  projection$edges <- .rc_remap_projection_metadata(projection$edges, group_meta, sample_col, display_cols)
  core <- rc_map_meta_module_core_reactions(projection$nodes, gem$gpr_table)
  if (nrow(core)) {
    core <- merge(core, unique(projection$nodes[, module_cols, drop = FALSE]), by = c("sample_id", "module_id"), all.x = TRUE, sort = FALSE)
    core <- core[, c(display_cols, setdiff(colnames(core), display_cols)), drop = FALSE]
  }
  expanded <- rc_expand_meta_module_reactions(gem, core, subsystem_table = subsystem_table, expansion_mode = expansion_mode)
  if (nrow(expanded$reaction_membership)) {
    expanded$reaction_membership <- merge(expanded$reaction_membership, unique(core[, module_cols, drop = FALSE]), by = c("sample_id", "module_id"), all.x = TRUE, sort = FALSE)
    expanded$reaction_membership <- expanded$reaction_membership[, c(display_cols, setdiff(colnames(expanded$reaction_membership), display_cols)), drop = FALSE]
  }
  if (nrow(expanded$summary)) {
    expanded$summary <- merge(expanded$summary, unique(core[, module_cols, drop = FALSE]), by = c("sample_id", "module_id"), all.x = TRUE, sort = FALSE)
    expanded$summary <- expanded$summary[, c(display_cols, setdiff(colnames(expanded$summary), display_cols)), drop = FALSE]
  }
  .rc_mm_write_tsv_gz(projection$nodes, file.path(outdir, "metabolic_gene_nodes.tsv.gz"))
  .rc_mm_write_tsv_gz(projection$edges, file.path(outdir, "metabolic_gene_edges.tsv.gz"))
  .rc_mm_write_tsv_gz(core, file.path(outdir, "core_gene_reaction.tsv.gz"))
  .rc_mm_write_tsv_gz(expanded$reaction_membership, file.path(outdir, "meta_module_reactions.tsv.gz"))
  .rc_mm_write_tsv_gz(expanded$summary, file.path(outdir, "meta_module_summary.tsv.gz"))

  out <- list(schema_version = "regcompass_pando_meta_module_v1.2",
    pando_installed_version = installed, pando_remote_username = pando_install$remote_username,
    pando_remote_repo = pando_install$remote_repo, pando_remote_ref = pando_install$remote_ref,
    pando_remote_sha = pando_install$remote_sha, target_metabolic_genes = target_genes,
    sample_status = status_table, tf_peak_gene_all = all_table,
    tf_peak_gene_significant = sig_table, metabolic_gene_nodes = projection$nodes,
    metabolic_gene_edges = projection$edges, core_gene_reaction = core,
    reaction_membership = expanded$reaction_membership, meta_module_summary = expanded$summary,
    crossref_maps = expanded$crossref_maps,
    normalization_policy = list(
      rna = "Step 1 normalized RNA data reused without group renormalization",
      atac = "cell-type-shared TF-IDF reused without group renormalization",
      zero_count_peaks = "excluded globally before TF-IDF and again within each Pando group"))
  saveRDS(out, file.path(outdir, "pando_meta_modules.rds"))
  out
}
