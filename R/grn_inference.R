.rc_require_normalized_assay <- function(object, assay, context) {
  value <- .rc_pando_assay_data(object, assay)
  if (!identical(colnames(value), colnames(object))) stop(context, " normalized assay data are not aligned to the analysis units.", call. = FALSE)
  invisible(TRUE)
}

.rc_run_condition_single_cell_grns <- function(
    object, gem, outdir, pfm, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list(method = "glm", tf_cor = 0.1, peak_cor = 0.01, adjust_method = "fdr", parallel = FALSE),
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = TRUE,
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    on_group_error = c("record", "stop")) {
  on_group_error <- match.arg(on_group_error)
  if (!is.list(pando_infer_args)) stop("`pando_infer_args` must be a list.", call. = FALSE)
  pando_infer_args <- modifyList(
    list(method = "glm", tf_cor = 0.1, peak_cor = 0.01, adjust_method = "fdr", parallel = FALSE),
    pando_infer_args
  )
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  if (!requireNamespace("Pando", quietly = TRUE)) stop("Install the pinned Pando fork before running GRN inference.", call. = FALSE)
  pando_install <- .rc_validate_pando_repository()
  group_cols <- c(condition_col, celltype_col)
  missing <- setdiff(group_cols, colnames(object@meta.data))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  .rc_require_normalized_assay(object, rna_assay, "RNA")
  .rc_require_normalized_assay(object, atac_assay, "ATAC")
  normalization <- object@misc$regcompass_atac_normalization %||% list()
  if (!identical(normalization$scope, "cell_type_across_conditions")) stop("Pando requires cell-type-shared ATAC TF-IDF across conditions.", call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "pando_objects"), recursive = TRUE, showWarnings = FALSE)
  metabolic_genes <- gem$metabolic_genes %||% rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  target_upper <- intersect(toupper(.rc_mm_trim_unique(rna_genes)), toupper(.rc_mm_trim_unique(metabolic_genes)))
  target_genes <- .rc_mm_trim_unique(rna_genes[toupper(rna_genes) %in% target_upper])
  if (!length(target_genes)) stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)
  meta <- object@meta.data
  meta$.rc_pando_group_id <- rc_make_stratum_id(meta, group_cols)
  group_ids <- unique(as.character(meta$.rc_pando_group_id))
  run_one_group <- function(group_id) {
    cells <- rownames(meta)[as.character(meta$.rc_pando_group_id) == group_id]
    vals <- meta[match(cells[[1L]], rownames(meta)), group_cols, drop = FALSE]
    status <- data.frame(
      group_id = group_id,
      condition = as.character(vals[[condition_col]][[1L]]),
      cell_type = as.character(vals[[celltype_col]][[1L]]),
      n_cells = length(cells), n_target_genes = length(target_genes),
      n_atac_peaks_input = NA_integer_, n_zero_count_peaks_excluded = NA_integer_, n_atac_peaks_used = NA_integer_,
      status = "pending", n_edges = 0L, n_significant_edges = 0L,
      error_class = NA_character_, error_message = NA_character_, stringsAsFactors = FALSE
    )
    names(status)[2:3] <- group_cols
    if (length(cells) < as.integer(min_cells)) {
      status$status <- "skipped_too_few_cells"
      return(list(status = status, all = data.frame(), significant = data.frame()))
    }
    one <- tryCatch({
      obj <- subset(object, cells = cells)
      filtered <- .rc_drop_zero_count_atac_features(obj, atac_assay, paste0("Pando group ", group_id))
      obj <- filtered$object
      init_defaults <- list(object = obj, peak_assay = atac_assay, rna_assay = rna_assay)
      init_defaults[names(pando_initiate_args)] <- NULL
      grn <- do.call(Pando::initiate_grn, c(init_defaults, pando_initiate_args))
      motif_defaults <- list(object = grn, pfm = pfm, genome = genome)
      motif_defaults[names(pando_motif_args)] <- NULL
      grn <- do.call(Pando::find_motifs, c(motif_defaults, pando_motif_args))
      infer_defaults <- list(object = grn, genes = target_genes)
      infer_defaults[names(pando_infer_args)] <- NULL
      grn <- do.call(Pando::infer_grn, c(infer_defaults, pando_infer_args))
      tab <- rc_extract_pando_tf_peak_gene(
        grn, sample_id = group_id, padj_threshold = padj_threshold,
        min_abs_estimate = min_abs_estimate, min_model_rsq = min_model_rsq,
        require_padj = require_padj
      )
      if (nrow(tab$significant)) {
        reliable_rsq <- if ("rsq" %in% colnames(tab$significant)) {
          value <- suppressWarnings(as.numeric(tab$significant$rsq))
          is.finite(value) & value >= min_model_rsq
        } else {
          rep(FALSE, nrow(tab$significant))
        }
        tab$significant <- tab$significant[reliable_rsq, , drop = FALSE]
      }
      add_meta <- function(x) {
        if (!nrow(x)) return(x)
        x$group_id <- group_id
        for (col in group_cols) x[[col]] <- as.character(vals[[col]][[1L]])
        x[, c("group_id", group_cols, setdiff(colnames(x), c("group_id", group_cols))), drop = FALSE]
      }
      tab$all <- add_meta(tab$all)
      tab$significant <- add_meta(tab$significant)
      tab$peak_diagnostics <- filtered$diagnostics
      if (isTRUE(save_pando_objects)) saveRDS(grn, file.path(outdir, "pando_objects", paste0(gsub("[^A-Za-z0-9_.-]+", "_", group_id), ".rds")))
      tab
    }, error = function(e) e)
    if (inherits(one, "error")) {
      status$status <- "failed"; status$error_class <- class(one)[[1L]]; status$error_message <- conditionMessage(one)
      if (identical(on_group_error, "stop")) stop(one)
      return(list(status = status, all = data.frame(), significant = data.frame()))
    }
    status$n_atac_peaks_input <- one$peak_diagnostics$n_input_peaks
    status$n_zero_count_peaks_excluded <- one$peak_diagnostics$n_zero_count_peaks_excluded
    status$n_atac_peaks_used <- one$peak_diagnostics$n_retained_peaks
    status$status <- "ok"; status$n_edges <- nrow(one$all); status$n_significant_edges <- nrow(one$significant)
    list(status = status, all = one$all, significant = one$significant)
  }
  results <- rc_parallel_lapply(group_ids, run_one_group, BPPARAM = BPPARAM)
  status <- do.call(rbind, lapply(results, `[[`, "status"))
  all_edges <- do.call(rbind, lapply(results, `[[`, "all")); if (is.null(all_edges)) all_edges <- data.frame()
  significant <- do.call(rbind, lapply(results, `[[`, "significant")); if (is.null(significant)) significant <- data.frame()
  .rc_mm_write_tsv_gz(status, file.path(outdir, "pando_group_status.tsv.gz"))
  .rc_mm_write_tsv_gz(all_edges, file.path(outdir, "pando_tf_peak_gene_all.tsv.gz"))
  .rc_mm_write_tsv_gz(significant, file.path(outdir, "pando_tf_peak_gene_significant.tsv.gz"))
  failed <- status$status != "ok"
  if (any(failed)) stop("Every condition-by-cell-type Pando GRN must complete successfully. Failed: ", paste(status$group_id[failed], collapse = "; "), call. = FALSE)
  if (!nrow(significant)) stop("No significant Pando TF-peak-gene edges were available.", call. = FALSE)
  answer <- list(
    schema_version = "regcompass_single_cell_grn_v1",
    pando_installed_version = pando_install$version,
    target_metabolic_genes = target_genes,
    sample_status = status,
    tf_peak_gene_all = all_edges,
    tf_peak_gene_significant = significant,
    normalization_policy = list(
      rna = "global single-cell NormalizeData before condition splitting",
      atac = "cell-type-shared TF-IDF across conditions before condition splitting",
      zero_count_peaks = "excluded globally before TF-IDF and within each Pando group",
      pando_peak_cor = pando_infer_args$peak_cor %||% 0.01,
      pando_rsq = paste0("finite rsq >= ", min_model_rsq)
    ),
    group_cols = group_cols
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_build_condition_meta_modules <- function(grn_result, gem, outdir, layer1_args = list()) {
  if (!is.list(grn_result) || !is.data.frame(grn_result$tf_peak_gene_significant)) stop("`grn_result` is not a valid single-cell GRN result.", call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  group_cols <- grn_result$group_cols
  display_cols <- c("group_id", group_cols)
  module_cols <- unique(c(display_cols, "sample_id", "module_id"))
  metabolic_genes <- gem$metabolic_genes %||% rc_metabolic_gpr_genes(gem$gpr_table)
  sig <- grn_result$tf_peak_gene_significant
  sig$sample_id <- sig$group_id
  projection <- rc_project_metabolic_grn(
    sig, metabolic_genes = metabolic_genes,
    top_k = layer1_args$top_k_neighbors %||% 5L,
    min_shared_tfs = layer1_args$min_shared_tfs %||% 1L,
    min_tf_jaccard = layer1_args$min_tf_jaccard %||% 0,
    max_targets_per_tf = layer1_args$max_targets_per_tf %||% 200L,
    include_direct_metabolic_tf = TRUE
  )
  group_meta <- unique(grn_result$sample_status[, display_cols, drop = FALSE])
  group_meta$analysis_unit_id <- group_meta$group_id
  projection$nodes <- .rc_remap_projection_metadata(projection$nodes, group_meta, "analysis_unit_id", display_cols)
  projection$edges <- .rc_remap_projection_metadata(projection$edges, group_meta, "analysis_unit_id", display_cols)
  core <- rc_map_meta_module_core_reactions(projection$nodes, gem$gpr_table)
  if (nrow(core)) {
    core <- merge(core, unique(projection$nodes[, module_cols, drop = FALSE]), by = c("sample_id", "module_id"), all.x = TRUE, sort = FALSE)
    core <- core[, c(display_cols, setdiff(colnames(core), display_cols)), drop = FALSE]
  }
  expanded <- rc_expand_meta_module_reactions(
    gem, core,
    subsystem_table = layer1_args$subsystem_table %||% NULL,
    expansion_mode = layer1_args$expansion_mode %||% "ordered_once"
  )
  if (nrow(expanded$reaction_membership)) {
    expanded$reaction_membership <- merge(expanded$reaction_membership, unique(core[, module_cols, drop = FALSE]), by = c("sample_id", "module_id"), all.x = TRUE, sort = FALSE)
    expanded$reaction_membership <- expanded$reaction_membership[, c(display_cols, setdiff(colnames(expanded$reaction_membership), display_cols)), drop = FALSE]
  }
  out <- c(grn_result, list(
    metabolic_gene_nodes = projection$nodes,
    metabolic_gene_edges = projection$edges,
    core_gene_reaction = core,
    reaction_membership = expanded$reaction_membership,
    meta_module_summary = expanded$summary,
    crossref_maps = expanded$crossref_maps
  ))
  local_args <- layer1_args$local_fastcore_args %||% list()
  local_args$enabled <- layer1_args$local_fastcore %||% local_args$enabled %||% TRUE
  completion <- .rc_complete_stratum_meta_modules(out, gem, outdir = file.path(outdir, "local_fastcore"), local_fastcore_args = local_args)
  out$biological_reaction_membership <- out$reaction_membership
  out$local_completed_reaction_membership <- completion$completed_reaction_membership
  out$local_fastcore_summary <- completion$summary
  out$local_fastcore_diagnostics <- completion$diagnostics
  out$local_fastcore_completion_iterations <- completion$completion_iterations
  out$local_fastcore_parent_scope <- completion$parent_scope
  out$analysis_group_unit <- "condition_x_celltype_single_cell_grn"
  .rc_mm_write_tsv_gz(projection$nodes, file.path(outdir, "metabolic_gene_nodes.tsv.gz"))
  .rc_mm_write_tsv_gz(projection$edges, file.path(outdir, "metabolic_gene_edges.tsv.gz"))
  .rc_mm_write_tsv_gz(core, file.path(outdir, "core_gene_reaction.tsv.gz"))
  .rc_mm_write_tsv_gz(expanded$reaction_membership, file.path(outdir, "meta_module_reactions.tsv.gz"))
  saveRDS(out, file.path(outdir, "condition_meta_modules.rds"))
  out
}
