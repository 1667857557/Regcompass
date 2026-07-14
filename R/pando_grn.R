#' Extract and filter a Pando TF-peak-gene coefficient table
#' @export
rc_extract_pando_tf_peak_gene <- function(grn_object,
                                          sample_id,
                                          padj_threshold = 0.05,
                                          min_abs_estimate = 0,
                                          min_model_rsq = 0.1,
                                          require_padj = TRUE) {
  if (!requireNamespace("Pando", quietly = TRUE)) stop("Package 'Pando' is required.", call. = FALSE)
  coefs <- as.data.frame(stats::coef(grn_object), stringsAsFactors = FALSE)
  if (!nrow(coefs)) {
    empty <- data.frame(sample_id = character(), tf = character(), target = character(), region = character(), stringsAsFactors = FALSE)
    return(list(all = empty, significant = empty))
  }
  required <- c("tf", "target", "region")
  missing <- setdiff(required, colnames(coefs))
  if (length(missing)) stop("Pando coefficient table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  fit <- tryCatch(as.data.frame(Pando::gof(grn_object), stringsAsFactors = FALSE), error = function(e) data.frame())
  if (nrow(fit) && "target" %in% colnames(fit)) {
    keep_fit <- setdiff(colnames(fit), intersect(colnames(fit), setdiff(colnames(coefs), "target")))
    coefs <- merge(coefs, fit[, keep_fit, drop = FALSE], by = "target", all.x = TRUE, sort = FALSE)
  }
  coefs$sample_id <- as.character(sample_id)
  coefs$tf <- toupper(as.character(coefs$tf))
  coefs$target <- toupper(as.character(coefs$target))
  coefs$region <- as.character(coefs$region)
  coefs <- coefs[, c("sample_id", setdiff(colnames(coefs), "sample_id")), drop = FALSE]

  keep <- rep(TRUE, nrow(coefs))
  if ("estimate" %in% colnames(coefs)) keep <- keep & is.finite(coefs$estimate) & abs(coefs$estimate) >= min_abs_estimate
  if ("rsq" %in% colnames(coefs)) keep <- keep & (is.na(coefs$rsq) | coefs$rsq >= min_model_rsq)
  if ("padj" %in% colnames(coefs)) {
    keep <- keep & !is.na(coefs$padj) & coefs$padj <= padj_threshold
  } else if (isTRUE(require_padj)) {
    stop("Pando network does not contain `padj`; use a p-value-producing model such as `method = 'glm'`, or set `require_padj = FALSE`.", call. = FALSE)
  }
  list(all = coefs, significant = coefs[keep, , drop = FALSE])
}

#' Project significant Pando edges onto a metabolic gene-gene network
#' @export
rc_project_metabolic_grn <- function(tf_peak_gene,
                                     metabolic_genes,
                                     top_k = 5L,
                                     min_shared_tfs = 1L,
                                     min_tf_jaccard = 0,
                                     max_targets_per_tf = 200L,
                                     include_direct_metabolic_tf = TRUE) {
  if (!is.data.frame(tf_peak_gene)) stop("`tf_peak_gene` must be a data.frame.", call. = FALSE)
  required <- c("sample_id", "tf", "target")
  missing <- setdiff(required, colnames(tf_peak_gene))
  if (length(missing)) stop("`tf_peak_gene` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  metabolic_genes <- unique(toupper(.rc_mm_trim_unique(metabolic_genes)))
  x <- tf_peak_gene
  x$tf <- toupper(as.character(x$tf))
  x$target <- toupper(as.character(x$target))
  x <- x[x$target %in% metabolic_genes, , drop = FALSE]
  if (!nrow(x)) {
    return(list(nodes = data.frame(sample_id = character(), gene = character(), node_role = character(), module_id = character(), stringsAsFactors = FALSE),
                edges = .rc_mm_empty_edges()))
  }
  strength <- if ("estimate" %in% colnames(x)) abs(as.numeric(x$estimate)) else rep(1, nrow(x))
  strength[!is.finite(strength)] <- 0
  x$.strength <- strength
  samples <- unique(as.character(x$sample_id))
  node_rows <- list()
  edge_rows <- list()

  for (sample in samples) {
    xs <- x[as.character(x$sample_id) == sample, , drop = FALSE]
    targets <- unique(xs$target)
    metabolic_tfs <- if (isTRUE(include_direct_metabolic_tf)) unique(xs$tf[xs$tf %in% metabolic_genes]) else character()
    nodes <- unique(c(targets, metabolic_tfs))
    tf_target <- stats::aggregate(xs$.strength,
                                  by = list(tf = xs$tf, target = xs$target),
                                  FUN = sum, na.rm = TRUE)
    colnames(tf_target)[[3L]] <- "strength"
    n_tfs <- table(tf_target$target)
    pair_acc <- list()
    pair_index <- 0L
    for (tf in unique(tf_target$tf)) {
      z <- tf_target[tf_target$tf == tf, , drop = FALSE]
      z <- z[order(z$strength, decreasing = TRUE), , drop = FALSE]
      if (nrow(z) > max_targets_per_tf) z <- z[seq_len(max_targets_per_tf), , drop = FALSE]
      if (nrow(z) < 2L) next
      cmb <- utils::combn(seq_len(nrow(z)), 2L)
      pair_index <- pair_index + 1L
      pair_acc[[pair_index]] <- data.frame(
        gene_a = pmin(z$target[cmb[1L, ]], z$target[cmb[2L, ]]),
        gene_b = pmax(z$target[cmb[1L, ]], z$target[cmb[2L, ]]),
        tf = tf,
        weight = pmin(z$strength[cmb[1L, ]], z$strength[cmb[2L, ]]),
        stringsAsFactors = FALSE
      )
    }
    shared <- if (length(pair_acc)) do.call(rbind, pair_acc) else data.frame()
    edges <- .rc_mm_empty_edges()[0, , drop = FALSE]
    if (nrow(shared)) {
      key <- paste(shared$gene_a, shared$gene_b, sep = "\001")
      split_rows <- split(seq_len(nrow(shared)), key)
      edges <- do.call(rbind, lapply(split_rows, function(ii) {
        a <- shared$gene_a[[ii[[1L]]]]
        b <- shared$gene_b[[ii[[1L]]]]
        shared_n <- length(unique(shared$tf[ii]))
        union_n <- as.integer(n_tfs[[a]]) + as.integer(n_tfs[[b]]) - shared_n
        data.frame(sample_id = sample, gene_a = a, gene_b = b,
                   edge_type = "shared_tf", shared_tf_count = shared_n,
                   projection_weight = sum(shared$weight[ii], na.rm = TRUE),
                   tf_jaccard = if (union_n > 0) shared_n / union_n else 0,
                   direct_regulatory = FALSE, module_id = NA_character_, stringsAsFactors = FALSE)
      }))
    }
    if (isTRUE(include_direct_metabolic_tf)) {
      direct <- unique(xs[xs$tf %in% metabolic_genes & xs$tf != xs$target,
                          c("tf", "target", ".strength"), drop = FALSE])
      if (nrow(direct)) {
        direct_edges <- stats::aggregate(direct$.strength,
                                         by = list(gene_a = pmin(direct$tf, direct$target),
                                                   gene_b = pmax(direct$tf, direct$target)),
                                         FUN = sum, na.rm = TRUE)
        colnames(direct_edges)[[3L]] <- "projection_weight"
        direct_edges$sample_id <- sample
        direct_edges$edge_type <- "direct_metabolic_tf"
        direct_edges$shared_tf_count <- 0L
        direct_edges$tf_jaccard <- 0
        direct_edges$direct_regulatory <- TRUE
        direct_edges$module_id <- NA_character_
        direct_edges <- direct_edges[, colnames(.rc_mm_empty_edges()), drop = FALSE]
        edges <- rbind(edges, direct_edges)
      }
    }
    if (nrow(edges)) {
      key <- paste(edges$gene_a, edges$gene_b, sep = "\001")
      edges <- do.call(rbind, lapply(split(seq_len(nrow(edges)), key), function(ii) {
        z <- edges[ii, , drop = FALSE]
        data.frame(sample_id = sample, gene_a = z$gene_a[[1L]], gene_b = z$gene_b[[1L]],
                   edge_type = paste(sort(unique(z$edge_type)), collapse = ";"),
                   shared_tf_count = max(z$shared_tf_count, na.rm = TRUE),
                   projection_weight = sum(z$projection_weight, na.rm = TRUE),
                   tf_jaccard = max(z$tf_jaccard, na.rm = TRUE),
                   direct_regulatory = any(z$direct_regulatory), module_id = NA_character_, stringsAsFactors = FALSE)
      }))
      edges <- edges[(edges$direct_regulatory |
                      (edges$shared_tf_count >= as.integer(min_shared_tfs) & edges$tf_jaccard >= min_tf_jaccard)), , drop = FALSE]
    }
    if (nrow(edges) && is.finite(top_k) && top_k > 0L) {
      selected <- rep(FALSE, nrow(edges))
      for (gene in nodes) {
        ii <- which(edges$gene_a == gene | edges$gene_b == gene)
        if (!length(ii)) next
        ord <- order(edges$direct_regulatory[ii], edges$projection_weight[ii], edges$shared_tf_count[ii], decreasing = TRUE)
        selected[ii[utils::head(ord, as.integer(top_k))]] <- TRUE
      }
      edges <- edges[selected, , drop = FALSE]
    }
    comps <- .rc_mm_components(nodes, edges)
    module_ids <- paste0(sample, "::GRN", sprintf("%04d", comps$component))
    role <- ifelse(comps$gene %in% targets & comps$gene %in% metabolic_tfs, "target_and_metabolic_tf",
                   ifelse(comps$gene %in% targets, "significant_target", "metabolic_tf_neighbor"))
    node_rows[[sample]] <- data.frame(sample_id = sample, gene = comps$gene,
                                      node_role = role, module_id = module_ids,
                                      stringsAsFactors = FALSE)
    if (nrow(edges)) {
      edges$module_id <- module_ids[match(edges$gene_a, comps$gene)]
      edge_rows[[sample]] <- edges
    }
  }
  list(nodes = do.call(rbind, node_rows),
       edges = if (length(edge_rows)) do.call(rbind, edge_rows) else .rc_mm_empty_edges())
}


#' Run sample-specific Pando GRNs and construct reaction meta-modules
#' @export
rc_run_pando_meta_modules <- function(metacell_object,
                                      gem,
                                      outdir,
                                      pfm,
                                      genome,
                                      sample_col = "sample_id",
                                      single_cell_genes = NULL,
                                      rna_assay = "RNA",
                                      atac_assay = "ATAC",
                                      pando_version = NULL,
                                      pando_remote_username = "1667857557",
                                      pando_remote_repo = "Pando_regcompass",
                                      require_pando_remote = TRUE,
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
                                      on_sample_error = c("record", "stop")) {
  expansion_mode <- match.arg(expansion_mode)
  on_sample_error <- match.arg(on_sample_error)
  if (!inherits(metacell_object, "Seurat")) stop("`metacell_object` must inherit from Seurat.", call. = FALSE)
  if (!requireNamespace("Pando", quietly = TRUE)) stop("Install the pinned Pando fork before running v1.2 meta-modules.", call. = FALSE)
  pando_install <- .rc_validate_pando_install(
    pando_version = pando_version,
    pando_remote_username = pando_remote_username,
    pando_remote_repo = pando_remote_repo,
    require_pando_remote = require_pando_remote
  )
  installed <- pando_install$version
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Signac", quietly = TRUE)) {
    stop("Packages 'Seurat' and 'Signac' are required.", call. = FALSE)
  }
  if (!sample_col %in% colnames(metacell_object@meta.data)) stop("Missing sample metadata column: ", sample_col, call. = FALSE)
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

  sample_ids <- .rc_mm_trim_unique(metacell_object@meta.data[[sample_col]])
  all_edges <- list()
  sig_edges <- list()
  status_rows <- list()

  for (sample in sample_ids) {
    cells <- rownames(metacell_object@meta.data)[as.character(metacell_object@meta.data[[sample_col]]) == sample]
    status <- data.frame(sample_id = sample, n_metacells = length(cells), n_target_genes = length(target_genes),
                         status = "pending", n_edges = 0L, n_significant_edges = 0L,
                         error_class = NA_character_, error_message = NA_character_, stringsAsFactors = FALSE)
    if (length(cells) < as.integer(min_metacells)) {
      status$status <- "skipped_too_few_metacells"
      status_rows[[sample]] <- status
      next
    }
    one <- tryCatch({
      obj <- subset(metacell_object, cells = cells)
      obj <- Seurat::NormalizeData(obj, assay = rna_assay, verbose = FALSE)
      obj <- Signac::RunTFIDF(obj, assay = atac_assay)
      sample_file_id <- gsub("[^A-Za-z0-9_.-]+", "_", sample)
      if (isTRUE(save_sample_metacell_objects)) {
        saveRDS(obj, file.path(outdir, "sample_metacell_objects", paste0(sample_file_id, ".rds")))
      }
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
                                           padj_threshold = padj_threshold,
                                           min_abs_estimate = min_abs_estimate,
                                           min_model_rsq = min_model_rsq,
                                           require_padj = require_padj)
      if (isTRUE(save_pando_objects)) saveRDS(grn, file.path(outdir, "pando_objects", paste0(sample_file_id, ".rds")))
      list(grn = grn, tables = tab)
    }, error = function(e) e)
    if (inherits(one, "error")) {
      status$status <- "failed"
      status$error_class <- class(one)[[1L]]
      status$error_message <- conditionMessage(one)
      status_rows[[sample]] <- status
      if (identical(on_sample_error, "stop")) stop(one)
      next
    }
    all_edges[[sample]] <- one$tables$all
    sig_edges[[sample]] <- one$tables$significant
    status$status <- "ok"
    status$n_edges <- nrow(one$tables$all)
    status$n_significant_edges <- nrow(one$tables$significant)
    status_rows[[sample]] <- status
  }

  status_table <- do.call(rbind, status_rows)
  all_table <- if (length(all_edges)) do.call(rbind, all_edges) else data.frame()
  sig_table <- if (length(sig_edges)) do.call(rbind, sig_edges) else data.frame()
  .rc_mm_write_tsv_gz(status_table, file.path(outdir, "pando_sample_status.tsv.gz"))
  .rc_mm_write_tsv_gz(all_table, file.path(outdir, "pando_tf_peak_gene_all.tsv.gz"))
  .rc_mm_write_tsv_gz(sig_table, file.path(outdir, "pando_tf_peak_gene_significant.tsv.gz"))
  if (!nrow(sig_table)) stop("No significant Pando TF-peak-gene edges were available across samples.", call. = FALSE)

  projection <- rc_project_metabolic_grn(sig_table, metabolic_genes = metabolic_genes,
                                          top_k = top_k_neighbors,
                                          min_shared_tfs = min_shared_tfs,
                                          min_tf_jaccard = min_tf_jaccard,
                                          max_targets_per_tf = max_targets_per_tf,
                                          include_direct_metabolic_tf = TRUE)
  core <- rc_map_meta_module_core_reactions(projection$nodes, gem$gpr_table)
  expanded <- rc_expand_meta_module_reactions(gem, core,
                                               subsystem_table = subsystem_table,
                                               expansion_mode = expansion_mode)
  .rc_mm_write_tsv_gz(projection$nodes, file.path(outdir, "metabolic_gene_nodes.tsv.gz"))
  .rc_mm_write_tsv_gz(projection$edges, file.path(outdir, "metabolic_gene_edges.tsv.gz"))
  .rc_mm_write_tsv_gz(core, file.path(outdir, "core_gene_reaction.tsv.gz"))
  .rc_mm_write_tsv_gz(expanded$reaction_membership, file.path(outdir, "meta_module_reactions.tsv.gz"))
  .rc_mm_write_tsv_gz(expanded$summary, file.path(outdir, "meta_module_summary.tsv.gz"))

  out <- list(schema_version = "regcompass_pando_meta_module_v1.2",
              pando_version = installed,
              pando_remote_username = pando_install$remote_username,
              pando_remote_repo = pando_install$remote_repo,
              pando_remote_ref = pando_install$remote_ref,
              pando_remote_sha = pando_install$remote_sha,
              target_metabolic_genes = target_genes,
              sample_status = status_table,
              tf_peak_gene_all = all_table,
              tf_peak_gene_significant = sig_table,
              metabolic_gene_nodes = projection$nodes,
              metabolic_gene_edges = projection$edges,
              core_gene_reaction = core,
              reaction_membership = expanded$reaction_membership,
              meta_module_summary = expanded$summary,
              crossref_maps = expanded$crossref_maps)
  saveRDS(out, file.path(outdir, "pando_meta_modules.rds"))
  out
}

