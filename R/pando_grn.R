#' Extract and filter a Pando TF-peak-gene coefficient table
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
.rc_project_metabolic_grn_base <- function(tf_peak_gene,
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


.rc_signed_relation <- function(values) {
  values <- values[is.finite(values) & values != 0]
  if (!length(values)) return(NA_character_)
  if (all(values > 0)) return("concordant")
  if (all(values < 0)) return("discordant")
  "mixed"
}

.rc_project_metabolic_grn_signed_metadata <- function(tf_peak_gene, metabolic_genes,
                                     top_k = 5L, min_shared_tfs = 1L,
                                     min_tf_jaccard = 0,
                                     max_targets_per_tf = 200L,
                                     include_direct_metabolic_tf = TRUE) {
  answer <- .rc_project_metabolic_grn_base(
    tf_peak_gene, metabolic_genes,
    top_k = Inf,
    min_shared_tfs = min_shared_tfs,
    min_tf_jaccard = min_tf_jaccard,
    max_targets_per_tf = max_targets_per_tf,
    include_direct_metabolic_tf = include_direct_metabolic_tf
  )
  edges <- answer$edges
  edges$regulator_set <- NA_character_
  edges$direct_regulator <- NA_character_
  edges$direct_target <- NA_character_
  edges$regulatory_relation <- NA_character_
  edges$signed_projection_weight <- NA_real_
  edges$direction_and_sign_preserved <- FALSE
  if (!nrow(edges) || !is.data.frame(tf_peak_gene) ||
      !all(c("sample_id", "tf", "target") %in% colnames(tf_peak_gene))) {
    answer$edges <- edges
    return(answer)
  }
  x <- tf_peak_gene
  x$sample_id <- as.character(x$sample_id)
  x$tf <- toupper(trimws(as.character(x$tf)))
  x$target <- toupper(trimws(as.character(x$target)))
  x$.estimate <- if ("estimate" %in% colnames(x))
    suppressWarnings(as.numeric(x$estimate)) else rep(NA_real_, nrow(x))
  x$.strength <- abs(x$.estimate)

  for (i in seq_len(nrow(edges))) {
    sample <- as.character(edges$sample_id[[i]])
    a <- toupper(as.character(edges$gene_a[[i]]))
    b <- toupper(as.character(edges$gene_b[[i]]))
    xs <- x[x$sample_id == sample, , drop = FALSE]
    direct <- xs[(xs$tf == a & xs$target == b) |
                   (xs$tf == b & xs$target == a), , drop = FALSE]
    if (nrow(direct)) {
      edges$direct_regulator[[i]] <- paste(unique(direct$tf), collapse = ";")
      edges$direct_target[[i]] <- paste(unique(direct$target), collapse = ";")
      edges$regulator_set[[i]] <- paste(unique(direct$tf), collapse = ";")
      edges$regulatory_relation[[i]] <- .rc_signed_relation(direct$.estimate)
      edges$signed_projection_weight[[i]] <- sum(direct$.estimate, na.rm = TRUE)
      edges$direction_and_sign_preserved[[i]] <- TRUE
      next
    }
    xa <- xs[xs$target == a, c("tf", ".estimate", ".strength"), drop = FALSE]
    xb <- xs[xs$target == b, c("tf", ".estimate", ".strength"), drop = FALSE]
    shared <- intersect(unique(xa$tf), unique(xb$tf))
    if (!length(shared)) next
    contributions <- vapply(shared, function(tf) {
      ea <- sum(xa$.estimate[xa$tf == tf], na.rm = TRUE)
      eb <- sum(xb$.estimate[xb$tf == tf], na.rm = TRUE)
      sa <- sum(xa$.strength[xa$tf == tf], na.rm = TRUE)
      sb <- sum(xb$.strength[xb$tf == tf], na.rm = TRUE)
      sign(ea) * sign(eb) * min(sa, sb)
    }, numeric(1))
    edges$regulator_set[[i]] <- paste(shared, collapse = ";")
    edges$regulatory_relation[[i]] <- .rc_signed_relation(contributions)
    edges$signed_projection_weight[[i]] <- sum(contributions, na.rm = TRUE)
    edges$direction_and_sign_preserved[[i]] <- TRUE
  }
  answer$edges <- edges
  answer
}

rc_project_metabolic_grn <- function(tf_peak_gene, metabolic_genes,
                                     top_k = 5L, min_shared_tfs = 1L,
                                     min_tf_jaccard = 0,
                                     max_targets_per_tf = 200L,
                                     include_direct_metabolic_tf = TRUE) {
  answer <- .rc_project_metabolic_grn_signed_metadata(
    tf_peak_gene = tf_peak_gene,
    metabolic_genes = metabolic_genes,
    top_k = top_k,
    min_shared_tfs = min_shared_tfs,
    min_tf_jaccard = min_tf_jaccard,
    max_targets_per_tf = max_targets_per_tf,
    include_direct_metabolic_tf = include_direct_metabolic_tf
  )
  edges <- answer$edges
  nodes <- answer$nodes
  if (!nrow(edges) || !nrow(nodes)) {
    edges$used_for_component <- logical(nrow(edges))
    answer$edges <- edges
    return(answer)
  }

  relation <- as.character(edges$regulatory_relation)
  component_candidate <- edges$direct_regulatory %in% TRUE |
    (!edges$direct_regulatory %in% TRUE & relation == "concordant")
  component_candidate[is.na(component_candidate)] <- FALSE
  edges$used_for_component <- component_candidate
  if (nrow(edges) && is.finite(top_k) && top_k > 0L) {
    selected <- rep(FALSE, nrow(edges))
    for (sample in unique(as.character(nodes$sample_id))) {
      sample_edge <- as.character(edges$sample_id) == sample
      sample_genes <- as.character(
        nodes$gene[as.character(nodes$sample_id) == sample]
      )
      for (gene in sample_genes) {
        index <- which(
          sample_edge & component_candidate &
            (edges$gene_a == gene | edges$gene_b == gene)
        )
        if (!length(index)) next
        order_index <- order(
          edges$direct_regulatory[index],
          edges$projection_weight[index],
          edges$shared_tf_count[index],
          decreasing = TRUE,
          na.last = TRUE
        )
        selected[index[utils::head(order_index, as.integer(top_k))]] <- TRUE
      }
    }
    edges$used_for_component <- selected
  }

  for (sample in unique(as.character(nodes$sample_id))) {
    node_index <- as.character(nodes$sample_id) == sample
    edge_index <- as.character(edges$sample_id) == sample
    component_edges <- edges[
      edge_index & edges$used_for_component,
      , drop = FALSE
    ]
    components <- .rc_mm_components(
      as.character(nodes$gene[node_index]),
      component_edges
    )
    module_ids <- paste0(
      sample, "::GRN", sprintf("%04d", components$component)
    )
    nodes$module_id[node_index] <- module_ids[
      match(as.character(nodes$gene[node_index]), components$gene)
    ]
    if (any(edge_index)) {
      sample_nodes <- nodes[node_index, , drop = FALSE]
      edges$module_id[edge_index] <- sample_nodes$module_id[
        match(
          as.character(edges$gene_a[edge_index]),
          as.character(sample_nodes$gene)
        )
      ]
    }
  }

  answer$nodes <- nodes
  answer$edges <- edges
  answer$component_policy <- paste(
    "direct regulatory edges plus concordant signed shared-TF edges;",
    "discordant and mixed shared-TF edges are diagnostic only"
  )
  answer
}

rc_run_pando_meta_modules <- function(metacell_object,
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

  meta <- metacell_object@meta.data
  meta$.rc_pando_group_id <- rc_make_stratum_id(meta, group_cols)
  group_ids <- unique(as.character(meta$.rc_pando_group_id))

  run_one_group <- function(group_id) {
    cells <- rownames(meta)[as.character(meta$.rc_pando_group_id) == group_id]
    vals <- meta[match(cells[[1L]], rownames(meta)), group_cols, drop = FALSE]
    sample <- as.character(vals[[sample_col]][[1L]])
    status <- data.frame(
      group_id = group_id,
      vals,
      n_metacells = length(cells),
      n_target_genes = length(target_genes),
      status = "pending",
      n_edges = 0L,
      n_significant_edges = 0L,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (length(cells) < as.integer(min_metacells)) {
      status$status <- "skipped_too_few_metacells"
      return(list(status = status, all = data.frame(), significant = data.frame()))
    }
    one <- tryCatch({
      obj <- subset(metacell_object, cells = cells)
      obj <- Seurat::NormalizeData(obj, assay = rna_assay, verbose = FALSE)
      obj <- Signac::RunTFIDF(obj, assay = atac_assay)
      group_file_id <- gsub("[^A-Za-z0-9_.-]+", "_", group_id)
      if (isTRUE(save_sample_metacell_objects)) {
        saveRDS(obj, file.path(outdir, "sample_metacell_objects", paste0(group_file_id, ".rds")))
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
      add_group_meta <- function(x) {
        if (!nrow(x)) return(x)
        x$group_id <- group_id
        for (col in group_cols) x[[col]] <- vals[[col]][[1L]]
        x[, c(display_cols, setdiff(colnames(x), display_cols)), drop = FALSE]
      }
      tab$all <- add_group_meta(tab$all)
      tab$significant <- add_group_meta(tab$significant)
      if (isTRUE(save_pando_objects)) saveRDS(grn, file.path(outdir, "pando_objects", paste0(group_file_id, ".rds")))
      tab
    }, error = function(e) e)
    if (inherits(one, "error")) {
      status$status <- "failed"
      status$error_class <- class(one)[[1L]]
      status$error_message <- conditionMessage(one)
      if (identical(on_sample_error, "stop")) stop(one)
      return(list(status = status, all = data.frame(), significant = data.frame()))
    }
    status$status <- "ok"
    status$n_edges <- nrow(one$all)
    status$n_significant_edges <- nrow(one$significant)
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

  sig_for_projection <- sig_table
  sig_for_projection$sample_id <- sig_for_projection$group_id
  projection <- rc_project_metabolic_grn(sig_for_projection, metabolic_genes = metabolic_genes,
                                          top_k = top_k_neighbors,
                                          min_shared_tfs = min_shared_tfs,
                                          min_tf_jaccard = min_tf_jaccard,
                                          max_targets_per_tf = max_targets_per_tf,
                                          include_direct_metabolic_tf = TRUE)
  group_meta <- unique(status_table[, c("group_id", group_cols), drop = FALSE])
  remap_projection <- function(x) {
    if (!nrow(x)) return(x)
    x$group_id <- as.character(x$sample_id)
    x <- merge(x, group_meta, by = "group_id", all.x = TRUE, sort = FALSE)
    x$sample_id <- as.character(x[[sample_col]])
    x[, c(display_cols, setdiff(colnames(x), display_cols)), drop = FALSE]
  }
  projection$nodes <- remap_projection(projection$nodes)
  projection$edges <- remap_projection(projection$edges)
  core <- rc_map_meta_module_core_reactions(projection$nodes, gem$gpr_table)
  if (nrow(core)) {
    core <- merge(core, unique(projection$nodes[, module_cols, drop = FALSE]),
                  by = c("sample_id", "module_id"), all.x = TRUE, sort = FALSE)
    core <- core[, c(display_cols, setdiff(colnames(core), display_cols)), drop = FALSE]
  }
  expanded <- rc_expand_meta_module_reactions(gem, core,
                                               subsystem_table = subsystem_table,
                                               expansion_mode = expansion_mode)
  if (nrow(expanded$reaction_membership)) {
    expanded$reaction_membership <- merge(expanded$reaction_membership, unique(core[, module_cols, drop = FALSE]),
                                          by = c("sample_id", "module_id"), all.x = TRUE, sort = FALSE)
    expanded$reaction_membership <- expanded$reaction_membership[, c(display_cols, setdiff(colnames(expanded$reaction_membership), display_cols)), drop = FALSE]
  }
  if (nrow(expanded$summary)) {
    expanded$summary <- merge(expanded$summary, unique(core[, module_cols, drop = FALSE]),
                              by = c("sample_id", "module_id"), all.x = TRUE, sort = FALSE)
    expanded$summary <- expanded$summary[, c(display_cols, setdiff(colnames(expanded$summary), display_cols)), drop = FALSE]
  }
  .rc_mm_write_tsv_gz(projection$nodes, file.path(outdir, "metabolic_gene_nodes.tsv.gz"))
  .rc_mm_write_tsv_gz(projection$edges, file.path(outdir, "metabolic_gene_edges.tsv.gz"))
  .rc_mm_write_tsv_gz(core, file.path(outdir, "core_gene_reaction.tsv.gz"))
  .rc_mm_write_tsv_gz(expanded$reaction_membership, file.path(outdir, "meta_module_reactions.tsv.gz"))
  .rc_mm_write_tsv_gz(expanded$summary, file.path(outdir, "meta_module_summary.tsv.gz"))

  out <- list(schema_version = "regcompass_pando_meta_module_v1.2",
              pando_installed_version = installed,
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
