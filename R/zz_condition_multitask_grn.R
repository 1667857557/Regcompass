# Unified condition-aware Pando integration and TF-by-ATAC penalty projection.

.rc_condition_network_slot <- function(grn_object, network_id) {
  grn <- methods::slot(grn_object, "grn")
  networks <- methods::slot(grn, "networks")
  if (!network_id %in% names(networks)) {
    stop("Pando condition network was not found: ", network_id, call. = FALSE)
  }
  networks[[network_id]]
}

.rc_condition_network_tables <- function(grn_object, network_id) {
  network <- .rc_condition_network_slot(grn_object, network_id)
  coefs <- as.data.frame(methods::slot(network, "coefs"),
                         stringsAsFactors = FALSE)
  fit <- as.data.frame(methods::slot(network, "fit"),
                       stringsAsFactors = FALSE)
  required <- c("tf", "target", "region", "estimate")
  missing <- setdiff(required, colnames(coefs))
  if (length(missing)) {
    stop(
      "Pando condition network `", network_id,
      "` lacks coefficient columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"term" %in% colnames(coefs)) {
    coefs$term <- paste(coefs$region, coefs$tf, sep = ":")
  }
  if (!"corr" %in% colnames(coefs)) coefs$corr <- NA_real_
  if (!all(c("target", "rsq") %in% colnames(fit))) {
    stop(
      "Pando condition network `", network_id,
      "` lacks target-level `target` and `rsq` fit fields.",
      call. = FALSE
    )
  }
  coefs$tf <- toupper(trimws(as.character(coefs$tf)))
  coefs$target <- toupper(trimws(as.character(coefs$target)))
  coefs$region <- trimws(as.character(coefs$region))
  coefs$term <- as.character(coefs$term)
  coefs$estimate <- suppressWarnings(as.numeric(coefs$estimate))
  coefs$corr <- suppressWarnings(as.numeric(coefs$corr))
  fit$target <- toupper(trimws(as.character(fit$target)))
  fit$rsq <- suppressWarnings(as.numeric(fit$rsq))
  list(coefs = coefs, fit = fit)
}

.rc_condition_edge_key <- function(x) {
  paste(x$tf, x$target, x$region, x$term, sep = "\001")
}

.rc_build_condition_effect_table <- function(
    condition_table, universal_table, condition_fit, universal_fit) {
  condition_table <- as.data.frame(condition_table, stringsAsFactors = FALSE)
  universal_table <- as.data.frame(universal_table, stringsAsFactors = FALSE)
  key_condition <- .rc_condition_edge_key(condition_table)
  key_universal <- .rc_condition_edge_key(universal_table)
  if (anyDuplicated(key_condition) || anyDuplicated(key_universal)) {
    stop("Pando condition/universal networks contain duplicated TF-peak-target edges.",
         call. = FALSE)
  }
  if (!setequal(key_condition, key_universal)) {
    stop(
      "Pando condition and universal networks do not share one edge dictionary.",
      call. = FALSE
    )
  }
  universal_table <- universal_table[
    match(key_condition, key_universal), , drop = FALSE
  ]
  condition_estimate <- suppressWarnings(as.numeric(condition_table$estimate))
  universal_estimate <- suppressWarnings(as.numeric(universal_table$estimate))
  out <- condition_table
  out$condition_estimate <- condition_estimate
  out$universal_estimate <- universal_estimate
  out$condition_effect <- condition_estimate - universal_estimate
  out$condition_corr <- suppressWarnings(as.numeric(condition_table$corr))
  out$universal_corr <- suppressWarnings(as.numeric(universal_table$corr))
  out$estimate <- out$condition_estimate
  condition_fit <- condition_fit[!duplicated(condition_fit$target), , drop = FALSE]
  universal_fit <- universal_fit[!duplicated(universal_fit$target), , drop = FALSE]
  out$rsq <- condition_fit$rsq[match(out$target, condition_fit$target)]
  out$universal_rsq <- universal_fit$rsq[
    match(out$target, universal_fit$target)
  ]
  out
}

.rc_extract_condition_multitask_grn <- function(
    grn_object, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  grn <- methods::slot(grn_object, "grn")
  params <- methods::slot(grn, "params")
  index <- params$condition_network_index
  diagnostics <- params$condition_fit_diagnostics
  if (!is.data.frame(index) ||
      !all(c("network_id", "cell_type", "network_level", "condition") %in%
           colnames(index))) {
    stop(
      "Installed Pando did not return `condition_network_index`; install ",
      "1667857557/Pando_regcompass at or after d6cce000.",
      call. = FALSE
    )
  }
  condition_rows <- index[index$network_level == "condition", , drop = FALSE]
  universal_rows <- index[index$network_level == "universal", , drop = FALSE]
  if (!nrow(condition_rows) || !nrow(universal_rows)) {
    stop("Pando returned no condition or universal networks.", call. = FALSE)
  }
  active_tol <- max(
    suppressWarnings(as.numeric(min_abs_estimate)),
    1e-8,
    na.rm = TRUE
  )
  all_rows <- vector("list", nrow(condition_rows))
  universal_output <- vector("list", nrow(universal_rows))
  for (i in seq_len(nrow(universal_rows))) {
    one <- universal_rows[i, , drop = FALSE]
    table <- .rc_condition_network_tables(
      grn_object, as.character(one$network_id[[1L]])
    )
    tab <- table$coefs
    tab$rsq <- table$fit$rsq[match(tab$target, table$fit$target)]
    tab$cell_type <- as.character(one$cell_type[[1L]])
    tab$network_id <- as.character(one$network_id[[1L]])
    tab$network_level <- "universal"
    universal_output[[i]] <- tab
  }
  for (i in seq_len(nrow(condition_rows))) {
    one <- condition_rows[i, , drop = FALSE]
    cell_type <- as.character(one$cell_type[[1L]])
    condition <- as.character(one$condition[[1L]])
    universal_row <- universal_rows[
      as.character(universal_rows$cell_type) == cell_type, , drop = FALSE
    ]
    if (nrow(universal_row) != 1L) {
      stop(
        "Pando must return exactly one universal network per cell type; cell type ",
        cell_type, " had ", nrow(universal_row), ".", call. = FALSE
      )
    }
    condition_id <- as.character(one$network_id[[1L]])
    universal_id <- as.character(universal_row$network_id[[1L]])
    condition_data <- .rc_condition_network_tables(grn_object, condition_id)
    universal_data <- .rc_condition_network_tables(grn_object, universal_id)
    tab <- .rc_build_condition_effect_table(
      condition_data$coefs, universal_data$coefs,
      condition_data$fit, universal_data$fit
    )
    group_values <- data.frame(
      condition_value = condition,
      celltype_value = cell_type,
      stringsAsFactors = FALSE
    )
    names(group_values) <- c(condition_col, celltype_col)
    tab$group_id <- rc_make_stratum_id(group_values, c(condition_col, celltype_col))
    tab[[condition_col]] <- condition
    tab[[celltype_col]] <- cell_type
    tab$condition_network_id <- condition_id
    tab$universal_network_id <- universal_id
    tab$network_level <- "condition"
    tab$fit_engine <- "multitask_sparse_group_glmnet"
    tab <- tab[, c(
      "group_id", condition_col, celltype_col,
      setdiff(colnames(tab), c("group_id", condition_col, celltype_col))
    ), drop = FALSE]
    all_rows[[i]] <- tab
  }
  all_edges <- do.call(rbind, all_rows)
  rownames(all_edges) <- NULL
  condition_active <- all_edges[
    is.finite(all_edges$condition_estimate) &
      abs(all_edges$condition_estimate) >= active_tol &
      is.finite(all_edges$rsq) & all_edges$rsq >= min_model_rsq,
    , drop = FALSE
  ]
  effect_all <- all_edges
  effect_all$estimate <- effect_all$condition_effect
  effect_active <- effect_all[
    is.finite(effect_all$condition_effect) &
      abs(effect_all$condition_effect) >= active_tol &
      is.finite(effect_all$rsq) & effect_all$rsq >= min_model_rsq,
    , drop = FALSE
  ]
  universal <- do.call(rbind, universal_output)
  rownames(universal) <- NULL
  list(
    network_index = index,
    fit_diagnostics = diagnostics %||% data.frame(),
    universal = universal,
    condition_all = all_edges,
    condition_active = condition_active,
    condition_effect_all = effect_all,
    condition_effect_active = effect_active,
    active_tol = active_tol
  )
}

.rc_run_condition_single_cell_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list(
      method = "multitask_glmnet",
      candidate_screen = "condition_union",
      tf_cor = 0.1,
      peak_cor = 0.01,
      alpha = 0.5,
      condition_mix = 0.5,
      condition_weight = "equal",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = FALSE,
      parallel = FALSE
    ),
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = FALSE,
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    on_group_error = c("record", "stop"),
    species = c("auto", "human", "mouse")) {
  on_group_error <- match.arg(on_group_error)
  species <- .rc_infer_gem_species(gem, species)
  if (!is.numeric(min_cells) || length(min_cells) != 1L ||
      !is.finite(min_cells) || min_cells < 1 ||
      abs(min_cells - round(min_cells)) > sqrt(.Machine$double.eps)) {
    stop("`min_cells` must be one positive integer.", call. = FALSE)
  }
  min_cells <- as.integer(min_cells)
  if (!is.list(pando_initiate_args) || !is.list(pando_motif_args) ||
      !is.list(pando_infer_args)) {
    stop("Pando initiate, motif, and inference arguments must be lists.",
         call. = FALSE)
  }
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop(
      "Install 1667857557/Pando_regcompass before condition-aware GRN inference.",
      call. = FALSE
    )
  }
  if (!exists("infer_condition_grn", envir = asNamespace("Pando"),
              inherits = FALSE)) {
    stop(
      "Installed Pando lacks `infer_condition_grn`; install commit d6cce000 or later.",
      call. = FALSE
    )
  }
  pando_install <- .rc_validate_pando_repository()
  motif_policy <- "user_supplied"
  if (is.null(pfm)) {
    pfm <- .rc_default_pando_motifs()
    motif_policy <- "Pando::motifs"
  }
  if (!length(pfm)) stop("`pfm` must be non-empty.", call. = FALSE)
  group_cols <- c(condition_col, celltype_col)
  missing <- setdiff(group_cols, colnames(object@meta.data))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  .rc_require_normalized_assay(object, rna_assay, "RNA")
  .rc_require_normalized_assay(object, atac_assay, "ATAC")
  normalization <- object@misc$regcompass_atac_normalization %||% list()
  if (!identical(normalization$scope, "cell_type_across_conditions")) {
    stop(
      "Pando condition inference requires cell-type-shared ATAC TF-IDF across conditions.",
      call. = FALSE
    )
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "pando_objects"),
             recursive = TRUE, showWarnings = FALSE)

  region_policy <- "user_supplied"
  if (!"regions" %in% names(pando_initiate_args) ||
      is.null(pando_initiate_args$regions)) {
    pando_initiate_args$regions <- .rc_default_pando_regions(species)
    region_policy <- if (identical(species, "human")) {
      paste(
        "union(Pando::phastConsElements20Mammals.UCSC.hg38,",
        "Pando::SCREEN.ccRE.UCSC.hg38)"
      )
    } else {
      "Pando::phastConsElements20Mammals.UCSC.hg38"
    }
  }

  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  target_upper <- intersect(
    toupper(.rc_mm_trim_unique(rna_genes)),
    toupper(.rc_mm_trim_unique(metabolic_genes))
  )
  target_genes <- .rc_mm_trim_unique(
    rna_genes[toupper(rna_genes) %in% target_upper]
  )
  if (!length(target_genes)) {
    stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)
  }

  filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Pando unified condition GRN"
  )
  object <- filtered$object
  init_defaults <- list(
    object = object,
    peak_assay = atac_assay,
    rna_assay = rna_assay
  )
  init_defaults[names(pando_initiate_args)] <- NULL
  grn <- do.call(Pando::initiate_grn, c(init_defaults, pando_initiate_args))
  motif_defaults <- list(object = grn, pfm = pfm, genome = genome)
  motif_defaults[names(pando_motif_args)] <- NULL
  grn <- do.call(Pando::find_motifs, c(motif_defaults, pando_motif_args))

  infer_defaults <- list(
    object = grn,
    cell_type_col = celltype_col,
    condition_col = condition_col,
    genes = target_genes,
    network_name = "regcompass_condition_grn",
    min_cells_per_condition = min_cells,
    on_small_condition = if (identical(on_group_error, "stop")) "error" else
      "skip_cell_type",
    BPPARAM = if (identical(BPPARAM, FALSE)) NULL else BPPARAM
  )
  infer_defaults[names(pando_infer_args)] <- NULL
  grn <- do.call(
    Pando::infer_condition_grn,
    c(infer_defaults, pando_infer_args)
  )

  extracted <- .rc_extract_condition_multitask_grn(
    grn,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq
  )
  meta <- object@meta.data
  expected <- unique(meta[, group_cols, drop = FALSE])
  expected$n_cells <- vapply(seq_len(nrow(expected)), function(i) {
    sum(
      as.character(meta[[condition_col]]) ==
        as.character(expected[[condition_col]][[i]]) &
      as.character(meta[[celltype_col]]) ==
        as.character(expected[[celltype_col]][[i]])
    )
  }, integer(1))
  expected$group_id <- rc_make_stratum_id(expected, group_cols)
  condition_index <- extracted$network_index[
    extracted$network_index$network_level == "condition", , drop = FALSE
  ]
  index_key <- paste(
    as.character(condition_index$condition),
    as.character(condition_index$cell_type),
    sep = "\001"
  )
  expected_key <- paste(
    as.character(expected[[condition_col]]),
    as.character(expected[[celltype_col]]),
    sep = "\001"
  )
  index_match <- match(expected_key, index_key)
  status <- expected
  status$n_target_genes <- length(target_genes)
  status$n_atac_peaks_input <- filtered$diagnostics$n_input_peaks
  status$n_zero_count_peaks_excluded <-
    filtered$diagnostics$n_zero_count_peaks_excluded
  status$n_atac_peaks_used <- filtered$diagnostics$n_retained_peaks
  status$status <- ifelse(is.na(index_match), "failed_missing_condition_network", "ok")
  status$n_edges <- vapply(status$group_id, function(id) {
    sum(extracted$condition_all$group_id == id)
  }, integer(1))
  status$n_significant_edges <- vapply(status$group_id, function(id) {
    sum(extracted$condition_active$group_id == id)
  }, integer(1))
  status$n_condition_effect_edges <- vapply(status$group_id, function(id) {
    sum(extracted$condition_effect_active$group_id == id)
  }, integer(1))
  status$error_class <- ifelse(status$status == "ok", NA_character_,
                               "missing_condition_network")
  status$error_message <- ifelse(
    status$status == "ok", NA_character_,
    "Pando did not emit a condition network for this condition-by-cell-type group."
  )
  status <- status[, c(
    "group_id", group_cols,
    setdiff(colnames(status), c("group_id", group_cols))
  ), drop = FALSE]
  if (any(status$status != "ok")) {
    stop(
      "Unified Pando condition GRN coverage is incomplete: ",
      paste(status$group_id[status$status != "ok"], collapse = "; "),
      call. = FALSE
    )
  }

  .rc_mm_write_tsv_gz(
    status, file.path(outdir, "pando_group_status.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$condition_all,
    file.path(outdir, "pando_tf_peak_gene_condition_all.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$condition_active,
    file.path(outdir, "pando_tf_peak_gene_condition_active.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$condition_effect_all,
    file.path(outdir, "pando_tf_peak_gene_condition_effect_all.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$condition_effect_active,
    file.path(outdir, "pando_tf_peak_gene_condition_effect_active.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$universal,
    file.path(outdir, "pando_tf_peak_gene_universal.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$network_index,
    file.path(outdir, "pando_condition_network_index.tsv.gz")
  )
  if (is.data.frame(extracted$fit_diagnostics) &&
      nrow(extracted$fit_diagnostics)) {
    .rc_mm_write_tsv_gz(
      extracted$fit_diagnostics,
      file.path(outdir, "pando_condition_fit_diagnostics.tsv.gz")
    )
  }
  if (isTRUE(save_pando_objects)) {
    saveRDS(grn, file.path(outdir, "pando_objects",
                           "condition_multitask_grn.rds"))
  }

  tf_metabolic_target_overlap <- intersect(
    unique(toupper(extracted$condition_effect_active$tf)),
    unique(toupper(target_genes))
  )
  answer <- list(
    schema_version = "regcompass_condition_multitask_grn_v1",
    pando_installed_version = pando_install$version,
    pando_installation = pando_install,
    target_metabolic_genes = target_genes,
    sample_status = status,
    pando_network_index = extracted$network_index,
    pando_fit_diagnostics = extracted$fit_diagnostics,
    tf_peak_gene_universal = extracted$universal,
    tf_peak_gene_condition_all = extracted$condition_all,
    tf_peak_gene_condition = extracted$condition_active,
    tf_peak_gene_condition_effect_all = extracted$condition_effect_all,
    tf_peak_gene_condition_effect = extracted$condition_effect_active,
    tf_metabolic_target_overlap = tf_metabolic_target_overlap,
    tf_peak_gene_all = extracted$condition_all,
    tf_peak_gene_significant = extracted$condition_active,
    normalization_policy = list(
      rna = "global single-cell normalized RNA shared by all conditions",
      atac = "cell-type-shared TF-IDF across conditions",
      grn_fit = paste(
        "one sparse-group multi-task Pando fit per cell type with a shared",
        "TF-peak-target edge dictionary and jointly estimated condition coefficients"
      ),
      universal_coefficient = "row mean of jointly estimated condition coefficients",
      condition_effect = "condition coefficient minus universal coefficient",
      core_reaction_evidence =
        "active condition-level TF-peak-target coefficients",
      penalty_regulatory_evidence = paste(
        "active condition-effect coefficients projected through metacell",
        "TF RNA by peak ATAC interaction activity"
      ),
      pando_motifs = motif_policy,
      pando_regions = region_policy,
      pando_padj = paste(
        "not applicable to the regularized multi-task solver;",
        "active coefficients and target-level R-squared are used"
      ),
      legacy_padj_threshold = padj_threshold,
      legacy_require_padj = require_padj,
      active_tolerance = extracted$active_tol,
      min_model_rsq = min_model_rsq
    ),
    group_cols = group_cols
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_tf_peak_interaction <- function(tf_activity, peak_activity) {
  tf_activity <- as.matrix(tf_activity)
  peak_activity <- as.matrix(peak_activity)
  if (!identical(dim(tf_activity), dim(peak_activity))) {
    stop("TF RNA and peak ATAC edge matrices must have identical dimensions.",
         call. = FALSE)
  }
  if (!identical(colnames(tf_activity), colnames(peak_activity))) {
    stop("TF RNA and peak ATAC matrices must contain identically ordered units.",
         call. = FALSE)
  }
  tf_activity * peak_activity
}

.rc_condition_gene_regulatory_modifier <- function(
    significant_edges, object, unit_meta,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    target_genes = NULL, min_scale = 0.05) {
  if (!is.data.frame(significant_edges)) {
    stop("`significant_edges` must be a data.frame.", call. = FALSE)
  }
  required_edges <- c(
    "target", "region", "tf", "estimate", condition_col, celltype_col
  )
  missing_edges <- setdiff(required_edges, colnames(significant_edges))
  if (length(missing_edges)) {
    stop(
      "Condition-effect edge table is missing columns: ",
      paste(missing_edges, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.data.frame(unit_meta) ||
      !all(c("pool_id", condition_col, celltype_col) %in% colnames(unit_meta))) {
    stop("`unit_meta` is incomplete for condition-pooled regulatory scoring.",
         call. = FALSE)
  }
  units <- colnames(object)
  unit_meta <- unit_meta[
    match(units, as.character(unit_meta$pool_id)), , drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Metacell metadata do not align to the scoring object.", call. = FALSE)
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
  attr(modifier, "reliability_policy") <- paste(
    "finite condition-network R-squared is mapped through sqrt(clamp(median R2,0,1));",
    "targets without finite R-squared receive regulatory reliability zero"
  )
  attr(modifier, "regulatory_predictor") <- paste(
    "condition_effect(TF-peak-target coefficient) multiplied by robust",
    "cell-type-referenced deviation of TF_RNA x peak_ATAC activity"
  )
  if (!nrow(significant_edges) || !length(genes)) return(modifier)

  edges <- significant_edges
  rsq <- if ("rsq" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$rsq))
  } else {
    rep(NA_real_, nrow(edges))
  }
  edges <- edges[is.finite(rsq), , drop = FALSE]
  if (!nrow(edges)) return(modifier)
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

  rna <- .rc_pando_assay_data(object, rna_assay)
  atac <- .rc_pando_assay_data(object, atac_assay)
  if (!setequal(colnames(rna), units) || !setequal(colnames(atac), units)) {
    stop("Normalized RNA/ATAC data and metacell units do not match.",
         call. = FALSE)
  }
  rna <- rna[, units, drop = FALSE]
  atac <- atac[, units, drop = FALSE]
  tf_keys <- toupper(trimws(rownames(rna)))
  tf_keep <- !is.na(tf_keys) & nzchar(tf_keys) & !duplicated(tf_keys)
  tf_lookup <- stats::setNames(rownames(rna)[tf_keep], tf_keys[tf_keep])
  peak_keys <- toupper(.rc_pando_region_key(rownames(atac)))
  peak_keep <- !is.na(peak_keys) & nzchar(peak_keys) & !duplicated(peak_keys)
  peak_lookup <- stats::setNames(rownames(atac)[peak_keep], peak_keys[peak_keep])
  edges$.tf_id <- unname(tf_lookup[edges$tf])
  edges$.peak_id <- unname(
    peak_lookup[toupper(.rc_pando_region_key(edges$region))]
  )
  edges <- edges[
    !is.na(edges$.tf_id) & nzchar(edges$.tf_id) &
      !is.na(edges$.peak_id) & nzchar(edges$.peak_id),
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
    celltype_value <- as.character(group_edges[[celltype_col]][[1L]])
    reference_units <- units[
      as.character(unit_meta[[celltype_col]]) == celltype_value
    ]
    if (!length(reference_units)) next
    for (target in unique(group_edges$target)) {
      selected <- group_edges[group_edges$target == target, , drop = FALSE]
      gene_id <- tolower(target)
      if (!gene_id %in% rownames(modifier) || !nrow(selected)) next
      tf_activity <- as.matrix(
        rna[selected$.tf_id, units, drop = FALSE]
      )
      peak_activity <- as.matrix(
        atac[selected$.peak_id, units, drop = FALSE]
      )
      edge_activity <- .rc_tf_peak_interaction(tf_activity, peak_activity)
      edge_deviation_reference <- .rc_edge_activity_deviation(
        edge_activity[, reference_units, drop = FALSE],
        min_scale = min_scale
      )
      edge_deviation <- edge_deviation_reference[
        , group_units, drop = FALSE
      ]
      weight <- abs(selected$estimate)
      weight[!is.finite(weight)] <- 0
      if (!any(weight > 0)) next
      weight <- weight / sum(weight)
      model_rsq <- suppressWarnings(as.numeric(selected$rsq))
      model_rsq <- model_rsq[is.finite(model_rsq)]
      reliability <- if (length(model_rsq)) {
        sqrt(min(max(stats::median(model_rsq), 0), 1))
      } else {
        0
      }
      signed_weight <- weight * sign(selected$estimate)
      value <- reliability * as.numeric(
        crossprod(signed_weight, edge_deviation)
      )
      modifier[gene_id, group_units] <- pmax(pmin(value, 1), -1)
    }
  }
  attr(modifier, "tf_target_overlap") <- intersect(
    unique(edges$tf), unique(edges$target)
  )
  modifier
}

.rc_integrate_regulatory_support <- function(
    rna_support, regulatory_modifier, alpha = 1) {
  rna_support <- as.matrix(rna_support)
  regulatory_modifier <- as.matrix(regulatory_modifier)
  if (!identical(dim(rna_support), dim(regulatory_modifier)) ||
      !identical(dimnames(rna_support), dimnames(regulatory_modifier))) {
    stop("RNA support and regulatory modifier matrices must align exactly.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0) {
    stop("`alpha` must be one finite non-negative number.", call. = FALSE)
  }
  C <- pmin(pmax(rna_support, 0), 1)
  R <- pmin(pmax(regulatory_modifier, -1), 1)
  multiplier <- 2^(alpha * R)
  numerator <- C * multiplier
  denominator <- 1 - C + numerator
  out <- numerator / denominator
  out[C <= 0] <- 0
  out[C >= 1] <- 1
  out[!is.finite(out)] <- NA_real_
  dimnames(out) <- dimnames(C)
  attr(out, "integration_formula") <- paste(
    "C_multiome = C_RNA * 2^(alpha * R_condition_TFxATAC) /",
    "(1 - C_RNA + C_RNA * 2^(alpha * R_condition_TFxATAC))"
  )
  attr(out, "score_semantics") <- paste(
    "zero-preserving bounded target-gene RNA support with a signed",
    "condition-effect TF-by-ATAC modifier on the support log-odds scale"
  )
  out
}

.rc_build_condition_pooled_layer1 <- function(
    metacell_object, meta_modules, gem, metacell_meta,
    sample_col = "sample_id", condition_col = "condition",
    celltype_col = "cell_type", rna_assay = "RNA", atac_assay = "ATAC",
    regulatory_alpha = 1,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE, BPPARAM = NULL) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  gpr_and_method <- match.arg(gpr_and_method)
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  counts <- .rc_get_assay_counts(metacell_object, rna_assay)
  full_library_size <- Matrix::colSums(counts)
  keep <- tolower(rownames(counts)) %in% gpr_genes
  rna_counts <- counts[keep, , drop = FALSE]
  rna_logcpm <- .rc_metacell_logcpm(
    rna_counts,
    library_size = full_library_size[colnames(rna_counts)]
  )
  rownames(rna_logcpm) <- tolower(rownames(rna_logcpm))
  if (anyDuplicated(rownames(rna_logcpm))) {
    stop("Duplicated GPR gene identifiers after case normalization.",
         call. = FALSE)
  }

  unit_meta <- metacell_meta
  id_col <- if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else {
    stop("Pooled metacell metadata lack metacell_id/pool_id.", call. = FALSE)
  }
  unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  unit_meta$unit_id <- unit_meta$pool_id
  unit_meta[[sample_col]] <- paste0(
    as.character(unit_meta[[condition_col]]), "__pooled"
  )
  unit_meta <- unit_meta[
    match(colnames(rna_logcpm), unit_meta$pool_id), , drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Pooled metacell metadata do not align with RNA counts.",
         call. = FALSE)
  }

  gene_rna_support <- rc_gene_score(
    rna_logcpm,
    mode = "absolute",
    half_saturation = gene_half_saturation
  )
  regulatory_edges <- meta_modules$tf_peak_gene_condition_effect
  if (!is.data.frame(regulatory_edges)) {
    stop(
      "Stage 3 lacks condition-effect TF-peak-target edges. Rerun GRN and ",
      "meta-module stages with the unified condition-aware Pando workflow.",
      call. = FALSE
    )
  }
  modifier <- .rc_condition_gene_regulatory_modifier(
    significant_edges = regulatory_edges,
    object = metacell_object,
    unit_meta = unit_meta,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    target_genes = rownames(gene_rna_support)
  )
  modifier <- modifier[
    rownames(gene_rna_support),
    colnames(gene_rna_support),
    drop = FALSE
  ]
  gene_multiome_support <- .rc_integrate_regulatory_support(
    gene_rna_support,
    modifier,
    alpha = regulatory_alpha
  )
  reaction_expression <- rc_reaction_capacity(
    parsed,
    gene_multiome_support,
    promiscuity_mode = "none",
    and_method = gpr_and_method,
    or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )

  list(
    schema_version = "regcompass_condition_multitask_layer1_v1",
    reaction_expression = reaction_expression,
    rna_metacell_logcpm = rna_logcpm,
    gene_support_rna = gene_rna_support,
    gene_regulatory_modifier = modifier,
    gene_support_multiome = gene_multiome_support,
    parsed_gpr = parsed,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, rownames(rna_logcpm)),
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "condition_only_metacell_with_posthoc_celltype",
    capacity_params = list(
      regulatory_alpha = regulatory_alpha,
      gene_half_saturation = gene_half_saturation,
      regulatory_mode =
        "condition_effect_coefficient_x_TF_RNA_x_peak_ATAC",
      promiscuity_mode = "none",
      and_method = gpr_and_method,
      or_method = "sum",
      parallel = parallel,
      bpparam_class = if (is.null(BPPARAM)) {
        "auto_or_sequential"
      } else if (identical(BPPARAM, FALSE)) {
        "sequential"
      } else {
        class(BPPARAM)[[1L]]
      }
    ),
    evidence_formula = paste(
      "condition GRN active edges define metabolic targets/core reactions;",
      "condition-effect coefficients weight TF RNA x peak ATAC activity;",
      "the signed modifier updates target-gene RNA support on the log-odds scale;",
      paste0("COMPASS-compatible ", gpr_and_method, " GPR-AND"),
      "and additive isozyme OR"
    ),
    evidence_inputs = c(
      "target_metabolic_gene_RNA",
      "regulatory_TF_RNA",
      "regulatory_peak_ATAC",
      "condition_level_Pando_coefficient",
      "universal_celltype_Pando_coefficient",
      "condition_effect_coefficient"
    )
  )
}
