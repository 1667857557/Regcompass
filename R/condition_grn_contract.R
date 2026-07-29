# Authoritative Pando ConditionGRNFit integration.

.rc_reference_contrast <- function(beta, reference_condition) {
  beta <- as.matrix(beta)
  if (!is.character(reference_condition) ||
      length(reference_condition) != 1L ||
      !reference_condition %in% colnames(beta)) {
    stop("Reference condition is absent from the coefficient matrix.",
         call. = FALSE)
  }
  sweep(beta, 1L, beta[, reference_condition], "-")
}

.rc_installed_package_file_fingerprint <- function(package) {
  root <- system.file(package = package)
  if (!nzchar(root) || !dir.exists(root)) return(NA_character_)
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) return(NA_character_)
  relative <- substring(files, nchar(root) + 2L)
  md5 <- unname(as.character(tools::md5sum(files)))
  .rc_condition_metacell_md5(data.frame(
    file = relative, md5 = md5, stringsAsFactors = FALSE
  ))
}

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  fits <- Pando::condition_grn_fit(
    grn_object, network_name = "regcompass_condition_grn"
  )
  if (inherits(fits, "ConditionGRNFit")) {
    fits <- list(fits)
  }
  if (!is.list(fits) || !length(fits) ||
      !all(vapply(fits, inherits, logical(1), "ConditionGRNFit"))) {
    stop("Pando did not return complete ConditionGRNFit contracts.",
         call. = FALSE)
  }
  active_tol <- max(
    suppressWarnings(as.numeric(min_abs_estimate)), 1e-8, na.rm = TRUE
  )
  condition_rows <- list()
  universal_rows <- list()
  row_index <- 1L
  for (fit in fits) {
    if (!identical(fit$schema_version, "pando_condition_grn_fit_v4") ||
        !identical(
          fit$fit_engine, "condition_sparse_common_scale_refit"
        ) ||
        !identical(
          fit$coefficient_scale,
          "pooled_cell_type_standardized_refit"
        )) {
      stop(
        "RegCompass requires Pando's pooled-scale condition-sparse ",
        "ConditionGRNFit v4.", call. = FALSE
      )
    }
    edge <- as.data.frame(fit$edge_table, stringsAsFactors = FALSE)
    transform <- as.data.frame(
      fit$predictor_transform, stringsAsFactors = FALSE
    )
    response <- as.data.frame(
      fit$response_transform, stringsAsFactors = FALSE
    )
    required_edge <- c("edge_id", "tf", "target", "region", "term")
    if (!all(required_edge %in% colnames(edge)) ||
        anyDuplicated(edge$edge_id)) {
      stop("ConditionGRNFit edge dictionary is invalid.", call. = FALSE)
    }
    if (!all(c("edge_id", "center", "scale") %in% colnames(transform)) ||
        anyDuplicated(transform$edge_id)) {
      stop("ConditionGRNFit predictor transform is invalid.",
           call. = FALSE)
    }
    transform <- transform[
      match(edge$edge_id, transform$edge_id), , drop = FALSE
    ]
    if (anyNA(transform$edge_id) ||
        any(!is.finite(transform$center)) ||
        any(!is.finite(transform$scale) | transform$scale <= 0)) {
      stop("ConditionGRNFit predictor transform is incomplete.",
           call. = FALSE)
    }
    beta <- as.matrix(fit$beta_condition_std)
    contrast <- as.matrix(fit$contrast)
    mask <- as.matrix(fit$estimability_mask)
    structural <- as.matrix(fit$structural_candidate_mask)
    screening <- as.matrix(fit$screening_mask)
    support <- as.matrix(fit$support_mask)
    active <- as.matrix(fit$active_mask)
    comparison <- as.matrix(fit$comparison_mask)
    if (!identical(dim(beta), c(nrow(edge), length(fit$condition_levels))) ||
        !identical(dim(contrast), dim(beta)) ||
        !identical(dim(mask), dim(beta)) ||
        !identical(colnames(beta), fit$condition_levels) ||
        !identical(colnames(contrast), fit$condition_levels) ||
        !identical(colnames(mask), fit$condition_levels) ||
        !identical(rownames(beta), edge$edge_id) ||
        !identical(rownames(contrast), edge$edge_id) ||
        !identical(rownames(mask), edge$edge_id) ||
        any(!is.finite(beta[mask])) || any(!is.na(beta[!mask])) ||
        !is.logical(mask) || anyNA(mask) ||
        !fit$reference_condition %in% colnames(beta)) {
      stop("ConditionGRNFit coefficient matrices are misaligned.",
           call. = FALSE)
    }
    masks <- list(
      structural_candidate_mask = structural,
      screening_mask = screening,
      support_mask = support,
      active_mask = active,
      comparison_mask = comparison
    )
    invalid_mask <- vapply(masks, function(value) {
      !is.logical(value) || anyNA(value) ||
        !identical(dim(value), dim(beta)) ||
        !identical(dimnames(value), dimnames(beta))
    }, logical(1))
    if (any(invalid_mask) || any(support & !mask) || any(active & !mask)) {
      stop("ConditionGRNFit v4 support masks are invalid.", call. = FALSE)
    }
    beta_reference <- beta[, fit$reference_condition]
    expected_contrast <- .rc_reference_contrast(
      beta, fit$reference_condition
    )
    if (!isTRUE(all.equal(
      unname(contrast), unname(expected_contrast), tolerance = 1e-10
    ))) {
      stop("ConditionGRNFit reference contrasts are inconsistent.",
           call. = FALSE)
    }
    target_rsq <- as.matrix(fit$condition_rsq_train)
    target_rsq_oof <- suppressWarnings(as.numeric(
      fit$target_rsq_oof_pooled
    ))
    names(target_rsq_oof) <- names(fit$target_rsq_oof_pooled)
    sample_blocked_oof_available <- fit$sample_blocked_oof_available
    if (!is.logical(sample_blocked_oof_available) ||
        anyNA(sample_blocked_oof_available) ||
        is.null(names(sample_blocked_oof_available))) {
      stop(
        "ConditionGRNFit v4 lacks valid sample-blocked OOF availability.",
        call. = FALSE
      )
    }
    if (!all(c("target", "center", "scale") %in% colnames(response))) {
      stop("ConditionGRNFit target transforms are incomplete.",
           call. = FALSE)
    }
    response$target <- toupper(trimws(as.character(response$target)))
    if (anyDuplicated(response$target) ||
        any(!is.finite(response$center)) ||
        any(!is.finite(response$scale) | response$scale <= 0) ||
        is.null(rownames(target_rsq)) ||
        !identical(colnames(target_rsq), fit$condition_levels) ||
        !identical(names(target_rsq_oof), rownames(target_rsq)) ||
        !identical(
          names(sample_blocked_oof_available), rownames(target_rsq)
        )) {
      stop("ConditionGRNFit target transforms or diagnostics are invalid.",
           call. = FALSE)
    }
    for (condition in fit$condition_levels) {
      tab <- edge
      tab$tf <- toupper(trimws(as.character(tab$tf)))
      tab$target <- toupper(trimws(as.character(tab$target)))
      tab$region <- trimws(as.character(tab$region))
      tab$condition_estimate <- as.numeric(beta[, condition])
      tab$reference_estimate <- as.numeric(beta_reference)
      tab$condition_effect <- as.numeric(contrast[, condition])
      tab$estimate <- tab$condition_estimate
      tab$corr <- NA_real_
      tab$eligible_in_condition <- as.logical(mask[, condition])
      tab$structural_candidate <- as.logical(structural[, condition])
      tab$screened_in_condition <- as.logical(screening[, condition])
      tab$selected_support <- as.logical(support[, condition])
      tab$active_in_condition <- as.logical(active[, condition])
      tab$comparable_to_reference <- as.logical(comparison[, condition])
      tab$reference_condition <- fit$reference_condition
      tab$predictor_center <- transform$center
      tab$predictor_scale <- transform$scale
      response_index <- match(tab$target, response$target)
      rsq_index <- match(tab$target, toupper(rownames(target_rsq)))
      if (anyNA(response_index) || anyNA(rsq_index)) {
        stop(
          "ConditionGRNFit edge targets do not align with target transforms ",
          "and diagnostics.", call. = FALSE
        )
      }
      tab$response_center <- response$center[response_index]
      tab$response_scale <- response$scale[response_index]
      tab$rsq_train <- target_rsq[rsq_index, condition]
      tab$rsq_oof_pooled <- target_rsq_oof[rsq_index]
      tab$sample_blocked_oof_available <-
        sample_blocked_oof_available[rsq_index]
      tab$response_independent_candidate_graph <- identical(
        fit$candidate_screen, "motif_domain"
      )
      tab$confirmatory_oof_available <-
        tab$sample_blocked_oof_available &
        tab$response_independent_candidate_graph &
        is.finite(tab$rsq_oof_pooled)
      tab$reliability_status <- ifelse(
        tab$confirmatory_oof_available,
        "confirmatory_sample_blocked_oof",
        ifelse(
          !tab$sample_blocked_oof_available,
          "unavailable_insufficient_biological_samples",
          "exploratory_response_dependent_screen_or_nonfinite_oof"
        )
      )
      tab$rsq <- tab$rsq_oof_pooled
      tab[[condition_col]] <- condition
      tab[[celltype_col]] <- fit$cell_type
      group_values <- tab[1L, c(condition_col, celltype_col), drop = FALSE]
      tab$group_id <- rc_make_stratum_id(
        group_values, c(condition_col, celltype_col)
      )
      tab$fit_engine <- fit$fit_engine
      tab$coefficient_scale <- fit$coefficient_scale
      tab <- tab[, c(
        "group_id", condition_col, celltype_col,
        setdiff(colnames(tab), c("group_id", condition_col, celltype_col))
      ), drop = FALSE]
      condition_rows[[row_index]] <- tab
      row_index <- row_index + 1L
    }
    summary <- edge
    summary$estimate <- rowMeans(beta, na.rm = TRUE)
    summary$corr <- NA_real_
    summary[[celltype_col]] <- fit$cell_type
    summary$reference_condition <- fit$reference_condition
    summary$summary_only <- TRUE
    universal_rows[[length(universal_rows) + 1L]] <- summary
  }
  all_edges <- do.call(rbind, condition_rows)
  rownames(all_edges) <- NULL
  reliable_or_unavailable <-
    !all_edges$confirmatory_oof_available |
    (
      is.finite(all_edges$rsq_oof_pooled) &
      all_edges$rsq_oof_pooled >= min_model_rsq
    )
  condition_active <- all_edges[
    is.finite(all_edges$condition_estimate) &
      abs(all_edges$condition_estimate) >= active_tol &
      all_edges$eligible_in_condition %in% TRUE &
      reliable_or_unavailable,
    , drop = FALSE
  ]
  effect_all <- all_edges
  effect_all$estimate <- effect_all$condition_effect
  effect_reliable_or_unavailable <-
    !effect_all$confirmatory_oof_available |
    (
      is.finite(effect_all$rsq_oof_pooled) &
      effect_all$rsq_oof_pooled >= min_model_rsq
    )
  effect_active <- effect_all[
    is.finite(effect_all$condition_effect) &
      abs(effect_all$condition_effect) >= active_tol &
      effect_all$comparable_to_reference %in% TRUE &
      effect_reliable_or_unavailable,
    , drop = FALSE
  ]
  params <- methods::slot(methods::slot(grn_object, "grn"), "params")
  list(
    network_index = params$condition_network_index %||% data.frame(),
    fit_diagnostics = params$condition_fit_diagnostics %||% data.frame(),
    fit_contracts = fits,
    universal = do.call(rbind, universal_rows),
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
      method = "shared_baseline_condition_sparse",
      candidate_screen = "motif_domain",
      tf_cor = 0.1,
      peak_cor = 0,
      alpha = 0.5,
      condition_mix = 0.5,
      condition_weight = "equal",
      cv_block_col = "sample_id",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE,
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
  .rc_validate_pando_evidence_filters(
    padj_threshold = padj_threshold,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq,
    require_padj = require_padj
  )
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
      "Installed Pando lacks `infer_condition_grn`; install Pando_regcompass >= 1.4.0.",
      call. = FALSE
    )
  }
  if (!exists("condition_grn_fit", envir = asNamespace("Pando"),
              inherits = FALSE)) {
    stop(
      "Installed Pando lacks the ConditionGRNFit v4 accessor.",
      call. = FALSE
    )
  }
  pando_infer_args <- utils::modifyList(
    list(
      method = "shared_baseline_condition_sparse",
      candidate_screen = "motif_domain",
      tf_cor = 0.1,
      peak_cor = 0,
      alpha = 0.5,
      condition_mix = 0.5,
      condition_weight = "equal",
      cv_block_col = "sample_id",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE,
      parallel = FALSE
    ),
    pando_infer_args
  )
  required_infer <- list(
    method = "shared_baseline_condition_sparse",
    condition_weight = "equal",
    scale = TRUE
  )
  invalid_infer <- names(required_infer)[!vapply(
    names(required_infer),
    function(name) {
      if (identical(name, "scale")) {
        return(isTRUE(pando_infer_args[[name]]))
      }
      identical(
        as.character(pando_infer_args[[name]]),
        as.character(required_infer[[name]])
      )
    },
    logical(1)
  )]
  if (length(invalid_infer)) {
    stop(
      "RegCompass condition comparability requires Pando settings: ",
      paste(
        paste0(names(required_infer), "=", unlist(required_infer)),
        collapse = ", "
      ),
      ". Incompatible fields: ", paste(invalid_infer, collapse = ", "),
      call. = FALSE
    )
  }
  if (isTRUE(require_padj)) {
    warning(
      "`require_padj = TRUE` is ignored by the regularized condition solver; ",
      "active coefficients and target-level R-squared are used.",
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

  extracted <- .rc_extract_condition_grn_contract(
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
  cv_block_col <- as.character(
    pando_infer_args$cv_block_col %||% ""
  )
  status$n_biological_samples <- vapply(
    seq_len(nrow(status)), function(i) {
      if (!nzchar(cv_block_col) || !cv_block_col %in% colnames(meta)) {
        return(NA_integer_)
      }
      rows <-
        as.character(meta[[condition_col]]) ==
          as.character(status[[condition_col]][[i]]) &
        as.character(meta[[celltype_col]]) ==
          as.character(status[[celltype_col]][[i]])
      length(unique(as.character(meta[[cv_block_col]][rows])))
    }, integer(1)
  )
  status$sample_blocked_oof_available <- vapply(
    seq_len(nrow(status)), function(i) {
      same_celltype <- as.character(status[[celltype_col]]) ==
        as.character(status[[celltype_col]][[i]])
      counts <- status$n_biological_samples[same_celltype]
      length(counts) > 0L && all(!is.na(counts) & counts >= 2L)
    }, logical(1)
  )
  status$response_independent_candidate_graph <- identical(
    pando_infer_args$candidate_screen, "motif_domain"
  )
  status$confirmatory_oof_available <-
    status$sample_blocked_oof_available &
    status$response_independent_candidate_graph
  status$grn_evidence_role <- ifelse(
    status$confirmatory_oof_available,
    "confirmatory_sample_blocked",
    "exploratory_zero_layer1_reliability"
  )
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
  predictor_transforms <- do.call(rbind, lapply(
    extracted$fit_contracts,
    function(fit) {
      out <- merge(
        fit$edge_table,
        fit$predictor_transform,
        by = "edge_id",
        all.x = TRUE,
        sort = FALSE
      )
      out$cell_type <- fit$cell_type
      out$reference_condition <- fit$reference_condition
      out$coefficient_scale <- fit$coefficient_scale
      out
    }
  ))
  .rc_mm_write_tsv_gz(
    predictor_transforms,
    file.path(outdir, "pando_edge_predictor_transforms.tsv.gz")
  )
  saveRDS(
    extracted$fit_contracts,
    file.path(outdir, "pando_condition_grn_fits.rds")
  )
  if (isTRUE(save_pando_objects)) {
    saveRDS(grn, file.path(outdir, "pando_objects",
                           "condition_grn_fit_v4.rds"))
  }

  tf_metabolic_target_overlap <- intersect(
    unique(toupper(extracted$condition_effect_active$tf)),
    unique(toupper(target_genes))
  )
  answer <- list(
    schema_version = "regcompass_condition_grn_fit_v4",
    pando_installed_version = pando_install$version,
    pando_installation = pando_install,
    pando_file_fingerprint =
      .rc_installed_package_file_fingerprint("Pando"),
    pando_grn_data = grn,
    paired_cell_ids = colnames(object),
    paired_cell_metadata = data.frame(
      cell_id = colnames(object),
      condition = as.character(object@meta.data[
        colnames(object), condition_col
      ]),
      cell_type = as.character(object@meta.data[
        colnames(object), celltype_col
      ]),
      sample_id = if ("sample_id" %in% colnames(object@meta.data)) {
        as.character(object@meta.data[colnames(object), "sample_id"])
      } else {
        NA_character_
      },
      cv_block = if (nzchar(cv_block_col) &&
                     cv_block_col %in% colnames(object@meta.data)) {
        as.character(object@meta.data[colnames(object), cv_block_col])
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    ),
    cv_block_col = if (nzchar(cv_block_col)) cv_block_col else NULL,
    assay_contract = list(
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      rna_layer = "data",
      atac_layer = "data",
      paired_cell_order = colnames(object)
    ),
    target_metabolic_genes = target_genes,
    condition_fit_status = status,
    sample_status = status,
    pando_network_index = extracted$network_index,
    pando_fit_diagnostics = extracted$fit_diagnostics,
    condition_grn_fits = extracted$fit_contracts,
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
        "one Pando fit contract per cell type with one TF-peak-target",
        "edge dictionary, condition-sparse selection, and common-metric refit"
      ),
      universal_coefficient = "visualization-only row mean; never used as a contrast baseline",
      reference_condition = unique(vapply(
        extracted$fit_contracts, `[[`, character(1), "reference_condition"
      )),
      condition_effect = "condition coefficient minus explicit reference-condition coefficient",
      coefficient_scale =
        "pooled cell-type TF-by-ATAC edge and target standardization",
      core_reaction_evidence =
        "active condition-level TF-peak-target coefficients",
      penalty_regulatory_evidence = paste(
        "absolute condition coefficients applied to Pando cell-first",
        "TF RNA by peak ATAC projections, then aggregated by",
        "SuperCell membership"
      ),
      pando_motifs = motif_policy,
      pando_regions = region_policy,
      pando_padj = paste(
        "not applicable to the regularized condition solver;",
        "active coefficients and target-level R-squared are used"
      ),
      pando_peak_cor = pando_infer_args$peak_cor,
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
