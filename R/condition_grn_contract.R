# Direct integration of the canonical Pando condition-GRN contract.

.RC_PANDO_CONDITION_GRN_FIT_SCHEMA <- "pando_condition_grn_fit"

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

.rc_require_pando_condition_grn_fit <- function(fit) {
  if (!inherits(fit, "ConditionGRNFit") ||
      !identical(fit$schema_version, .RC_PANDO_CONDITION_GRN_FIT_SCHEMA)) {
    stop(
      "RegCompass requires the canonical pando_condition_grn_fit contract.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  fits <- Pando::condition_grn_fit(
    grn_object, network_name = "regcompass_condition_grn"
  )
  if (inherits(fits, "ConditionGRNFit")) fits <- list(fits)
  if (!is.list(fits) || !length(fits)) {
    stop("Pando did not return condition GRN fit contracts.", call. = FALSE)
  }
  invisible(lapply(fits, .rc_require_pando_condition_grn_fit))
  active_tol <- max(as.numeric(min_abs_estimate), 1e-8)
  rows <- list()
  universal <- list()
  index <- 1L
  for (fit in fits) {
    if (!identical(
          fit$fit_engine, "condition_sparse_within_cell_type_oof_refit"
        ) || !identical(
          fit$coefficient_scale,
          "equal_condition_within_variance_standardized_refit"
        )) {
      stop("Pando condition fit engine or coefficient scale is incompatible.",
           call. = FALSE)
    }
    edge <- as.data.frame(fit$edge_table, stringsAsFactors = FALSE)
    required_edge <- c("edge_id", "tf", "target", "region", "term")
    if (!all(required_edge %in% colnames(edge)) || anyDuplicated(edge$edge_id)) {
      stop("ConditionGRNFit edge dictionary is invalid.", call. = FALSE)
    }
    beta <- as.matrix(fit$beta_condition_std)
    matrices <- list(
      estimability = as.matrix(fit$estimability_mask),
      structural = as.matrix(fit$structural_candidate_mask),
      screening = as.matrix(fit$screening_mask),
      support = as.matrix(fit$support_mask),
      active = as.matrix(fit$active_mask)
    )
    expected <- c(nrow(edge), length(fit$condition_levels))
    invalid <- vapply(c(list(beta = beta), matrices), function(value) {
      !identical(dim(value), expected) ||
        !identical(rownames(value), edge$edge_id) ||
        !identical(colnames(value), fit$condition_levels)
    }, logical(1))
    if (any(invalid) || !all(vapply(matrices, is.logical, logical(1))) ||
        any(vapply(matrices, anyNA, logical(1))) ||
        any(!is.finite(beta[matrices$estimability])) ||
        any(!is.na(beta[!matrices$estimability])) ||
        any(matrices$support & !matrices$estimability) ||
        any(matrices$active & !matrices$estimability)) {
      stop("ConditionGRNFit coefficient and support matrices are invalid.",
           call. = FALSE)
    }
    transform <- as.data.frame(fit$predictor_transform)
    response <- as.data.frame(fit$response_transform)
    transform <- transform[match(edge$edge_id, transform$edge_id), , drop = FALSE]
    if (anyNA(transform$edge_id) ||
        any(!is.finite(transform$center)) ||
        any(!is.finite(transform$scale) | transform$scale <= 0)) {
      stop("ConditionGRNFit predictor transforms are invalid.", call. = FALSE)
    }
    response$target <- toupper(as.character(response$target))
    target_rsq <- as.matrix(fit$condition_rsq_train)
    pooled_oof <- fit$target_rsq_oof_pooled
    predictive <- fit$predictive_oof_available
    if (!identical(fit$projection_origin,
                   "outer_condition_stratified_cell_oof") ||
        !isTRUE(fit$projection_used_for_penalty) ||
        !identical(fit$full_fit_projection_used_for_penalty, FALSE) ||
        any(fit$oof_cell_coverage != 1) ||
        any(fit$oof_assignment_count != 1L)) {
      stop("ConditionGRNFit lacks complete outer-heldout projections.",
           call. = FALSE)
    }
    fit_cell_type <- as.character(fit$cell_type)
    for (condition in fit$condition_levels) {
      tab <- edge
      tab$tf <- toupper(as.character(tab$tf))
      tab$target <- toupper(as.character(tab$target))
      tab$region <- as.character(tab$region)
      tab$condition_estimate <- as.numeric(beta[, condition])
      tab$condition_effect <- tab$condition_estimate
      tab$effect_definition <- "absolute_condition_coefficient"
      tab$estimate <- tab$condition_estimate
      tab$corr <- NA_real_
      tab$eligible_in_condition <- matrices$estimability[, condition]
      tab$structural_candidate <- matrices$structural[, condition]
      tab$screened_in_condition <- matrices$screening[, condition]
      tab$selected_support <- matrices$support[, condition]
      tab$active_in_condition <- matrices$active[, condition]
      tab$predictor_center <- transform$center
      tab$predictor_scale <- transform$scale
      response_index <- match(tab$target, response$target)
      rsq_index <- match(tab$target, toupper(rownames(target_rsq)))
      if (anyNA(response_index) || anyNA(rsq_index)) {
        stop("ConditionGRNFit target diagnostics are misaligned.", call. = FALSE)
      }
      tab$response_center <- response$center[response_index]
      tab$response_scale <- response$scale[response_index]
      tab$rsq_train <- target_rsq[rsq_index, condition]
      tab$rsq_oof_pooled <- as.numeric(pooled_oof[rsq_index])
      tab$predictive_oof_available <- as.logical(predictive[rsq_index])
      tab$oof_reliability_available <-
        tab$predictive_oof_available & is.finite(tab$rsq_oof_pooled)
      tab$rsq <- tab$rsq_oof_pooled
      tab[[condition_col]] <- condition
      tab[[celltype_col]] <- fit_cell_type
      tab$group_id <- rc_make_stratum_id(
        tab[1L, c(condition_col, celltype_col), drop = FALSE],
        c(condition_col, celltype_col)
      )
      tab$fit_engine <- fit$fit_engine
      tab$coefficient_scale <- fit$coefficient_scale
      tab$coefficient_contract <- "absolute_condition_effects_only"
      tab <- tab[, c(
        "group_id", condition_col, celltype_col,
        setdiff(colnames(tab), c("group_id", condition_col, celltype_col))
      ), drop = FALSE]
      rows[[index]] <- tab
      index <- index + 1L
    }
    summary <- edge
    summary$estimate <- rowMeans(beta, na.rm = TRUE)
    summary$corr <- NA_real_
    summary[[celltype_col]] <- fit_cell_type
    summary$summary_only <- TRUE
    summary$coefficient_contract <- "absolute_condition_effects_only"
    universal[[length(universal) + 1L]] <- summary
  }
  all_edges <- do.call(rbind, rows)
  rownames(all_edges) <- NULL
  reliable <- all_edges$oof_reliability_available &
    is.finite(all_edges$rsq_oof_pooled) &
    all_edges$rsq_oof_pooled >= min_model_rsq
  active <- all_edges[
    is.finite(all_edges$condition_estimate) &
      abs(all_edges$condition_estimate) >= active_tol &
      all_edges$eligible_in_condition %in% TRUE & reliable,
    , drop = FALSE
  ]
  params <- methods::slot(methods::slot(grn_object, "grn"), "params")
  list(
    network_index = params$condition_network_index %||% data.frame(),
    fit_diagnostics = params$condition_fit_diagnostics %||% data.frame(),
    fit_contracts = fits,
    universal = do.call(rbind, universal),
    condition_all = all_edges,
    condition_active = active,
    condition_effect_all = all_edges,
    condition_effect_active = active,
    active_tol = active_tol,
    coefficient_contract = "absolute_condition_effects_only",
    pando_fit_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
  )
}

.rc_fit_condition_grns_by_cell_type <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    cell_type = NULL, rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list(
      candidate_screen = "motif_domain", tf_cor = 0.1, peak_cor = 0,
      alpha = 0.5, condition_mix = 0.5, condition_weight = "equal",
      nlambda = 50L, outer_nfolds = 5L, inner_nfolds = 5L,
      lambda_selection = "lambda.1se", scale = TRUE, parallel = FALSE
    ),
    min_abs_estimate = 0, min_model_rsq = 0.1,
    save_pando_objects = TRUE, BPPARAM = NULL,
    progress_monitor = NULL,
    species = c("auto", "human", "mouse")) {
  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  .rc_step_monitor_event(
    progress_monitor, "pando_contract_check",
    "validating Pando runtime and paired cell metadata", current = 5L,
    context = list(species = species)
  )
  .rc_validate_condition_celltype_metadata(
    object@meta.data, condition_col, celltype_col,
    require_multiple_conditions = TRUE
  )
  pando_runtime <- .rc_require_pando_hybrid_runtime()
  .rc_step_monitor_event(
    progress_monitor, "pando_runtime",
    "validated Pando fused native condition runtime", current = 5L,
    context = list(
      pando_version = pando_runtime$version,
      native_sparse_abi = pando_runtime$native_sparse_abi,
      target_engine = pando_runtime$target_engine_backend,
      inner_cv = pando_runtime$inner_cv_backend,
      refit = pando_runtime$refit_backend
    )
  )
  pando_infer_args <- modifyList(list(
    candidate_screen = "motif_domain", tf_cor = 0.1, peak_cor = 0,
    alpha = 0.5, condition_mix = 0.5, condition_weight = "equal",
    nlambda = 50L, outer_nfolds = 5L, inner_nfolds = 5L,
    lambda_selection = "lambda.1se", scale = TRUE, parallel = FALSE
  ), pando_infer_args)
  if (!identical(pando_infer_args$candidate_screen, "motif_domain") ||
      !identical(pando_infer_args$condition_weight, "equal") ||
      !isTRUE(pando_infer_args$scale)) {
    stop(
      "Condition comparability requires candidate_screen='motif_domain', condition_weight='equal', and scale=TRUE.",
      call. = FALSE
    )
  }
  .rc_step_monitor_event(
    progress_monitor, "condition_design",
    "condition-comparable fit controls validated", current = 5L,
    context = list(
      conditions = length(unique(as.character(
        object@meta.data[[condition_col]]
      ))),
      cell_types = length(unique(as.character(
        object@meta.data[[celltype_col]]
      ))),
      outer_folds = pando_infer_args$outer_nfolds,
      inner_folds = pando_infer_args$inner_nfolds,
      nlambda = pando_infer_args$nlambda,
      lambda_selection = pando_infer_args$lambda_selection
    )
  )
  if (is.null(pfm)) pfm <- .rc_default_pando_motifs()
  if (!"regions" %in% names(pando_initiate_args) ||
      is.null(pando_initiate_args$regions)) {
    pando_initiate_args$regions <- .rc_default_pando_regions(species)
  }
  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  target_upper <- intersect(toupper(rna_genes), toupper(metabolic_genes))
  target_genes <- rna_genes[toupper(rna_genes) %in% target_upper]
  if (!length(target_genes)) {
    stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)
  }
  .rc_step_monitor_event(
    progress_monitor, "target_selection",
    "resolved metabolic target genes shared by RNA and GEM", current = 6L,
    context = list(targets = length(target_genes))
  )
  filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Pando condition GRNs"
  )
  object <- filtered$object
  .rc_step_monitor_event(
    progress_monitor, "atac_feature_filter",
    "removed globally zero ATAC features", current = 6L,
    context = list(
      cells = ncol(object),
      removed_atac_features = filtered$n_removed %||% NA_integer_
    )
  )
  init <- list(object = object, peak_assay = atac_assay, rna_assay = rna_assay)
  init[names(pando_initiate_args)] <- NULL
  .rc_step_monitor_event(
    progress_monitor, "candidate_initialization",
    "initializing Pando regulatory candidate space", current = 7L,
    context = list(targets = length(target_genes))
  )
  grn <- do.call(Pando::initiate_grn, c(init, pando_initiate_args))
  .rc_step_monitor_event(
    progress_monitor, "candidate_initialization_complete",
    "Pando regulatory candidate space initialized", current = 7L,
    context = list(targets = length(target_genes))
  )
  pando_motif_args <- .rc_regcompass_motif_args(pando_motif_args)
  motif <- list(object = grn, pfm = pfm, genome = genome)
  motif[names(pando_motif_args)] <- NULL
  .rc_step_monitor_event(
    progress_monitor, "motif_mapping",
    "mapping binary peak-by-motif candidates", current = 8L
  )
  grn <- do.call(Pando::find_motifs, c(motif, pando_motif_args))
  .rc_step_monitor_event(
    progress_monitor, "motif_mapping_complete",
    "completed motif-to-peak and TF mapping", current = 8L
  )
  infer <- list(
    object = grn,
    cell_type_col = celltype_col,
    condition_col = condition_col,
    cell_type = cell_type,
    genes = target_genes,
    network_name = "regcompass_condition_grn",
    min_cells_per_condition = as.integer(min_cells),
    small_condition_action = "error",
    BPPARAM = if (identical(BPPARAM, FALSE)) NULL else BPPARAM
  )
  infer[names(pando_infer_args)] <- NULL
  .rc_step_monitor_event(
    progress_monitor, "nested_cv",
    "running fused per-target outer/inner CV, path, refit and validation",
    current = 9L,
    context = list(
      targets = length(target_genes),
      cell_types = length(unique(as.character(
        object@meta.data[[celltype_col]]
      ))),
      conditions = length(unique(as.character(
        object@meta.data[[condition_col]]
      ))),
      outer_folds = pando_infer_args$outer_nfolds,
      inner_folds = pando_infer_args$inner_nfolds,
      nlambda = pando_infer_args$nlambda,
      solver = "hybrid_gram_or_sparse_matrix_free",
      validation = "exact_sufficient_statistics",
      oof = "outer_selected_model_only"
    )
  )
  grn <- do.call(Pando::infer_condition_grn, c(infer, pando_infer_args))
  .rc_step_monitor_event(
    progress_monitor, "nested_cv_complete",
    "Pando fused target engine completed", current = 10L
  )
  .rc_step_monitor_event(
    progress_monitor, "contract_extraction",
    "extracting and validating ConditionGRNFit contracts", current = 10L
  )
  extracted <- .rc_extract_condition_grn_contract(
    grn, condition_col, celltype_col,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq
  )
  .rc_step_monitor_event(
    progress_monitor, "contract_extraction_complete",
    "validated and extracted ConditionGRNFit contracts", current = 10L,
    context = list(
      fitted_cell_types = length(extracted$fit_contracts),
      all_edges = nrow(extracted$condition_all),
      active_edges = nrow(extracted$condition_active)
    )
  )
  meta <- object@meta.data
  fitted_cell_types <- unique(vapply(
    extracted$fit_contracts, `[[`, character(1), "cell_type"
  ))
  expected <- unique(meta[
    as.character(meta[[celltype_col]]) %in% fitted_cell_types,
    c(condition_col, celltype_col), drop = FALSE
  ])
  expected$n_cells <- vapply(seq_len(nrow(expected)), function(i) {
    sum(as.character(meta[[condition_col]]) == expected[[condition_col]][[i]] &
      as.character(meta[[celltype_col]]) == expected[[celltype_col]][[i]])
  }, integer(1))
  expected$group_id <- rc_make_stratum_id(
    expected, c(condition_col, celltype_col)
  )
  expected$status <- "ok"
  expected$n_target_genes <- length(target_genes)
  expected$n_edges <- vapply(expected$group_id, function(id) {
    sum(extracted$condition_all$group_id == id)
  }, integer(1))
  expected$n_active_edges <- vapply(expected$group_id, function(id) {
    sum(extracted$condition_active$group_id == id)
  }, integer(1))
  expected$predictive_oof_available <- TRUE
  expected$oof_validation_level <-
    "outer_condition_stratified_heldout_cells"
  expected$grn_evidence_role <- "within_cell_type_multiomic_condition_grn"
  status <- expected[, c(
    "group_id", condition_col, celltype_col,
    setdiff(colnames(expected), c("group_id", condition_col, celltype_col))
  ), drop = FALSE]
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_step_monitor_event(
    progress_monitor, "grn_artifacts",
    "writing Stage 1 GRN tables and fit contracts", current = 11L,
    context = list(groups = nrow(status))
  )
  .rc_write_tsv_gz(status, file.path(outdir, "pando_group_status.tsv.gz"))
  .rc_write_tsv_gz(extracted$condition_all,
    file.path(outdir, "pando_tf_peak_gene_condition_all.tsv.gz"))
  .rc_write_tsv_gz(extracted$condition_active,
    file.path(outdir, "pando_tf_peak_gene_condition_active.tsv.gz"))
  saveRDS(extracted$fit_contracts,
    file.path(outdir, "pando_condition_grn_fits.rds"))
  if (isTRUE(save_pando_objects)) {
    dir.create(file.path(outdir, "pando_objects"), recursive = TRUE,
               showWarnings = FALSE)
    saveRDS(grn, file.path(outdir, "pando_objects", "condition_grn_fit.rds"))
  }
  selected_cells <- rownames(meta)[
    as.character(meta[[celltype_col]]) %in% fitted_cell_types
  ]
  answer <- list(
    schema_version = "regcompass_condition_grn_fit_v1",
    analysis_mode = "condition_grn",
    condition_coefficients_calculated = TRUE,
    pando_fit_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    pando_installed_version = as.character(utils::packageVersion("Pando")),
    pando_file_fingerprint = .rc_installed_package_file_fingerprint("Pando"),
    pando_grn_data = grn,
    paired_cell_ids = selected_cells,
    paired_cell_metadata = data.frame(
      cell_id = selected_cells,
      condition = as.character(meta[selected_cells, condition_col]),
      cell_type = as.character(meta[selected_cells, celltype_col]),
      stringsAsFactors = FALSE
    ),
    target_metabolic_genes = target_genes,
    condition_fit_status = status,
    pando_network_index = extracted$network_index,
    pando_fit_diagnostics = extracted$fit_diagnostics,
    condition_grn_fits = extracted$fit_contracts,
    tf_peak_gene_universal = extracted$universal,
    tf_peak_gene_condition_all = extracted$condition_all,
    tf_peak_gene_condition = extracted$condition_active,
    tf_peak_gene_condition_effect_all = extracted$condition_all,
    tf_peak_gene_condition_effect = extracted$condition_active,
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF across conditions",
      grn_fit = "condition-aware Pando nested OOF within broad cell type",
      condition_effect =
        "absolute condition coefficient on the shared equal-condition coordinate",
      coefficient_contract = "absolute_condition_effects_only",
      penalty_regulatory_evidence =
        "outer-heldout cell-first TF-by-ATAC projection"
    ),
    group_cols = c(condition_col, celltype_col)
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
