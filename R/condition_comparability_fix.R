# Condition-comparison safeguards loaded after condition_grn_contract.R.

.rc_condition_fit_comparison_mask <- function(fit) {
  beta <- as.matrix(fit$beta)
  eligibility <- as.matrix(fit$eligibility_mask)
  if (!is.logical(eligibility) || anyNA(eligibility) ||
      !identical(dim(eligibility), dim(beta)) ||
      !identical(dimnames(eligibility), dimnames(beta))) {
    stop("ConditionGRNFit eligibility mask is invalid.", call. = FALSE)
  }
  comparison <- fit$comparison_mask
  if (is.null(comparison)) {
    if (!fit$reference_condition %in% colnames(eligibility)) {
      stop("ConditionGRNFit reference condition is absent.", call. = FALSE)
    }
    reference_eligible <- eligibility[, fit$reference_condition]
    comparison <- eligibility & matrix(
      reference_eligible,
      nrow = nrow(eligibility),
      ncol = ncol(eligibility)
    )
    dimnames(comparison) <- dimnames(eligibility)
  }
  comparison <- as.matrix(comparison)
  if (!is.logical(comparison) || anyNA(comparison) ||
      !identical(dim(comparison), dim(beta)) ||
      !identical(dimnames(comparison), dimnames(beta))) {
    stop("ConditionGRNFit comparison mask is invalid.", call. = FALSE)
  }
  expected <- eligibility & matrix(
    eligibility[, fit$reference_condition],
    nrow = nrow(eligibility),
    ncol = ncol(eligibility)
  )
  dimnames(expected) <- dimnames(eligibility)
  if (!identical(comparison, expected)) {
    stop(
      "ConditionGRNFit comparison mask is inconsistent with eligibility.",
      call. = FALSE
    )
  }
  comparison
}

.rc_apply_condition_comparison_semantics <- function(
    extracted, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  if (!is.list(extracted) || !is.list(extracted$fit_contracts)) {
    stop("Condition GRN extraction lacks fit contracts.", call. = FALSE)
  }
  lookup_rows <- lapply(extracted$fit_contracts, function(fit) {
    comparison <- .rc_condition_fit_comparison_mask(fit)
    do.call(rbind, lapply(colnames(comparison), function(condition) {
      data.frame(
        edge_id = rownames(comparison),
        condition = condition,
        cell_type = fit$cell_type,
        comparable_to_reference = as.logical(comparison[, condition]),
        stringsAsFactors = FALSE
      )
    }))
  })
  lookup <- do.call(rbind, lookup_rows)
  if (!nrow(lookup) || anyDuplicated(paste(
    lookup$cell_type, lookup$condition, lookup$edge_id, sep = "\001"
  ))) {
    stop("Condition comparison-mask lookup is empty or duplicated.",
         call. = FALSE)
  }
  annotate <- function(table) {
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    if (!nrow(table)) {
      table$comparable_to_reference <- logical()
      return(table)
    }
    required <- c("edge_id", condition_col, celltype_col)
    missing <- setdiff(required, colnames(table))
    if (length(missing)) {
      stop(
        "Condition edge table lacks comparison keys: ",
        paste(missing, collapse = ", "), call. = FALSE
      )
    }
    table_key <- paste(
      as.character(table[[celltype_col]]),
      as.character(table[[condition_col]]),
      as.character(table$edge_id), sep = "\001"
    )
    lookup_key <- paste(
      lookup$cell_type, lookup$condition, lookup$edge_id, sep = "\001"
    )
    index <- match(table_key, lookup_key)
    if (anyNA(index)) {
      stop("Condition edge rows do not align with comparison masks.",
           call. = FALSE)
    }
    table$comparable_to_reference <-
      as.logical(lookup$comparable_to_reference[index])
    table
  }

  condition_all <- annotate(extracted$condition_all)
  effect_all <- annotate(extracted$condition_effect_all)
  active_tol <- max(
    suppressWarnings(as.numeric(min_abs_estimate)), 1e-8, na.rm = TRUE
  )
  eligible <- if ("eligible_in_condition" %in% colnames(condition_all)) {
    condition_all$eligible_in_condition %in% TRUE
  } else {
    rep(TRUE, nrow(condition_all))
  }
  if (!"sample_blocked_oof_available" %in% colnames(condition_all) ||
      !"sample_blocked_oof_available" %in% colnames(effect_all)) {
    stop(
      "Condition edge tables lack sample-blocked OOF availability.",
      call. = FALSE
    )
  }
  condition_oof_available <- if (
      "confirmatory_oof_available" %in% colnames(condition_all)) {
    condition_all$confirmatory_oof_available
  } else {
    condition_all$sample_blocked_oof_available
  }
  effect_oof_available <- if (
      "confirmatory_oof_available" %in% colnames(effect_all)) {
    effect_all$confirmatory_oof_available
  } else {
    effect_all$sample_blocked_oof_available
  }
  condition_reliable_or_unavailable <-
    !condition_oof_available |
    (
      is.finite(condition_all$rsq) &
      condition_all$rsq >= min_model_rsq
    )
  effect_reliable_or_unavailable <-
    !effect_oof_available |
    (
      is.finite(effect_all$rsq) &
      effect_all$rsq >= min_model_rsq
    )
  condition_active <- condition_all[
    eligible &
      is.finite(condition_all$condition_estimate) &
      abs(condition_all$condition_estimate) >= active_tol &
      condition_reliable_or_unavailable,
    , drop = FALSE
  ]
  effect_active <- effect_all[
    effect_all$comparable_to_reference %in% TRUE &
      is.finite(effect_all$condition_effect) &
      abs(effect_all$condition_effect) >= active_tol &
      effect_reliable_or_unavailable,
    , drop = FALSE
  ]

  extracted$condition_all <- condition_all
  extracted$condition_active <- condition_active
  extracted$condition_effect_all <- effect_all
  extracted$condition_effect_active <- effect_active
  extracted$comparison_policy <- paste(
    "condition effects require eligibility in both the condition and",
    "the explicit reference condition"
  )
  extracted$reliability_policy <- paste(
    "min_model_rsq is applied only where sample-blocked OOF is estimable;",
    "single-sample conditions retain exploratory coefficients but receive",
    "zero regulatory reliability in Layer 1"
  )
  extracted
}

.rc_extract_condition_grn_contract_without_comparison_guard <-
  .rc_extract_condition_grn_contract

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  extracted <- .rc_extract_condition_grn_contract_without_comparison_guard(
    grn_object = grn_object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq
  )
  .rc_apply_condition_comparison_semantics(
    extracted = extracted,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq
  )
}

.rc_run_condition_single_cell_grns_without_safe_defaults <-
  .rc_run_condition_single_cell_grns

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
  if (!is.list(pando_infer_args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  if (is.null(pando_infer_args$candidate_screen)) {
    pando_infer_args$candidate_screen <- "motif_domain"
  }
  .rc_run_condition_single_cell_grns_without_safe_defaults(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    min_cells = min_cells,
    pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args,
    pando_infer_args = pando_infer_args,
    padj_threshold = padj_threshold,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq,
    require_padj = require_padj,
    save_pando_objects = save_pando_objects,
    BPPARAM = BPPARAM,
    on_group_error = on_group_error,
    species = species
  )
}

.rc_default_pando_regions_without_genome_guard <- .rc_default_pando_regions

.rc_default_pando_regions <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (identical(species, "mouse")) {
    stop(
      paste(
        "No mouse-coordinate regulatory-region set is bundled by the required",
        "Pando fork. Supply `pando_initiate_args$regions` as a GRanges object",
        "whose genome build matches both the ATAC peaks and `genome`",
        "(for example mm10 or mm39). The hg38 conserved-element set is not",
        "valid for mouse input."
      ),
      call. = FALSE
    )
  }
  .rc_default_pando_regions_without_genome_guard(species = species)
}
