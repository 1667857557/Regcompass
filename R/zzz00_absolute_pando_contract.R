# Consume absolute condition effects without a reference-condition contrast.

.rc_reference_contract_fields <- c(
  "reference_condition", "reference_estimate", "reference_beta",
  "contrast", "comparison_mask", "comparable_to_reference",
  "contrast_formula", "comparison_mask_formula"
)

.rc_strip_reference_contract <- function(x) {
  if (is.data.frame(x)) {
    return(x[, setdiff(colnames(x), .rc_reference_contract_fields),
             drop = FALSE])
  }
  if (!is.list(x)) return(x)
  x[intersect(names(x), .rc_reference_contract_fields)] <- NULL
  if (is.data.frame(x$response_transform)) {
    x$response_transform <- .rc_strip_reference_contract(
      x$response_transform
    )
  }
  if (is.list(x$direction_semantics)) {
    x$direction_semantics$pairwise <- NULL
  }
  x$coefficient_contract <- "absolute_condition_effects_only"
  x
}

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  fits <- Pando::condition_grn_fit(
    grn_object, network_name = "regcompass_condition_grn"
  )
  if (inherits(fits, "ConditionGRNFit")) fits <- list(fits)
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
  fits_internal <- vector("list", length(fits))

  for (fit_index in seq_along(fits)) {
    fit <- .rc_strip_reference_contract(fits[[fit_index]])
    if (!identical(fit$schema_version, "pando_condition_grn_fit_v5") ||
        !identical(
          fit$fit_engine, "condition_sparse_within_cell_type_oof_refit"
        ) ||
        !identical(
          fit$coefficient_scale,
          "equal_condition_within_variance_standardized_refit"
        )) {
      stop(
        "RegCompass requires one comparable ConditionGRNFit v5 per cell type.",
        call. = FALSE
      )
    }
    fit_cell_type <- trimws(as.character(fit$cell_type %||% ""))
    if (length(fit_cell_type) != 1L || is.na(fit_cell_type) ||
        !nzchar(fit_cell_type)) {
      stop("ConditionGRNFit must identify one fitted cell type.",
           call. = FALSE)
    }
    edge <- as.data.frame(fit$edge_table, stringsAsFactors = FALSE)
    required_edge <- c("edge_id", "tf", "target", "region", "term")
    if (!all(required_edge %in% colnames(edge)) ||
        anyDuplicated(edge$edge_id)) {
      stop("ConditionGRNFit edge dictionary is invalid.", call. = FALSE)
    }
    beta <- as.matrix(fit$beta_condition_std %||% fit$beta_condition)
    estimability <- as.matrix(fit$estimability_mask)
    structural <- as.matrix(
      fit$structural_candidate_mask %||% fit$topology_mask
    )
    screening <- as.matrix(fit$screening_mask)
    support <- as.matrix(fit$support_mask)
    active <- as.matrix(fit$active_mask)
    expected_dim <- c(nrow(edge), length(fit$condition_levels))
    matrices <- list(
      beta = beta,
      estimability = estimability,
      structural = structural,
      screening = screening,
      support = support,
      active = active
    )
    invalid <- vapply(matrices, function(value) {
      !identical(dim(value), expected_dim) ||
        !identical(rownames(value), edge$edge_id) ||
        !identical(colnames(value), fit$condition_levels)
    }, logical(1))
    if (any(invalid) || !is.logical(estimability) || anyNA(estimability) ||
        any(!is.finite(beta[estimability])) || any(!is.na(beta[!estimability])) ||
        !all(vapply(matrices[-1L], is.logical, logical(1))) ||
        any(vapply(matrices[-1L], anyNA, logical(1))) ||
        any(support & !estimability) || any(active & !estimability)) {
      stop("ConditionGRNFit absolute coefficient matrices are invalid.",
           call. = FALSE)
    }
    transform <- as.data.frame(
      fit$predictor_transform, stringsAsFactors = FALSE
    )
    response <- as.data.frame(
      fit$response_transform, stringsAsFactors = FALSE
    )
    if (!all(c("edge_id", "center", "scale") %in% colnames(transform)) ||
        !all(c("target", "center", "scale") %in% colnames(response))) {
      stop("ConditionGRNFit transforms are incomplete.", call. = FALSE)
    }
    transform <- transform[match(edge$edge_id, transform$edge_id), , drop = FALSE]
    if (anyNA(transform$edge_id) || any(!is.finite(transform$center)) ||
        any(!is.finite(transform$scale) | transform$scale <= 0)) {
      stop("ConditionGRNFit predictor transform is invalid.", call. = FALSE)
    }
    response$target <- toupper(trimws(as.character(response$target)))
    if (anyDuplicated(response$target) || any(!is.finite(response$center)) ||
        any(!is.finite(response$scale) | response$scale <= 0)) {
      stop("ConditionGRNFit response transform is invalid.", call. = FALSE)
    }
    target_rsq <- as.matrix(fit$condition_rsq_train)
    pooled_oof <- suppressWarnings(as.numeric(fit$target_rsq_oof_pooled))
    names(pooled_oof) <- names(fit$target_rsq_oof_pooled)
    predictive <- fit$predictive_oof_available
    if (is.null(names(predictive)) || is.null(names(pooled_oof)) ||
        !identical(names(predictive), names(pooled_oof))) {
      stop("ConditionGRNFit OOF reliability vectors are misaligned.",
           call. = FALSE)
    }
    provenance <- as.data.frame(
      fit$cell_provenance, stringsAsFactors = FALSE
    )
    if (!all(c("cell_id", "condition", "cell_type") %in%
             colnames(provenance)) ||
        anyNA(provenance[, c("condition", "cell_type"), drop = FALSE]) ||
        !all(as.character(provenance$cell_type) == fit_cell_type)) {
      stop("ConditionGRNFit cell provenance is invalid.", call. = FALSE)
    }

    for (condition in fit$condition_levels) {
      tab <- edge
      tab$tf <- toupper(trimws(as.character(tab$tf)))
      tab$target <- toupper(trimws(as.character(tab$target)))
      tab$region <- trimws(as.character(tab$region))
      tab$condition_estimate <- as.numeric(beta[, condition])
      tab$condition_effect <- tab$condition_estimate
      tab$effect_definition <- "absolute_condition_coefficient"
      tab$estimate <- tab$condition_estimate
      tab$corr <- NA_real_
      tab$eligible_in_condition <- as.logical(estimability[, condition])
      tab$structural_candidate <- as.logical(structural[, condition])
      tab$screened_in_condition <- as.logical(screening[, condition])
      tab$selected_support <- as.logical(support[, condition])
      tab$active_in_condition <- as.logical(active[, condition])
      tab$predictor_center <- transform$center
      tab$predictor_scale <- transform$scale
      response_index <- match(tab$target, response$target)
      rsq_index <- match(tab$target, toupper(rownames(target_rsq)))
      if (anyNA(response_index) || anyNA(rsq_index)) {
        stop("ConditionGRNFit target diagnostics are misaligned.",
             call. = FALSE)
      }
      tab$response_center <- response$center[response_index]
      tab$response_scale <- response$scale[response_index]
      tab$rsq_train <- target_rsq[rsq_index, condition]
      tab$rsq_oof_pooled <- pooled_oof[rsq_index]
      tab$predictive_oof_available <- as.logical(predictive[rsq_index])
      tab$oof_reliability_available <-
        tab$predictive_oof_available & is.finite(tab$rsq_oof_pooled)
      tab$reliability_status <- ifelse(
        tab$oof_reliability_available,
        "outer_condition_stratified_cell_oof",
        "unavailable_outer_condition_stratified_cell_oof"
      )
      tab$rsq <- tab$rsq_oof_pooled
      tab[[condition_col]] <- condition
      tab[[celltype_col]] <- fit_cell_type
      tab$group_id <- rc_make_stratum_id(
        tab[1L, c(condition_col, celltype_col), drop = FALSE],
        c(condition_col, celltype_col)
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
    summary[[celltype_col]] <- fit_cell_type
    summary$summary_only <- TRUE
    summary$coefficient_contract <- "absolute_condition_effects_only"
    universal_rows[[length(universal_rows) + 1L]] <- summary

    # The legacy writer still reads this field while constructing provenance.
    # It is removed from every returned and persisted artifact below.
    fit$reference_condition <- NA_character_
    fits_internal[[fit_index]] <- fit
  }

  all_edges <- do.call(rbind, condition_rows)
  rownames(all_edges) <- NULL
  reliable <- all_edges$oof_reliability_available &
    is.finite(all_edges$rsq_oof_pooled) &
    all_edges$rsq_oof_pooled >= min_model_rsq
  condition_active <- all_edges[
    is.finite(all_edges$condition_estimate) &
      abs(all_edges$condition_estimate) >= active_tol &
      all_edges$eligible_in_condition %in% TRUE & reliable,
    , drop = FALSE
  ]
  effect_all <- all_edges
  effect_all$estimate <- effect_all$condition_effect
  effect_active <- effect_all[
    is.finite(effect_all$condition_effect) &
      abs(effect_all$condition_effect) >= active_tol &
      effect_all$eligible_in_condition %in% TRUE & reliable,
    , drop = FALSE
  ]
  params <- methods::slot(methods::slot(grn_object, "grn"), "params")
  list(
    network_index = .rc_strip_reference_contract(
      params$condition_network_index %||% data.frame()
    ),
    fit_diagnostics = params$condition_fit_diagnostics %||% data.frame(),
    fit_contracts = fits_internal,
    universal = do.call(rbind, universal_rows),
    condition_all = all_edges,
    condition_active = condition_active,
    condition_effect_all = effect_all,
    condition_effect_active = effect_active,
    active_tol = active_tol,
    coefficient_contract = "absolute_condition_effects_only"
  )
}

.rc_fit_condition_grns_by_cell_type_reference <-
  .rc_fit_condition_grns_by_cell_type
.rc_fit_condition_grns_by_cell_type <- function(...) {
  call <- match.call(expand.dots = TRUE)
  outdir <- eval(call$outdir, parent.frame())
  answer <- .rc_fit_condition_grns_by_cell_type_reference(...)
  answer$condition_grn_fits <- lapply(
    answer$condition_grn_fits, .rc_strip_reference_contract
  )
  absolute_tables <- c(
    "tf_peak_gene_universal", "tf_peak_gene_condition_all",
    "tf_peak_gene_condition", "tf_peak_gene_condition_effect_all",
    "tf_peak_gene_condition_effect"
  )
  for (field in absolute_tables) {
    if (is.data.frame(answer[[field]])) {
      answer[[field]] <- .rc_strip_reference_contract(answer[[field]])
    }
  }
  answer$normalization_policy$reference_condition <- NULL
  answer$normalization_policy$condition_effect <-
    "absolute condition coefficient on the shared equal-condition coordinate"
  answer$normalization_policy$coefficient_contract <-
    "absolute_condition_effects_only"
  answer$reference_contrast <- NULL

  saveRDS(
    answer$condition_grn_fits,
    file.path(outdir, "pando_condition_grn_fits.rds")
  )
  predictor_transforms <- do.call(rbind, lapply(
    answer$condition_grn_fits,
    function(fit) merge(
      fit$edge_table,
      fit$predictor_transform,
      by = "edge_id", all.x = TRUE, sort = FALSE
    )
  ))
  .rc_mm_write_tsv_gz(
    predictor_transforms,
    file.path(outdir, "pando_edge_predictor_transforms.tsv.gz")
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
