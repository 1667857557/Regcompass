# Integration of the Pando E-star/JSE conditional-GRN contract.

.RC_PANDO_CONDITION_GRN_FIT_SCHEMA <-
  "pando_condition_grn_common_dictionary_v1"
.RC_PANDO_CONDITION_GRN_MODEL_SCHEMA <-
  "pando_condition_grn_Estar_jointse_v1"
.RC_PANDO_CONDITION_GRN_ENGINE <-
  "condition_union_Estar_z025_jointse"
.RC_PANDO_CONDITION_INFERENCE_SCHEMA <-
  "scheme_e_fusion_component_joint_refit_v1"
.RC_PANDO_CONDITION_PROJECTION_POLICY <-
  "any_condition_padj_exact_edge_union"
.RC_PANDO_CONDITION_DICTIONARY_POLICY <-
  "global_and_condition_union_pando_correlation_supported_frozen_dictionary"
.RC_PANDO_CONDITION_PENALTY_FAMILY <-
  "information_scaled_sparse_deviation"
.RC_PANDO_CONDITION_BH_SCOPE <- "condition_target_BH"
.RC_PANDO_CONDITION_SCHEME_E_Z <- 0.25

.rc_pando_execution_summary <- function(diagnostics = NULL, fits = NULL) {
  engines <- if (is.list(fits) && length(fits)) {
    unique(vapply(
      fits,
      function(fit) as.character(fit$fit_engine)[[1L]],
      character(1)
    ))
  } else {
    .RC_PANDO_CONDITION_GRN_ENGINE
  }
  list(
    fit_engine = paste(engines, collapse = ";"),
    targets_total = if (is.data.frame(diagnostics)) {
      length(unique(as.character(diagnostics$target)))
    } else 0L,
    targets_failed = if (is.data.frame(diagnostics) &&
      "fit_status" %in% colnames(diagnostics)) {
      sum(trimws(as.character(diagnostics$fit_status)) != "ok", na.rm = TRUE)
    } else 0L
  )
}

.rc_split_support_conditions <- function(value) {
  value <- as.character(value)
  lapply(value, function(one) {
    if (is.na(one) || !nzchar(one)) character() else
      sort(unique(strsplit(one, ";", fixed = TRUE)[[1L]]))
  })
}

.rc_condition_bh_expected <- function(coefficient) {
  expected <- rep(NA_real_, nrow(coefficient))
  key <- paste(
    as.character(coefficient$condition),
    toupper(trimws(as.character(coefficient$target))),
    sep = "\001"
  )
  for (family in unique(key)) {
    index <- which(key == family)
    valid <- index[
      coefficient$inference_estimable[index] %in% TRUE &
      is.finite(suppressWarnings(as.numeric(coefficient$pval[index])))
    ]
    if (length(valid)) {
      expected[valid] <- stats::p.adjust(
        as.numeric(coefficient$pval[valid]), method = "BH"
      )
    }
  }
  expected
}

.rc_require_pando_condition_grn_fit <- function(fit) {
  if (!inherits(fit, "ConditionGRNFit") ||
      !identical(fit$schema_version, .RC_PANDO_CONDITION_GRN_FIT_SCHEMA)) {
    stop(
      "RegCompass requires Pando common-dictionary condition GRN schema `",
      .RC_PANDO_CONDITION_GRN_FIT_SCHEMA, "`.", call. = FALSE
    )
  }
  required <- c(
    "model_schema", "fit_engine", "inference_schema", "cell_type",
    "condition_levels", "reference_condition", "condition_cell_ids",
    "edge_dictionary", "dictionary_support_table",
    "dictionary_support_summary", "coefficients", "contrasts", "fit",
    "network_names", "padj_threshold", "adjust_method", "scale",
    "interaction", "projection_effect_column", "projection_policy",
    "rna_layer", "peak_layer", "peak_value_type",
    "preprocessing_fingerprint", "condition_col", "cell_type_col",
    "coefficient_scale", "target_genes", "fit_dictionary_policy",
    "candidate_edge_count", "fit_dictionary_edge_count",
    "candidate_tf_cor", "candidate_peak_cor", "deviation_penalty",
    "target_solver", "target_scaling", "target_contrast_tree",
    "rsq_definition"
  )
  if (!all(required %in% names(fit))) {
    stop("Pando E-star/JSE condition fit contract is incomplete.",
         call. = FALSE)
  }
  threshold <- suppressWarnings(as.numeric(fit$padj_threshold))
  penalty <- fit$deviation_penalty
  penalty_ok <- is.list(penalty) &&
    identical(
      as.character(penalty$family),
      .RC_PANDO_CONDITION_PENALTY_FAMILY
    ) &&
    isTRUE(all.equal(
      as.numeric(penalty$z), .RC_PANDO_CONDITION_SCHEME_E_Z,
      tolerance = 1e-15
    ))
  if (!identical(fit$model_schema, .RC_PANDO_CONDITION_GRN_MODEL_SCHEMA) ||
      !identical(fit$fit_engine, .RC_PANDO_CONDITION_GRN_ENGINE) ||
      !identical(fit$inference_schema,
                 .RC_PANDO_CONDITION_INFERENCE_SCHEMA) ||
      !identical(fit$projection_policy,
                 .RC_PANDO_CONDITION_PROJECTION_POLICY) ||
      !identical(fit$fit_dictionary_policy,
                 .RC_PANDO_CONDITION_DICTIONARY_POLICY) ||
      !identical(fit$projection_effect_column, "penalty_effect") ||
      !identical(fit$scale, FALSE) ||
      !identical(fit$interaction, ":") ||
      !identical(toupper(as.character(fit$adjust_method)), "BH") ||
      length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1 || !penalty_ok) {
    stop("Pando E-star z=0.25/JSE model labels are inconsistent.",
         call. = FALSE)
  }

  levels <- as.character(fit$condition_levels)
  reference <- as.character(fit$reference_condition)
  if (length(levels) < 2L || anyNA(levels) || any(!nzchar(levels)) ||
      anyDuplicated(levels) || length(reference) != 1L ||
      !reference %in% levels) {
    stop("Pando condition levels/reference are invalid.", call. = FALSE)
  }
  cells <- fit$condition_cell_ids[levels]
  if (!is.list(cells) || any(lengths(cells) < 3L)) {
    stop("Every Pando condition must retain at least three paired cells.",
         call. = FALSE)
  }
  all_cells <- as.character(unlist(cells, use.names = FALSE))
  if (anyNA(all_cells) || any(!nzchar(all_cells)) || anyDuplicated(all_cells)) {
    stop("Pando condition cell IDs must be complete and disjoint.",
         call. = FALSE)
  }

  edge <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
  required_edge <- c(
    "edge_id", "target", "tf", "region", "atac_feature_id",
    "candidate_index", "source_global", "source_conditions", "n_sources"
  )
  if (!nrow(edge) || !all(required_edge %in% colnames(edge)) ||
      anyDuplicated(as.character(edge$edge_id))) {
    stop("Pando exact-edge dictionary is incomplete.", call. = FALSE)
  }
  expected_edge_id <- paste(edge$target, edge$tf, edge$region, sep = "||")
  if (!identical(as.character(edge$edge_id), expected_edge_id)) {
    stop("Pando edge_id does not match exact target-TF-region coordinates.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(as.integer(fit$candidate_edge_count), nrow(edge))) ||
      !isTRUE(all.equal(as.integer(fit$fit_dictionary_edge_count), nrow(edge)))) {
    stop("Pando fit changed the frozen exact-edge dictionary size.",
         call. = FALSE)
  }

  support <- as.data.frame(
    fit$dictionary_support_table, stringsAsFactors = FALSE
  )
  required_support <- c(
    "edge_id", "source_type", "condition",
    "peak_target_cor", "tf_target_cor"
  )
  if (!nrow(support) ||
      !all(required_support %in% colnames(support)) ||
      any(!as.character(support$edge_id) %in% as.character(edge$edge_id))) {
    stop("Pando candidate-support provenance is incomplete.", call. = FALSE)
  }
  support_key <- paste(
    as.character(support$edge_id),
    as.character(support$source_type),
    ifelse(
      as.character(support$source_type) == "global",
      "__global__", as.character(support$condition)
    ),
    sep = "\001"
  )
  if (anyDuplicated(support_key)) {
    stop("Pando candidate-support provenance contains duplicate sources.",
         call. = FALSE)
  }

  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required_coefficient <- c(
    "edge_id", "target", "tf", "region", "condition", "z",
    "estimate", "penalty_effect", "beta_shared", "condition_deviation",
    "contrast_identifiable", "shared_by_boundary", "boundary_condition",
    "fused_by_penalty", "fusion_component_id", "shared_edge",
    "raw_information_condition", "profile_information",
    "inference_schema", "inference_component_id",
    "inference_hypothesis_id", "inference_estimate", "inference_se",
    "inference_statistic", "inference_estimable", "pval", "padj",
    "bh_scope", "bh_family_size", "condition_significant",
    "statistically_supported", "significant", "pando_estimation_active",
    "active", "edge_union_supported", "supporting_conditions",
    "n_supporting_conditions", "all_conditions_fit_valid",
    "active_in_regcompass", "fit_status", "penalty_family",
    "penalty_value", "solver_status", "objective", "kkt_residual",
    "iterations"
  )
  if (!nrow(coefficient) ||
      !all(required_coefficient %in% colnames(coefficient))) {
    stop("Pando E-star/JSE coefficient rows are incomplete.", call. = FALSE)
  }
  dictionary_ids <- sort(as.character(edge$edge_id))
  for (condition in levels) {
    one <- coefficient[
      as.character(coefficient$condition) == condition, , drop = FALSE
    ]
    if (nrow(one) != nrow(edge) || anyDuplicated(one$edge_id) ||
        !identical(sort(as.character(one$edge_id)), dictionary_ids)) {
      stop("Every condition must contain every frozen exact edge once.",
           call. = FALSE)
    }
  }

  if (any(as.character(coefficient$inference_schema) !=
          .RC_PANDO_CONDITION_INFERENCE_SCHEMA) ||
      any(as.character(coefficient$bh_scope) !=
          .RC_PANDO_CONDITION_BH_SCOPE) ||
      any(as.character(coefficient$penalty_family) !=
          .RC_PANDO_CONDITION_PENALTY_FAMILY) ||
      any(abs(as.numeric(coefficient$penalty_value) -
              .RC_PANDO_CONDITION_SCHEME_E_Z) > 1e-15) ||
      any(abs(as.numeric(coefficient$z) -
              .RC_PANDO_CONDITION_SCHEME_E_Z) > 1e-15) ||
      any(as.character(coefficient$solver_status) != "ok") ||
      any(!is.finite(as.numeric(coefficient$kkt_residual)))) {
    stop("Pando coefficient rows are not converged E-star z=0.25/JSE results.",
         call. = FALSE)
  }

  estimate <- suppressWarnings(as.numeric(coefficient$estimate))
  effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  comparable <- is.finite(estimate) & is.finite(effect)
  if (any(is.finite(estimate) != is.finite(effect)) ||
      any(abs(estimate[comparable] - effect[comparable]) > 1e-12)) {
    stop("Pando penalty_effect does not equal production beta_E.",
         call. = FALSE)
  }
  if (!identical(
      as.logical(coefficient$active), is.finite(estimate)
  ) || !identical(
      as.logical(coefficient$pando_estimation_active), is.finite(estimate)
  )) {
    stop("Pando active flags must retain every finite production coefficient.",
         call. = FALSE)
  }

  expected_padj <- .rc_condition_bh_expected(coefficient)
  observed_padj <- suppressWarnings(as.numeric(coefficient$padj))
  comparable <- is.finite(expected_padj) & is.finite(observed_padj)
  if (any(is.finite(expected_padj) != is.finite(observed_padj)) ||
      any(abs(expected_padj[comparable] - observed_padj[comparable]) > 1e-10)) {
    stop("Pando padj is not BH within condition x target.", call. = FALSE)
  }
  expected_significant <-
    coefficient$inference_estimable %in% TRUE &
    is.finite(observed_padj) & observed_padj < threshold
  if (!identical(
      as.logical(coefficient$condition_significant), expected_significant
  ) || !identical(
      as.logical(coefficient$statistically_supported), expected_significant
  ) || !identical(
      as.logical(coefficient$significant), expected_significant
  )) {
    stop("Pando condition-significance flags do not match BH padj < threshold.",
         call. = FALSE)
  }

  for (edge_id in dictionary_ids) {
    index <- which(as.character(coefficient$edge_id) == edge_id)
    valid <- length(index) == length(levels) &&
      all(as.character(coefficient$fit_status[index]) == "ok") &&
      all(is.finite(effect[index]))
    supporting <- levels[levels %in% as.character(
      coefficient$condition[index][expected_significant[index]]
    )]
    union <- valid && length(supporting) > 0L
    if (any(as.logical(coefficient$all_conditions_fit_valid[index]) != valid) ||
        any(as.logical(coefficient$edge_union_supported[index]) != union) ||
        any(as.logical(coefficient$active_in_regcompass[index]) != union) ||
        any(as.integer(coefficient$n_supporting_conditions[index]) !=
            length(supporting)) ||
        any(as.character(coefficient$supporting_conditions[index]) !=
            paste(supporting, collapse = ";"))) {
      stop("Pando any-condition exact-edge union metadata are inconsistent.",
           call. = FALSE)
    }
  }

  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  required_fit <- c(
    "target", "condition", "rsq", "rsq_definition", "fit_status",
    "sigma2_common", "inference_sigma2", "inference_residual_df",
    "reference_condition", "deviation_z", "penalty_family",
    "solver_status", "kkt_residual", "iterations",
    "predictor_scale_reference", "inference_schema",
    "orthogonality_error", "dr_error"
  )
  if (!nrow(fit_table) || !all(required_fit %in% colnames(fit_table)) ||
      any(as.character(fit_table$fit_status) != "ok") ||
      any(as.character(fit_table$penalty_family) !=
          .RC_PANDO_CONDITION_PENALTY_FAMILY) ||
      any(as.character(fit_table$inference_schema) !=
          .RC_PANDO_CONDITION_INFERENCE_SCHEMA) ||
      any(as.character(fit_table$predictor_scale_reference) !=
          "equal_condition_within_condition_rms") ||
      any(abs(as.numeric(fit_table$deviation_z) -
              .RC_PANDO_CONDITION_SCHEME_E_Z) > 1e-15) ||
      any(as.numeric(fit_table$orthogonality_error) > 1e-6) ||
      any(as.numeric(fit_table$dr_error) > 1e-6)) {
    stop("Pando target-level E-star/JSE diagnostics are incomplete.",
         call. = FALSE)
  }
  fit_key <- paste(
    toupper(as.character(fit_table$target)),
    as.character(fit_table$condition), sep = "\001"
  )
  if (anyDuplicated(fit_key)) {
    stop("Pando target-condition diagnostics are duplicated.",
         call. = FALSE)
  }

  contrast <- as.data.frame(fit$contrasts, stringsAsFactors = FALSE)
  required_contrast <- c(
    "edge_id", "target", "condition_a", "condition_b",
    "estimate_a", "estimate_b", "contrast_estimate",
    "inference_contrast_estimate", "contrast_se", "contrast_pval",
    "contrast_padj", "contrast_estimable", "contrast_identifiable",
    "contrast_status", "contrast_significant", "shared_by_boundary",
    "fused_by_penalty", "shared_edge", "penalty_family",
    "penalty_value", "solver_status", "kkt_residual", "iterations"
  )
  expected_n <- choose(length(levels), 2L) * nrow(edge)
  if (!is.data.frame(contrast) || nrow(contrast) != expected_n ||
      !all(required_contrast %in% colnames(contrast))) {
    stop("Pando pairwise contrast table is incomplete.", call. = FALSE)
  }
  for (i in seq_len(nrow(contrast))) {
    edge_id <- as.character(contrast$edge_id[[i]])
    a <- as.character(contrast$condition_a[[i]])
    b <- as.character(contrast$condition_b[[i]])
    ia <- which(
      as.character(coefficient$edge_id) == edge_id &
      as.character(coefficient$condition) == a
    )
    ib <- which(
      as.character(coefficient$edge_id) == edge_id &
      as.character(coefficient$condition) == b
    )
    if (length(ia) != 1L || length(ib) != 1L) {
      stop("Pando contrast cannot be aligned to production coefficients.",
           call. = FALSE)
    }
    expected_delta <- estimate[[ib]] - estimate[[ia]]
    observed_delta <- as.numeric(contrast$contrast_estimate[[i]])
    if (!is.finite(observed_delta) ||
        abs(observed_delta - expected_delta) > 1e-9) {
      stop("Pando pairwise contrasts do not close on joint production beta_E.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col) {
  fits <- Pando::condition_grn_fit(grn_object)
  if (inherits(fits, "ConditionGRNFit")) fits <- list(fits)
  if (!is.list(fits) || !length(fits)) {
    stop("Pando did not return common-dictionary condition fits.",
         call. = FALSE)
  }
  invisible(lapply(fits, .rc_require_pando_condition_grn_fit))

  data_object <- methods::slot(grn_object, "data")
  metadata <- methods::slot(data_object, "meta.data")
  if (!is.data.frame(metadata) ||
      !all(c(condition_col, celltype_col) %in% colnames(metadata)) ||
      is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
    stop("Pando object metadata cannot validate condition-fit cell mappings.",
         call. = FALSE)
  }

  rows <- universal <- diagnostics <- contrasts <- gated_fits <- list()
  for (fit_index in seq_along(fits)) {
    fit <- .rc_apply_condition_penalty_gate(fits[[fit_index]])
    gated_fits[[fit_index]] <- fit
    for (condition in fit$condition_levels) {
      cells <- as.character(fit$condition_cell_ids[[condition]])
      if (any(!cells %in% rownames(metadata)) ||
          any(as.character(metadata[cells, condition_col]) != condition) ||
          any(as.character(metadata[cells, celltype_col]) !=
              as.character(fit$cell_type))) {
        stop("Pando fitted cell assignments disagree with object metadata.",
             call. = FALSE)
      }
    }

    threshold <- as.numeric(fit$padj_threshold)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    coefficient$tf <- toupper(as.character(coefficient$tf))
    coefficient$target <- toupper(as.character(coefficient$target))
    coefficient$region <- as.character(coefficient$region)
    coefficient$condition <- as.character(coefficient$condition)
    coefficient[[condition_col]] <- coefficient$condition
    coefficient[[celltype_col]] <- as.character(fit$cell_type)
    coefficient$group_id <- rc_make_stratum_id(
      coefficient[, c(condition_col, celltype_col), drop = FALSE],
      c(condition_col, celltype_col)
    )
    coefficient$condition_estimate <- as.numeric(coefficient$estimate)
    coefficient$condition_effect <- as.numeric(coefficient$penalty_effect)
    coefficient$effect_definition <-
      "E_star_z025_continuous_condition_coefficient_raw_tf_atac_units"
    coefficient$coefficient_contract <-
      "fixed_exact_edge_dictionary_Estar_z025_jointse"
    coefficient$fit_engine <- fit$fit_engine
    coefficient$coefficient_scale <- fit$coefficient_scale
    coefficient$padj_threshold <- threshold
    coefficient$eligible_in_condition <- coefficient$penalty_eligible %in% TRUE
    coefficient$biological_differential_eligible <-
      coefficient$contrast_identifiable %in% TRUE

    fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
    fit_table$target <- toupper(as.character(fit_table$target))
    fit_table$condition <- as.character(fit_table$condition)
    fit_key <- paste(fit_table$target, fit_table$condition, sep = "\001")
    coef_key <- paste(coefficient$target, coefficient$condition, sep = "\001")
    index <- match(coef_key, fit_key)
    if (anyNA(index)) {
      stop("Condition coefficients and target diagnostics are misaligned.",
           call. = FALSE)
    }
    coefficient$rsq <- as.numeric(fit_table$rsq[index])
    coefficient$fit_status <- as.character(fit_table$fit_status[index])
    coefficient$reliable_model <- coefficient$fit_status == "ok"
    rows[[length(rows) + 1L]] <- coefficient
    diagnostics[[length(diagnostics) + 1L]] <- fit_table

    edge <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
    one_condition <- fit$condition_levels[[1L]]
    first_rows <- coefficient[coefficient$condition == one_condition, , drop = FALSE]
    shared_by_id <- stats::setNames(
      as.numeric(first_rows$beta_shared), as.character(first_rows$edge_id)
    )
    summary <- edge
    summary$estimate <- unname(shared_by_id[as.character(summary$edge_id)])
    summary$corr <- NA_real_
    summary[[celltype_col]] <- as.character(fit$cell_type)
    summary$summary_only <- TRUE
    summary$coefficient_contract <-
      "E_star_unpenalized_shared_MLE_summary_only"
    universal[[length(universal) + 1L]] <- summary

    contrast <- as.data.frame(fit$contrasts, stringsAsFactors = FALSE)
    contrast$tf <- toupper(as.character(contrast$tf))
    contrast$target <- toupper(as.character(contrast$target))
    contrast[[celltype_col]] <- as.character(fit$cell_type)
    contrast$fit_engine <- fit$fit_engine
    contrast$padj_threshold <- threshold
    contrast$biological_differential_eligible <-
      contrast$contrast_identifiable %in% TRUE
    contrast$differential_claim_eligible <-
      contrast$contrast_identifiable %in% TRUE &
      contrast$contrast_status != "not_estimable"
    contrasts[[length(contrasts) + 1L]] <- contrast
  }
  names(gated_fits) <- names(fits)

  all_edges <- do.call(rbind, rows)
  rownames(all_edges) <- NULL
  active <- all_edges[all_edges$penalty_eligible %in% TRUE, , drop = FALSE]
  all_contrasts <- do.call(rbind, contrasts)
  rownames(all_contrasts) <- NULL
  params <- methods::slot(methods::slot(grn_object, "grn"), "params")
  list(
    network_index = params$condition_network_index %||% data.frame(),
    fit_diagnostics = do.call(rbind, diagnostics),
    fit_contracts = gated_fits,
    universal = do.call(rbind, universal),
    condition_all = all_edges,
    condition_active = active,
    condition_effect_all = all_edges,
    condition_effect_active = active,
    condition_contrasts = all_contrasts,
    active_tol = 0,
    penalty_filter = paste(
      "all conditions valid and any condition-target BH padj < threshold;",
      "retain every condition production beta_E for the admitted exact edge"
    ),
    coefficient_contract =
      "fixed_exact_edge_dictionary_Estar_z025_jointse",
    pando_fit_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    pando_memory_contract = "pando_native_condition_Estar_z025_jointse"
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
      tf_cor = 0.05, peak_cor = 0.05, adjust_method = "BH",
      padj_threshold = 0.05, rank_action = "mark",
      min_residual_df = 1L, reference_condition = NULL
    ),
    save_pando_objects = TRUE, BPPARAM = NULL,
    progress_monitor = NULL,
    species = c("auto", "human", "mouse")) {
  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  .rc_validate_condition_celltype_metadata(
    object@meta.data, condition_col, celltype_col,
    require_multiple_conditions = TRUE
  )
  if (!is.list(pando_infer_args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  allowed_infer_args <- c(
    "tf_cor", "peak_cor", "adjust_method", "padj_threshold",
    "rank_action", "min_residual_df", "reference_condition",
    "rna_layer", "peak_layer", "peak_value_type"
  )
  unknown <- setdiff(names(pando_infer_args), allowed_infer_args)
  if (length(unknown)) {
    stop(
      "Unsupported `pando_infer_args`: ", paste(unknown, collapse = ", "),
      ". Conditional E-star uses fixed z=0.25 and does not expose a ",
      "sensitivity/CV control surface.", call. = FALSE
    )
  }
  pando_infer_args <- utils::modifyList(list(
    tf_cor = 0.05, peak_cor = 0.05, adjust_method = "BH",
    padj_threshold = 0.05, rank_action = "mark", min_residual_df = 1L,
    reference_condition = NULL,
    rna_layer = "data", peak_layer = "data",
    peak_value_type = "normalized"
  ), pando_infer_args)
  threshold <- suppressWarnings(as.numeric(pando_infer_args$padj_threshold))
  tf_threshold <- suppressWarnings(as.numeric(pando_infer_args$tf_cor))
  peak_threshold <- suppressWarnings(as.numeric(pando_infer_args$peak_cor))
  if (!identical(toupper(as.character(pando_infer_args$adjust_method)), "BH") ||
      length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1 ||
      length(tf_threshold) != 1L || !is.finite(tf_threshold) ||
      tf_threshold < 0 || tf_threshold > 1 ||
      length(peak_threshold) != 1L || !is.finite(peak_threshold) ||
      peak_threshold < 0 || peak_threshold > 1) {
    stop(
      "Canonical conditional GRNs require BH, padj_threshold in (0,1), ",
      "and tf_cor/peak_cor in [0,1].", call. = FALSE
    )
  }
  pando_infer_args$padj_threshold <- threshold
  pando_infer_args$tf_cor <- tf_threshold
  pando_infer_args$peak_cor <- peak_threshold
  if (!is.null(pando_infer_args$reference_condition)) {
    reference <- as.character(pando_infer_args$reference_condition)
    if (length(reference) != 1L || is.na(reference) ||
        !nzchar(trimws(reference)) || reference != trimws(reference)) {
      stop(
        "`pando_infer_args$reference_condition` must be NULL or one complete ",
        "predefined condition label.", call. = FALSE
      )
    }
    pando_infer_args$reference_condition <- reference
  }

  condition_types <- if (is.null(cell_type)) {
    unique(as.character(object@meta.data[[celltype_col]]))
  } else unique(as.character(cell_type))
  missing_types <- setdiff(
    condition_types, unique(as.character(object@meta.data[[celltype_col]]))
  )
  if (length(missing_types)) {
    stop(
      "Requested condition-GRN cell type(s) were not found: ",
      paste(missing_types, collapse = ", "), call. = FALSE
    )
  }
  plans <- .rc_condition_parallel_plan(
    metadata = object@meta.data,
    condition_types = condition_types,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_cells = min_cells
  )
  condition_parallel <- !identical(BPPARAM, FALSE) && !is.null(BPPARAM)
  worker_limit <- if (condition_parallel) {
    .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  } else 1L

  .rc_step_monitor_event(
    progress_monitor, "condition_design",
    "resolved Pando E-star z=0.25/JSE conditional-GRN plan",
    current = 5L,
    context = list(
      cell_types = length(plans), tf_cor = tf_threshold,
      peak_cor = peak_threshold, padj_threshold = threshold,
      reference_condition = pando_infer_args$reference_condition %||%
        "<first-retained>",
      bh_scope = .RC_PANDO_CONDITION_BH_SCOPE,
      scheme_e_z = .RC_PANDO_CONDITION_SCHEME_E_Z,
      workers = worker_limit,
      estimator = .RC_PANDO_CONDITION_GRN_ENGINE
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
    stop("No overlap between RNA genes and GEM metabolic genes.",
         call. = FALSE)
  }

  prepare_tasks <- lapply(names(plans), function(type) {
    list(
      cell_type = type,
      object = subset(object, cells = plans[[type]]$global_cells)
    )
  })
  outer_parallel <- condition_parallel && length(prepare_tasks) > 1L
  prepared <- rc_parallel_lapply(
    prepare_tasks,
    .rc_condition_prepare_celltype_task,
    BPPARAM = if (outer_parallel) BPPARAM else FALSE,
    atac_assay = atac_assay,
    rna_assay = rna_assay,
    pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args,
    pfm = pfm,
    genome = genome
  )
  prepared_types <- vapply(prepared, `[[`, character(1), "cell_type")
  if (anyDuplicated(prepared_types) ||
      !setequal(prepared_types, names(plans))) {
    stop("Condition-GRN Pando initialization returned an invalid cell-type set.",
         call. = FALSE)
  }
  names(prepared) <- prepared_types
  prepared <- prepared[names(plans)]
  invisible(gc(verbose = FALSE, full = TRUE))

  fit_tasks <- lapply(names(plans), function(type) {
    list(cell_type = type, grn = prepared[[type]]$grn)
  })
  inner_parallel <- condition_parallel && length(fit_tasks) == 1L
  fit_results <- rc_parallel_lapply(
    fit_tasks,
    .rc_condition_ridge_fit_task,
    BPPARAM = if (outer_parallel) BPPARAM else FALSE,
    target_genes = target_genes,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_cells = min_cells,
    pando_infer_args = pando_infer_args,
    inner_parallel = inner_parallel,
    PANDO_BPPARAM = if (inner_parallel) BPPARAM else NULL
  )
  fit_types <- vapply(fit_results, `[[`, character(1), "cell_type")
  if (anyDuplicated(fit_types) || !setequal(fit_types, names(plans))) {
    stop("Pando E-star/JSE fits returned an invalid cell-type set.",
         call. = FALSE)
  }
  names(fit_results) <- fit_types
  fit_results <- fit_results[names(plans)]
  invisible(gc(verbose = FALSE, full = TRUE))

  parallel_plan <- list(
    scope = "cell_type_condition_Estar_JSE",
    cell_type_prepare_tasks = length(prepare_tasks),
    condition_Estar_fit_tasks = length(fit_tasks),
    workers = worker_limit,
    outer_celltype_parallel = outer_parallel,
    inner_target_parallel = inner_parallel,
    nested_parallel = FALSE,
    stage_barrier = paste(
      "global + condition candidate discovery -> exact union ->",
      "equal-condition RMS scaling -> pooled residual information ->",
      "Q-orthogonal E-star z=0.25 -> fusion-component joint refit ->",
      "condition-target BH -> any-condition exact-edge RegCompass union"
    )
  )

  results <- list()
  meta <- object@meta.data
  for (type in names(plans)) {
    grn <- fit_results[[type]]$grn
    fit_contract <- fit_results[[type]]$fit
    invisible(.rc_require_pando_condition_grn_fit(fit_contract))
    extracted <- .rc_extract_condition_grn_contract(
      grn, condition_col, celltype_col
    )
    execution_summary <- .rc_pando_execution_summary(
      extracted$fit_diagnostics, extracted$fit_contracts
    )
    execution_summary$parallel_plan <- parallel_plan

    status_rows <- list()
    for (condition in fit_contract$condition_levels) {
      cells <- fit_contract$condition_cell_ids[[condition]]
      all_rows <- extracted$condition_all[
        extracted$condition_all[[condition_col]] == condition &
          extracted$condition_all[[celltype_col]] == type,
        , drop = FALSE
      ]
      active_rows <- extracted$condition_active[
        extracted$condition_active[[condition_col]] == condition &
          extracted$condition_active[[celltype_col]] == type,
        , drop = FALSE
      ]
      key_frame <- data.frame(
        condition_value = condition,
        celltype_value = type,
        stringsAsFactors = FALSE
      )
      names(key_frame) <- c(condition_col, celltype_col)
      one_status <- data.frame(
        group_id = rc_make_stratum_id(
          key_frame, c(condition_col, celltype_col)
        ),
        condition_value = condition,
        celltype_value = type,
        reference_condition = fit_contract$reference_condition,
        n_cells = length(cells), status = "ok",
        n_target_genes = length(unique(all_rows$target)),
        n_edges = nrow(all_rows),
        n_active_edges = nrow(active_rows),
        n_condition_significant_edges =
          sum(all_rows$condition_significant %in% TRUE),
        n_regcompass_union_edges =
          sum(all_rows$active_in_regcompass %in% TRUE),
        n_contrast_identifiable_edges =
          sum(all_rows$contrast_identifiable %in% TRUE),
        n_shared_by_boundary =
          sum(all_rows$shared_by_boundary %in% TRUE),
        n_fused_by_penalty =
          sum(all_rows$fused_by_penalty %in% TRUE),
        padj_threshold = threshold,
        scheme_e_z = .RC_PANDO_CONDITION_SCHEME_E_Z,
        stringsAsFactors = FALSE
      )
      names(one_status)[names(one_status) == "condition_value"] <- condition_col
      names(one_status)[names(one_status) == "celltype_value"] <- celltype_col
      status_rows[[length(status_rows) + 1L]] <- one_status
    }
    status <- do.call(rbind, status_rows)
    selected_cells <- unique(unlist(
      fit_contract$condition_cell_ids, use.names = FALSE
    ))

    results[[type]] <- list(
      schema_version = "regcompass_condition_grn_common_dictionary",
      analysis_mode = "condition_grn",
      condition_coefficients_calculated = TRUE,
      pando_fit_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
      pando_model_schema = .RC_PANDO_CONDITION_GRN_MODEL_SCHEMA,
      pando_installed_version = as.character(utils::packageVersion("Pando")),
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
      pando_execution_summary = execution_summary,
      condition_grn_fits = extracted$fit_contracts,
      tf_peak_gene_universal = extracted$universal,
      tf_peak_gene_condition_all = extracted$condition_all,
      tf_peak_gene_condition = extracted$condition_active,
      tf_peak_gene_condition_effect_all = extracted$condition_all,
      tf_peak_gene_condition_effect = extracted$condition_active,
      tf_peak_gene_condition_contrasts = extracted$condition_contrasts,
      normalization_policy = list(
        rna = "global single-cell normalized RNA",
        atac = "cell-type-shared TF-IDF across conditions",
        grn_fit = paste(
          "Pando exact-edge common dictionary, Q-orthogonal E-star z=0.25,",
          "and fusion-component joint-refit inference"
        ),
        condition_effect =
          "continuous production beta_E retained in every condition",
        condition_significance =
          "BH within condition x target; strict padj < threshold",
        regcompass_handoff = paste(
          "all conditions fit-valid plus any condition significant;",
          "same admitted exact-edge set in every condition"
        ),
        target_rsq = "Pando full-data R2 is diagnostic only",
        parallel_contract = parallel_plan,
        penalty_regulatory_evidence = paste(
          "condition-specific production penalty_effect times the canonical",
          "RegCompass paired TF-ATAC metacell exposure"
        )
      ),
      group_cols = c(condition_col, celltype_col)
    )
  }

  answer <- .rc_merge_condition_job_results(results)
  answer$pando_execution_summary$parallel_plan <- parallel_plan
  answer$normalization_policy$parallel_contract <- parallel_plan
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(
    answer$condition_fit_status,
    file.path(outdir, "pando_group_status.tsv.gz")
  )
  .rc_write_tsv_gz(
    answer$tf_peak_gene_condition_all,
    file.path(outdir, "pando_tf_peak_gene_condition_all.tsv.gz")
  )
  .rc_write_tsv_gz(
    answer$tf_peak_gene_condition,
    file.path(outdir, "pando_tf_peak_gene_condition_active.tsv.gz")
  )
  .rc_write_tsv_gz(
    answer$tf_peak_gene_condition_contrasts,
    file.path(outdir, "pando_tf_peak_gene_condition_contrasts.tsv.gz")
  )
  .rc_write_tsv_gz(
    answer$tf_peak_gene_universal,
    file.path(outdir, "pando_tf_peak_gene_universal.tsv.gz")
  )
  saveRDS(
    answer$condition_grn_fits,
    file.path(outdir, "pando_condition_grn_fits.rds")
  )
  if (isTRUE(save_pando_objects)) {
    dir.create(
      file.path(outdir, "pando_objects"),
      recursive = TRUE, showWarnings = FALSE
    )
    for (type in names(answer$pando_grn_data_by_cell_type)) {
      saveRDS(
        answer$pando_grn_data_by_cell_type[[type]],
        file.path(
          outdir, "pando_objects",
          paste0("condition_grn_", .rc_safe_path_component(type), ".rds")
        )
      )
    }
  }
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
