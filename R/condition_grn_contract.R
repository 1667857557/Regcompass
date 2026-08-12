# Integration of the Pando global-plus-condition common-dictionary GRN contract.

.RC_PANDO_CONDITION_GRN_FIT_SCHEMA <-
  "pando_condition_grn_common_dictionary_v1"
.RC_PANDO_CONDITION_GRN_MODEL_SCHEMA <-
  "pando_condition_grn_multitask_ridge_v3"
.RC_PANDO_CONDITION_GRN_ENGINE <-
  "condition_union_single_no_fusion_common_lambda_ridge"
.RC_PANDO_CONDITION_PROJECTION_POLICY <-
  "active_global_or_local_pando_support_and_condition_bh_ridge_effects"
.RC_PANDO_CONDITION_DICTIONARY_POLICY <-
  "global_and_condition_union_pando_correlation_supported_frozen_dictionary"

.rc_pando_execution_summary <- function(diagnostics = NULL, fits = NULL) {
  engines <- if (is.list(fits) && length(fits)) {
    unique(vapply(fits, function(fit) as.character(fit$fit_engine)[[1L]],
                  character(1)))
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

.rc_require_pando_condition_grn_fit <- function(fit) {
  if (!inherits(fit, "ConditionGRNFit") ||
      !identical(fit$schema_version, .RC_PANDO_CONDITION_GRN_FIT_SCHEMA)) {
    stop(
      "RegCompass requires Pando common-dictionary condition GRN schema `",
      .RC_PANDO_CONDITION_GRN_FIT_SCHEMA, "`.", call. = FALSE
    )
  }
  required <- c(
    "model_schema", "cell_type", "condition_levels", "condition_cell_ids",
    "edge_dictionary", "dictionary_support_table",
    "dictionary_support_summary", "coefficients", "contrasts", "fit",
    "network_names", "padj_threshold", "adjust_method", "scale",
    "interaction", "projection_effect_column", "projection_policy",
    "rna_layer", "peak_layer", "peak_value_type",
    "preprocessing_fingerprint", "condition_col", "cell_type_col",
    "fit_engine", "coefficient_scale", "target_genes",
    "fit_dictionary_policy", "candidate_edge_count",
    "fit_dictionary_edge_count", "candidate_tf_cor", "candidate_peak_cor"
  )
  threshold <- suppressWarnings(as.numeric(fit$padj_threshold))
  threshold_valid <- length(threshold) == 1L && is.finite(threshold) &&
    threshold > 0 && threshold < 1
  if (!all(required %in% names(fit)) ||
      !identical(fit$model_schema, .RC_PANDO_CONDITION_GRN_MODEL_SCHEMA) ||
      !identical(fit$fit_engine, .RC_PANDO_CONDITION_GRN_ENGINE) ||
      !identical(fit$fit_dictionary_policy,
                 .RC_PANDO_CONDITION_DICTIONARY_POLICY) ||
      !identical(fit$scale, FALSE) ||
      !identical(fit$interaction, ":") ||
      !identical(fit$projection_effect_column, "penalty_effect") ||
      !identical(fit$projection_policy,
                 .RC_PANDO_CONDITION_PROJECTION_POLICY) ||
      !identical(toupper(as.character(fit$adjust_method)), "BH") ||
      !threshold_valid ||
      !isTRUE(attr(
        fit$edge_dictionary,
        "preprocessing_provenance_verified",
        exact = TRUE
      )) ||
      any(!nzchar(c(
        as.character(fit$rna_layer), as.character(fit$peak_layer),
        as.character(fit$peak_value_type),
        as.character(fit$preprocessing_fingerprint)
      )))) {
    stop("Pando no-fusion condition fit contract is incomplete.",
         call. = FALSE)
  }

  scalar_text <- function(value) {
    is.character(value) && length(value) == 1L &&
      !is.na(value) && nzchar(trimws(value))
  }
  if (!scalar_text(fit$cell_type) ||
      !scalar_text(fit$condition_col) ||
      !scalar_text(fit$cell_type_col) ||
      identical(fit$condition_col, fit$cell_type_col) ||
      !scalar_text(fit$coefficient_scale)) {
    stop("Pando condition fit identifiers and model labels are invalid.",
         call. = FALSE)
  }

  levels <- as.character(fit$condition_levels)
  if (length(levels) < 2L || anyNA(levels) || any(!nzchar(levels)) ||
      anyDuplicated(levels)) {
    stop("Pando condition levels must contain at least two unique labels.",
         call. = FALSE)
  }
  cells_by_condition <- fit$condition_cell_ids
  cell_list_names <- names(cells_by_condition)
  if (!is.list(cells_by_condition) || is.null(cell_list_names) ||
      anyNA(cell_list_names) || any(!nzchar(cell_list_names)) ||
      anyDuplicated(cell_list_names) || !all(levels %in% cell_list_names)) {
    stop("Pando condition cell IDs are not uniquely named for every condition.",
         call. = FALSE)
  }
  cells_by_condition <- cells_by_condition[levels]
  if (any(lengths(cells_by_condition) < 3L)) {
    stop("Every Pando fitted condition must contain at least three cells.",
         call. = FALSE)
  }
  fitted_cells <- as.character(unlist(cells_by_condition, use.names = FALSE))
  if (!length(fitted_cells) || anyNA(fitted_cells) ||
      any(!nzchar(fitted_cells)) || anyDuplicated(fitted_cells)) {
    stop("Pando fitted cells must be complete and condition-disjoint.",
         call. = FALSE)
  }

  targets <- unique(as.character(fit$target_genes))
  if (!length(targets) || anyNA(targets) || any(!nzchar(targets)) ||
      anyDuplicated(toupper(targets))) {
    stop("Pando fitted target genes are empty or case-ambiguous.",
         call. = FALSE)
  }

  edge <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
  support <- as.data.frame(fit$dictionary_support_table,
                           stringsAsFactors = FALSE)
  support_summary <- as.data.frame(
    fit$dictionary_support_summary, stringsAsFactors = FALSE
  )
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required_edge <- c(
    "edge_id", "target", "tf", "region", "atac_feature_id",
    "candidate_index", "source_global", "source_conditions", "n_sources"
  )
  required_support <- c(
    "edge_id", "target", "tf", "region", "atac_feature_id",
    "source_type", "condition", "peak_target_cor", "tf_target_cor"
  )
  required_summary <- c(
    "edge_id", "source_global", "n_support_conditions",
    "support_conditions", "n_sources"
  )
  required_coefficient <- c(
    "edge_id", "target", "tf", "region", "condition", "estimate",
    "shared_estimate", "condition_deviation", "std_err", "statistic",
    "pval", "padj", "statistically_supported", "global_support",
    "local_support", "active", "significant", "penalty_effect",
    "estimable", "zero_variance", "aliased"
  )
  if (!nrow(edge) || !nrow(support) || !nrow(support_summary) ||
      !all(required_edge %in% colnames(edge)) ||
      !all(required_support %in% colnames(support)) ||
      !all(required_summary %in% colnames(support_summary)) ||
      !all(required_coefficient %in% colnames(coefficient)) ||
      anyNA(edge$edge_id) || any(!nzchar(as.character(edge$edge_id))) ||
      anyDuplicated(edge$edge_id) ||
      any(!coefficient$condition %in% levels)) {
    stop("Pando global/condition common dictionary is incomplete.",
         call. = FALSE)
  }
  expected_edge_id <- paste(edge$target, edge$tf, edge$region, sep = "||")
  if (!identical(as.character(edge$edge_id), expected_edge_id)) {
    stop("Pando dictionary edge IDs are not exact target-TF-region keys.",
         call. = FALSE)
  }
  candidate_index <- suppressWarnings(as.integer(edge$candidate_index))
  if (anyNA(candidate_index) || any(candidate_index < 1L) ||
      anyDuplicated(candidate_index)) {
    stop("Pando dictionary candidate indices are invalid.", call. = FALSE)
  }
  if (!isTRUE(all.equal(as.integer(fit$candidate_edge_count), nrow(edge))) ||
      !isTRUE(all.equal(as.integer(fit$fit_dictionary_edge_count), nrow(edge)))) {
    stop(
      "Pando candidate and fitted dictionary counts must equal the deduplicated ",
      "common dictionary size.", call. = FALSE
    )
  }

  support_type <- as.character(support$source_type)
  if (any(!support_type %in% c("global", "condition"))) {
    stop("Pando dictionary support source_type is invalid.", call. = FALSE)
  }
  global_rows <- support_type == "global"
  local_rows <- support_type == "condition"
  if (any(!is.na(support$condition[global_rows]) &
          nzchar(as.character(support$condition[global_rows]))) ||
      anyNA(support$condition[local_rows]) ||
      any(!as.character(support$condition[local_rows]) %in% levels)) {
    stop("Pando global/local support condition labels are inconsistent.",
         call. = FALSE)
  }
  support_key <- paste(
    support$edge_id, support_type,
    ifelse(global_rows, "__global__", as.character(support$condition)),
    sep = "\001"
  )
  if (anyDuplicated(support_key) ||
      any(!as.character(support$edge_id) %in% as.character(edge$edge_id))) {
    stop("Pando dictionary support rows are duplicated or reference unknown edges.",
         call. = FALSE)
  }
  tf_threshold <- suppressWarnings(as.numeric(fit$candidate_tf_cor))
  peak_threshold <- suppressWarnings(as.numeric(fit$candidate_peak_cor))
  if (length(tf_threshold) != 1L || !is.finite(tf_threshold) ||
      tf_threshold < 0 || tf_threshold > 1 ||
      length(peak_threshold) != 1L || !is.finite(peak_threshold) ||
      peak_threshold < 0 || peak_threshold > 1 ||
      any(abs(as.numeric(support$tf_target_cor)) <= tf_threshold) ||
      any(abs(as.numeric(support$peak_target_cor)) <= peak_threshold)) {
    stop("Pando support rows do not satisfy the configured correlation gates.",
         call. = FALSE)
  }

  if (anyDuplicated(support_summary$edge_id) ||
      !setequal(as.character(support_summary$edge_id),
                as.character(edge$edge_id))) {
    stop("Pando dictionary support summary does not cover each exact edge once.",
         call. = FALSE)
  }
  summary_index <- match(edge$edge_id, support_summary$edge_id)
  if (anyNA(summary_index) ||
      !identical(as.logical(edge$source_global),
                 as.logical(support_summary$source_global[summary_index])) ||
      any(as.integer(support_summary$n_sources[summary_index]) < 1L)) {
    stop("Pando dictionary source provenance is inconsistent.", call. = FALSE)
  }
  edge_conditions <- .rc_split_support_conditions(edge$source_conditions)
  summary_conditions <- .rc_split_support_conditions(
    support_summary$support_conditions[summary_index]
  )
  if (any(vapply(seq_len(nrow(edge)), function(i) {
    !identical(edge_conditions[[i]], summary_conditions[[i]])
  }, logical(1)))) {
    stop("Pando condition-source provenance differs between dictionary tables.",
         call. = FALSE)
  }

  dictionary_ids <- sort(as.character(edge$edge_id))
  global_ids <- unique(as.character(support$edge_id[global_rows]))
  local_key <- paste(
    as.character(support$edge_id[local_rows]),
    as.character(support$condition[local_rows]), sep = "\001"
  )
  for (condition in levels) {
    one <- coefficient[
      as.character(coefficient$condition) == condition, , drop = FALSE
    ]
    if (nrow(one) != nrow(edge) || anyDuplicated(one$edge_id) ||
        !identical(sort(as.character(one$edge_id)), dictionary_ids)) {
      stop(
        "Every condition must contain every common-dictionary edge exactly once.",
        call. = FALSE
      )
    }
    valid <- one$estimable %in% TRUE & is.finite(as.numeric(one$pval))
    expected_padj <- rep(NA_real_, nrow(one))
    expected_padj[valid] <- stats::p.adjust(
      as.numeric(one$pval[valid]), method = "BH"
    )
    observed_padj <- as.numeric(one$padj)
    comparable <- is.finite(expected_padj) & is.finite(observed_padj)
    if (any(is.finite(expected_padj) != is.finite(observed_padj)) ||
        any(abs(expected_padj[comparable] - observed_padj[comparable]) > 1e-10)) {
      stop("Stored condition padj values do not equal BH adjustment.",
           call. = FALSE)
    }
    expected_statistical <- one$estimable %in% TRUE &
      is.finite(as.numeric(one$estimate)) &
      is.finite(observed_padj) & observed_padj < threshold
    expected_global <- as.character(one$edge_id) %in% global_ids
    expected_local <- paste(
      as.character(one$edge_id), condition, sep = "\001"
    ) %in% local_key
    expected_active <- expected_statistical & (expected_global | expected_local)
    if (!identical(as.logical(one$statistically_supported), expected_statistical) ||
        !identical(as.logical(one$global_support), expected_global) ||
        !identical(as.logical(one$local_support), expected_local) ||
        !identical(as.logical(one$active), expected_active) ||
        !identical(as.logical(one$significant), expected_active)) {
      stop(
        "Pando condition activity flags are inconsistent with ridge BH and ",
        "global/local candidate support.", call. = FALSE
      )
    }
    expected_effect <- ifelse(
      expected_active, as.numeric(one$estimate), 0
    )
    observed_effect <- as.numeric(one$penalty_effect)
    comparable_effect <- is.finite(expected_effect) & is.finite(observed_effect)
    if (any(is.finite(expected_effect) != is.finite(observed_effect)) ||
        any(abs(expected_effect[comparable_effect] -
                observed_effect[comparable_effect]) > 1e-12)) {
      stop("Pando penalty_effect does not match active condition coefficients.",
           call. = FALSE)
    }
  }

  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  required_fit <- c(
    "target", "condition", "rsq", "rsq_oof", "fit_status", "lambda",
    "predictor_scale_reference"
  )
  if (!all(required_fit %in% colnames(fit_table)) || !nrow(fit_table) ||
      any(as.character(fit_table$predictor_scale_reference) !=
          "equal_condition_within_condition_rms")) {
    stop("Pando target-level no-fusion ridge diagnostics are incomplete.",
         call. = FALSE)
  }
  fit_key <- paste(
    toupper(as.character(fit_table$target)),
    as.character(fit_table$condition), sep = "\001"
  )
  if (anyNA(fit_key) || anyDuplicated(fit_key) ||
      any(!as.character(fit_table$condition) %in% levels)) {
    stop("Pando target-level fit diagnostics are duplicated or mislabelled.",
         call. = FALSE)
  }

  contrast <- as.data.frame(fit$contrasts, stringsAsFactors = FALSE)
  required_contrast <- c(
    "edge_id", "target", "condition_a", "condition_b",
    "contrast_estimate", "contrast_se", "contrast_pval", "contrast_padj",
    "contrast_estimable", "contrast_significant"
  )
  if (!nrow(contrast) ||
      !all(required_contrast %in% colnames(contrast)) ||
      any(!contrast$condition_a %in% levels) ||
      any(!contrast$condition_b %in% levels) ||
      any(!contrast$edge_id %in% dictionary_ids)) {
    stop("Pando pairwise no-fusion condition-contrast table is incomplete.",
         call. = FALSE)
  }
  pair_key <- paste(contrast$condition_a, contrast$condition_b, sep = "\001")
  for (pair in unique(pair_key)) {
    index <- which(pair_key == pair)
    valid <- index[contrast$contrast_estimable[index] %in% TRUE &
                   is.finite(as.numeric(contrast$contrast_pval[index]))]
    expected <- rep(NA_real_, length(index))
    if (length(valid)) {
      expected[match(valid, index)] <- stats::p.adjust(
        as.numeric(contrast$contrast_pval[valid]), method = "BH"
      )
    }
    observed <- as.numeric(contrast$contrast_padj[index])
    comparable <- is.finite(expected) & is.finite(observed)
    if (any(is.finite(expected) != is.finite(observed)) ||
        any(abs(expected[comparable] - observed[comparable]) > 1e-10)) {
      stop("Stored condition-contrast padj values do not equal BH adjustment.",
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
    stop("Pando object metadata cannot validate condition fit cell mappings.",
         call. = FALSE)
  }
  for (fit in fits) {
    if (!identical(as.character(fit$condition_col), condition_col) ||
        !identical(as.character(fit$cell_type_col), celltype_col)) {
      stop("Pando fit metadata columns do not match the RegCompass request.",
           call. = FALSE)
    }
    for (condition in fit$condition_levels) {
      cells <- as.character(fit$condition_cell_ids[[condition]])
      missing <- setdiff(cells, rownames(metadata))
      if (length(missing)) {
        stop(
          "Pando fit references cells absent from its stored object; first ",
          "missing ID: ", missing[[1L]], ".", call. = FALSE
        )
      }
      observed_condition <- as.character(metadata[cells, condition_col])
      observed_celltype <- as.character(metadata[cells, celltype_col])
      if (anyNA(observed_condition) || anyNA(observed_celltype) ||
          any(observed_condition != condition) ||
          any(observed_celltype != as.character(fit$cell_type))) {
        stop(
          "Pando fit cell assignments disagree with stored object metadata ",
          "for cell type '", as.character(fit$cell_type),
          "' and condition '", condition, "'.", call. = FALSE
        )
      }
    }
  }

  rows <- list()
  universal <- list()
  diagnostics <- list()
  contrasts <- list()
  for (fit in fits) {
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
      "active_no_fusion_ridge_condition_coefficient_raw_tf_atac_units"
    coefficient$coefficient_contract <-
      "same_global_condition_union_dictionary_no_fusion_ridge_raw_units"
    coefficient$fit_engine <- fit$fit_engine
    coefficient$coefficient_scale <- fit$coefficient_scale
    coefficient$eligible_in_condition <- coefficient$active %in% TRUE
    coefficient$padj_threshold <- threshold

    fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
    fit_table$target <- toupper(as.character(fit_table$target))
    fit_table$condition <- as.character(fit_table$condition)
    fit_key <- paste(fit_table$target, fit_table$condition, sep = "\001")
    coef_key <- paste(coefficient$target, coefficient$condition, sep = "\001")
    fit_index <- match(coef_key, fit_key)
    if (anyNA(fit_index)) {
      stop("Condition coefficient and target fit diagnostics are misaligned.",
           call. = FALSE)
    }
    coefficient$rsq <- as.numeric(fit_table$rsq[fit_index])
    coefficient$fit_status <- as.character(fit_table$fit_status[fit_index])
    coefficient$reliable_model <- coefficient$fit_status == "ok"
    coefficient$penalty_eligible <-
      coefficient$active %in% TRUE & coefficient$fit_status == "ok" &
      is.finite(as.numeric(coefficient$penalty_effect))
    coefficient$active_in_condition <- coefficient$penalty_eligible
    rows[[length(rows) + 1L]] <- coefficient
    diagnostics[[length(diagnostics) + 1L]] <- fit_table

    edge <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
    shared_by_id <- stats::setNames(
      as.numeric(coefficient$shared_estimate[
        match(edge$edge_id, coefficient$edge_id)
      ]), edge$edge_id
    )
    summary <- edge
    summary$estimate <- unname(shared_by_id[summary$edge_id])
    summary$corr <- NA_real_
    summary[[celltype_col]] <- as.character(fit$cell_type)
    summary$summary_only <- TRUE
    summary$coefficient_contract <-
      "mean_condition_coefficient_summary_only_no_fusion_ridge"
    universal[[length(universal) + 1L]] <- summary

    contrast <- as.data.frame(fit$contrasts, stringsAsFactors = FALSE)
    contrast$tf <- toupper(as.character(contrast$tf))
    contrast$target <- toupper(as.character(contrast$target))
    contrast[[celltype_col]] <- as.character(fit$cell_type)
    contrast$fit_engine <- fit$fit_engine
    contrast$padj_threshold <- threshold
    contrasts[[length(contrasts) + 1L]] <- contrast
  }

  all_edges <- do.call(rbind, rows)
  rownames(all_edges) <- NULL
  active <- all_edges[all_edges$penalty_eligible %in% TRUE, , drop = FALSE]
  params <- methods::slot(methods::slot(grn_object, "grn"), "params")
  thresholds <- unique(as.numeric(all_edges$padj_threshold))
  list(
    network_index = params$condition_network_index %||% data.frame(),
    fit_diagnostics = do.call(rbind, diagnostics),
    fit_contracts = fits,
    universal = do.call(rbind, universal),
    condition_all = all_edges,
    condition_active = active,
    condition_effect_all = all_edges,
    condition_effect_active = active,
    condition_contrasts = do.call(rbind, contrasts),
    active_tol = 0,
    penalty_filter = paste0(
      "Pando active condition edge & fit_status == 'ok'; BH threshold=",
      paste(format(thresholds, trim = TRUE), collapse = ",")
    ),
    coefficient_contract =
      "same_global_condition_union_dictionary_no_fusion_ridge_raw_units",
    pando_fit_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    pando_memory_contract = "pando_native_condition_no_fusion_ridge"
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
      min_residual_df = 1L
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
    "rank_action", "min_residual_df", "rna_layer", "peak_layer",
    "peak_value_type", "condition_ridge_control"
  )
  unknown_infer_args <- setdiff(names(pando_infer_args), allowed_infer_args)
  if (length(unknown_infer_args)) {
    stop(
      "Unsupported `pando_infer_args`: ",
      paste(unknown_infer_args, collapse = ", "), call. = FALSE
    )
  }
  pando_infer_args <- utils::modifyList(list(
    tf_cor = 0.05, peak_cor = 0.05, adjust_method = "BH",
    padj_threshold = 0.05, rank_action = "mark", min_residual_df = 1L,
    rna_layer = "data", peak_layer = "data",
    peak_value_type = "normalized", condition_ridge_control = list()
  ), pando_infer_args)
  threshold <- suppressWarnings(as.numeric(pando_infer_args$padj_threshold))
  tf_threshold <- suppressWarnings(as.numeric(pando_infer_args$tf_cor))
  peak_threshold <- suppressWarnings(as.numeric(pando_infer_args$peak_cor))
  ridge_control <- pando_infer_args$condition_ridge_control
  if (!identical(toupper(as.character(pando_infer_args$adjust_method)), "BH") ||
      length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1 || !is.list(ridge_control) ||
      length(tf_threshold) != 1L || !is.finite(tf_threshold) ||
      tf_threshold < 0 || tf_threshold > 1 ||
      length(peak_threshold) != 1L || !is.finite(peak_threshold) ||
      peak_threshold < 0 || peak_threshold > 1) {
    stop(
      "Canonical RegCompass condition effects require BH adjustment with ",
      "padj_threshold in (0, 1), tf_cor/peak_cor in [0, 1], and ",
      "condition_ridge_control as a list.",
      call. = FALSE
    )
  }
  if ("fusion_ratio" %in% names(ridge_control)) {
    stop("`condition_ridge_control$fusion_ratio` is obsolete and unsupported.",
         call. = FALSE)
  }
  pando_infer_args$padj_threshold <- threshold
  pando_infer_args$tf_cor <- tf_threshold
  pando_infer_args$peak_cor <- peak_threshold

  condition_types <- if (is.null(cell_type)) {
    unique(as.character(object@meta.data[[celltype_col]]))
  } else {
    unique(as.character(cell_type))
  }
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
    "resolved Pando global-plus-condition common-dictionary ridge plan",
    current = 5L,
    context = list(
      cell_types = length(plans),
      tf_cor = pando_infer_args$tf_cor,
      peak_cor = pando_infer_args$peak_cor,
      padj_threshold = threshold,
      workers = worker_limit,
      estimator = .RC_PANDO_CONDITION_GRN_ENGINE,
      dictionary_policy = .RC_PANDO_CONDITION_DICTIONARY_POLICY
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
  .rc_step_monitor_event(
    progress_monitor, "condition_celltype_prepare",
    "initializing one Pando object per condition-GRN cell type",
    current = 6L,
    context = list(tasks = length(prepare_tasks))
  )
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
  .rc_step_monitor_event(
    progress_monitor, "condition_ridge_fit",
    "running Pando global-plus-condition union no-fusion ridge GRNs",
    current = 8L,
    context = list(
      tasks = length(fit_tasks),
      outer_celltype_parallel = outer_parallel,
      inner_target_parallel = inner_parallel,
      workers = worker_limit,
      padj_threshold = threshold
    )
  )
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
    stop("Pando condition ridge fits returned an invalid cell-type set.",
         call. = FALSE)
  }
  names(fit_results) <- fit_types
  fit_results <- fit_results[names(plans)]
  invisible(gc(verbose = FALSE, full = TRUE))

  parallel_plan <- list(
    scope = "cell_type_condition_ridge",
    cell_type_prepare_tasks = length(prepare_tasks),
    condition_ridge_fit_tasks = length(fit_tasks),
    workers = worker_limit,
    outer_celltype_parallel = outer_parallel,
    inner_target_parallel = inner_parallel,
    nested_parallel = FALSE,
    stage_barrier = paste(
      "Pando pooled/global candidate discovery + condition candidate discovery",
      "-> exact deduplicated union -> single common-lambda no-fusion ridge",
      "-> condition BH active-edge gate"
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
      key_frame <- data.frame(
        condition_value = condition,
        celltype_value = type,
        stringsAsFactors = FALSE
      )
      names(key_frame) <- c(condition_col, celltype_col)
      group_id <- rc_make_stratum_id(
        key_frame, c(condition_col, celltype_col)
      )
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
      status_rows[[length(status_rows) + 1L]] <- data.frame(
        group_id = group_id,
        condition_value = condition,
        celltype_value = type,
        n_cells = length(cells),
        status = "ok",
        n_target_genes = length(unique(all_rows$target)),
        n_edges = nrow(all_rows),
        n_active_edges = nrow(active_rows),
        padj_threshold = threshold,
        grn_evidence_role =
          "global_or_local_pando_candidate_single_no_fusion_ridge_active_edge",
        stringsAsFactors = FALSE
      )
      names(status_rows[[length(status_rows)]])[
        names(status_rows[[length(status_rows)]]) == "condition_value"
      ] <- condition_col
      names(status_rows[[length(status_rows)]])[
        names(status_rows[[length(status_rows)]]) == "celltype_value"
      ] <- celltype_col
    }
    status <- do.call(rbind, status_rows)
    rownames(status) <- NULL
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
          "Pando pooled/global plus condition-specific correlation candidate",
          "discovery, exact deduplicated common dictionary, single no-fusion",
          "common-lambda ridge fit"
        ),
        fit_dictionary = paste0(
          "exact edge union satisfying tf_cor=", pando_infer_args$tf_cor,
          " and peak_cor=", pando_infer_args$peak_cor,
          " in pooled/global cells or at least one condition"
        ),
        condition_effect = paste0(
          "raw-unit condition ridge coefficient active only when condition BH ",
          "padj < ", format(threshold, trim = TRUE),
          " and edge has global or same-condition Pando candidate support"
        ),
        coefficient_contract =
          "same_global_condition_union_dictionary_no_fusion_ridge_raw_units",
        ridge_inference = paste(
          "approximate BH-adjusted ridge-Wald coefficient and direct pairwise",
          "contrast inference conditional on candidate discovery and CV lambda"
        ),
        parallel_contract = parallel_plan,
        penalty_regulatory_evidence = paste(
          "active penalty_effect times metacell-mean TF times",
          "metacell-mean ATAC"
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
