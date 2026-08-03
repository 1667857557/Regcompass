# Integration of the Pando common-dictionary condition-GRN contract.

.RC_PANDO_CONDITION_GRN_FIT_SCHEMA <-
  "pando_condition_grn_common_dictionary_v1"

# Retained only because the Stage 1 dispatcher constructs this field before
# calling the condition implementation. The common-dictionary GLM has no native
# engine or memory planner.
.rc_pando_execution_summary <- function(diagnostics = NULL) {
  list(
    fit_engine = "two_stage_exact_edge_union_fixed_dictionary_glm",
    targets_total = if (is.data.frame(diagnostics)) {
      length(unique(as.character(diagnostics$target)))
    } else 0L,
    targets_failed = if (is.data.frame(diagnostics) &&
      "fit_status" %in% colnames(diagnostics)) {
      sum(!diagnostics$fit_status %in% c("ok", "rank_deficient"), na.rm = TRUE)
    } else 0L
  )
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
    "cell_type", "condition_levels", "condition_cell_ids",
    "edge_dictionary", "coefficients", "fit", "network_names",
    "padj_threshold", "adjust_method", "scale", "interaction",
    "projection_effect_column", "projection_policy", "rna_layer",
    "peak_layer", "peak_value_type", "preprocessing_fingerprint",
    "dictionary_preprocessing_provenance_verified"
  )
  if (!all(required %in% names(fit)) ||
      !identical(fit$scale, FALSE) ||
      !identical(fit$interaction, ":") ||
      !identical(fit$projection_effect_column, "penalty_effect") ||
      !identical(fit$projection_policy, "padj_significant_effects_only") ||
      !identical(toupper(as.character(fit$adjust_method)), "BH") ||
      !isTRUE(all.equal(as.numeric(fit$padj_threshold), 0.05)) ||
      !isTRUE(fit$dictionary_preprocessing_provenance_verified) ||
      any(!nzchar(c(
        as.character(fit$rna_layer), as.character(fit$peak_layer),
        as.character(fit$peak_value_type),
        as.character(fit$preprocessing_fingerprint)
      )))) {
    stop("Pando common-dictionary condition fit contract is incomplete.",
         call. = FALSE)
  }

  edge <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required_edge <- c(
    "edge_id", "target", "tf", "region", "atac_feature_id",
    "candidate_index"
  )
  required_coefficient <- c(
    "edge_id", "target", "tf", "region", "condition", "estimate",
    "std_err", "statistic", "pval", "padj", "significant",
    "penalty_effect", "estimable", "zero_variance", "aliased"
  )
  if (!nrow(edge) ||
      !all(required_edge %in% colnames(edge)) ||
      !all(required_coefficient %in% colnames(coefficient)) ||
      anyNA(edge$edge_id) || any(!nzchar(as.character(edge$edge_id))) ||
      anyDuplicated(edge$edge_id) ||
      any(!coefficient$condition %in% fit$condition_levels)) {
    stop("Pando common-dictionary coefficient table is incomplete.",
         call. = FALSE)
  }

  dictionary_ids <- sort(as.character(edge$edge_id))
  for (condition in fit$condition_levels) {
    one <- coefficient[
      as.character(coefficient$condition) == condition, , drop = FALSE
    ]
    if (nrow(one) != nrow(edge) || anyDuplicated(one$edge_id) ||
        !identical(sort(as.character(one$edge_id)), dictionary_ids)) {
      stop(
        "Every condition must contain every frozen dictionary edge exactly once.",
        call. = FALSE
      )
    }

    valid_p <- is.finite(as.numeric(one$pval))
    expected_padj <- rep(NA_real_, nrow(one))
    expected_padj[valid_p] <- stats::p.adjust(
      as.numeric(one$pval[valid_p]), method = "BH"
    )
    observed_padj <- as.numeric(one$padj)
    comparable <- is.finite(expected_padj) & is.finite(observed_padj)
    if (any(is.finite(expected_padj) != is.finite(observed_padj)) ||
        any(abs(expected_padj[comparable] - observed_padj[comparable]) > 1e-10)) {
      stop("Stored condition padj values do not equal BH adjustment.",
           call. = FALSE)
    }
  }

  expected_significant <- coefficient$estimable %in% TRUE &
    is.finite(as.numeric(coefficient$padj)) &
    as.numeric(coefficient$padj) < 0.05
  if (!identical(as.logical(coefficient$significant), expected_significant)) {
    stop("Pando significant-edge flags are not exactly estimable & padj < 0.05.",
         call. = FALSE)
  }

  expected_effect <- ifelse(
    expected_significant, as.numeric(coefficient$estimate), 0
  )
  observed_effect <- as.numeric(coefficient$penalty_effect)
  comparable_effect <- is.finite(expected_effect) & is.finite(observed_effect)
  if (any(is.finite(expected_effect) != is.finite(observed_effect)) ||
      any(abs(expected_effect[comparable_effect] -
              observed_effect[comparable_effect]) > 1e-12)) {
    stop("Pando penalty_effect does not match the strict BH edge contract.",
         call. = FALSE)
  }

  if ("direction" %in% colnames(coefficient)) {
    expected_direction <- ifelse(
      !coefficient$estimable, "undefined",
      ifelse(coefficient$estimate > 0, "positive",
             ifelse(coefficient$estimate < 0, "negative", "zero"))
    )
    if (!identical(as.character(coefficient$direction), expected_direction)) {
      stop("Pando coefficient directions are inconsistent with estimates.",
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

  rows <- list()
  universal <- list()
  diagnostics <- list()
  for (fit in fits) {
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
    coefficient$condition_effect <- coefficient$condition_estimate
    coefficient$effect_definition <-
      "fixed_dictionary_condition_glm_coefficient"
    coefficient$coefficient_contract <-
      "same_exact_edge_dictionary_unscaled_gaussian_glm"
    coefficient$fit_engine <- fit$fit_engine
    coefficient$coefficient_scale <- fit$coefficient_scale
    coefficient$eligible_in_condition <- coefficient$estimable

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
    coefficient$reliable_model <-
      coefficient$fit_status %in% c("ok", "rank_deficient")
    coefficient$penalty_eligible <-
      coefficient$estimable %in% TRUE &
      is.finite(as.numeric(coefficient$padj)) &
      as.numeric(coefficient$padj) < 0.05 &
      is.finite(as.numeric(coefficient$penalty_effect))
    coefficient$active_in_condition <- coefficient$penalty_eligible
    rows[[length(rows) + 1L]] <- coefficient
    diagnostics[[length(diagnostics) + 1L]] <- fit_table

    edge <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
    keys <- split(
      seq_len(nrow(coefficient)), as.character(coefficient$edge_id)
    )
    mean_or_na <- function(value) {
      value <- as.numeric(value)
      value <- value[is.finite(value)]
      if (length(value)) mean(value) else NA_real_
    }
    summary <- edge
    summary$estimate <- vapply(summary$edge_id, function(id) {
      mean_or_na(coefficient$condition_estimate[keys[[id]]])
    }, numeric(1))
    summary$corr <- NA_real_
    summary[[celltype_col]] <- as.character(fit$cell_type)
    summary$summary_only <- TRUE
    summary$coefficient_contract <-
      "same_exact_edge_dictionary_unscaled_gaussian_glm"
    universal[[length(universal) + 1L]] <- summary
  }

  all_edges <- do.call(rbind, rows)
  rownames(all_edges) <- NULL
  active <- all_edges[all_edges$penalty_eligible %in% TRUE, , drop = FALSE]
  params <- methods::slot(methods::slot(grn_object, "grn"), "params")
  list(
    network_index = params$condition_network_index %||% data.frame(),
    fit_diagnostics = do.call(rbind, diagnostics),
    fit_contracts = fits,
    universal = do.call(rbind, universal),
    condition_all = all_edges,
    condition_active = active,
    condition_effect_all = all_edges,
    condition_effect_active = active,
    active_tol = 0,
    penalty_filter =
      "estimable & BH-adjusted padj < 0.05; no effect-size or model-R2 gate",
    coefficient_contract =
      "same_exact_edge_dictionary_unscaled_gaussian_glm",
    pando_fit_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    pando_memory_contract = "not_applicable_no_native_condition_engine"
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
      tf_cor = 0.1, peak_cor = 0, adjust_method = "BH",
      padj_threshold = 0.05, rank_action = "mark",
      min_residual_df = 1L, parallel = FALSE
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
    "peak_value_type"
  )
  unknown_infer_args <- setdiff(names(pando_infer_args), allowed_infer_args)
  if (length(unknown_infer_args)) {
    stop(
      "Unsupported `pando_infer_args`: ",
      paste(unknown_infer_args, collapse = ", "), call. = FALSE
    )
  }
  defaults <- list(
    tf_cor = 0.1,
    peak_cor = 0,
    adjust_method = "BH",
    padj_threshold = 0.05,
    rank_action = "mark",
    min_residual_df = 1L
  )
  pando_infer_args <- utils::modifyList(defaults, pando_infer_args)
  if (!identical(toupper(as.character(pando_infer_args$adjust_method)), "BH") ||
      !isTRUE(all.equal(as.numeric(pando_infer_args$padj_threshold), 0.05))) {
    stop("Canonical RegCompass condition effects require BH padj < 0.05.",
         call. = FALSE)
  }
  .rc_step_monitor_event(
    progress_monitor, "condition_design",
    "configured exact-edge union and fixed-dictionary condition GLMs",
    current = 5L,
    context = list(
      tf_cor = pando_infer_args$tf_cor,
      peak_cor = pando_infer_args$peak_cor,
      adjust_method = "BH",
      padj_threshold = 0.05,
      scale = FALSE,
      interaction = ":"
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
  filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Pando condition GRNs"
  )
  object <- filtered$object
  init <- list(object = object, peak_assay = atac_assay, rna_assay = rna_assay)
  init[names(pando_initiate_args)] <- NULL
  grn <- do.call(Pando::initiate_grn, c(init, pando_initiate_args))
  pando_motif_args <- .rc_regcompass_motif_args(pando_motif_args)
  motif <- list(object = grn, pfm = pfm, genome = genome)
  motif[names(pando_motif_args)] <- NULL
  grn <- do.call(Pando::find_motifs, c(motif, pando_motif_args))

  infer <- list(
    object = grn,
    cell_type_col = celltype_col,
    condition_col = condition_col,
    cell_type = cell_type,
    genes = target_genes,
    network_name = "regcompass_condition_grn",
    min_cells_per_condition = as.integer(min_cells),
    small_condition_action = "error",
    tf_cor = pando_infer_args$tf_cor,
    peak_cor = pando_infer_args$peak_cor,
    adjust_method = "BH",
    padj_threshold = 0.05,
    rank_action = pando_infer_args$rank_action,
    min_residual_df = pando_infer_args$min_residual_df,
    parallel = FALSE,
    overwrite = TRUE,
    verbose = TRUE
  )
  .rc_step_monitor_event(
    progress_monitor, "fixed_dictionary_fit",
    "running global/condition discovery, exact union and condition GLMs",
    current = 9L, context = list(targets = length(target_genes))
  )
  grn <- do.call(Pando::infer_condition_grn, infer)
  extracted <- .rc_extract_condition_grn_contract(
    grn, condition_col, celltype_col
  )
  execution_summary <- .rc_pando_execution_summary(
    extracted$fit_diagnostics
  )

  meta <- object@meta.data
  status_rows <- list()
  for (fit in extracted$fit_contracts) {
    for (condition in fit$condition_levels) {
      cells <- fit$condition_cell_ids[[condition]]
      key_frame <- data.frame(
        condition_value = condition,
        celltype_value = fit$cell_type,
        stringsAsFactors = FALSE
      )
      names(key_frame) <- c(condition_col, celltype_col)
      id <- rc_make_stratum_id(
        key_frame, c(condition_col, celltype_col)
      )
      all_rows <- extracted$condition_all[
        extracted$condition_all[[condition_col]] == condition &
          extracted$condition_all[[celltype_col]] == fit$cell_type,
        , drop = FALSE
      ]
      active_rows <- extracted$condition_active[
        extracted$condition_active[[condition_col]] == condition &
          extracted$condition_active[[celltype_col]] == fit$cell_type,
        , drop = FALSE
      ]
      one <- data.frame(
        group_id = id,
        condition_value = condition,
        celltype_value = fit$cell_type,
        n_cells = length(cells),
        status = "ok",
        n_target_genes = length(unique(all_rows$target)),
        n_edges = nrow(all_rows),
        n_active_edges = nrow(active_rows),
        grn_evidence_role =
          "within_cell_type_common_dictionary_condition_glm",
        stringsAsFactors = FALSE
      )
      names(one)[names(one) == "condition_value"] <- condition_col
      names(one)[names(one) == "celltype_value"] <- celltype_col
      status_rows[[length(status_rows) + 1L]] <- one
    }
  }
  status <- do.call(rbind, status_rows)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(status, file.path(outdir, "pando_group_status.tsv.gz"))
  .rc_write_tsv_gz(extracted$condition_all,
    file.path(outdir, "pando_tf_peak_gene_condition_all.tsv.gz"))
  .rc_write_tsv_gz(extracted$condition_active,
    file.path(outdir, "pando_tf_peak_gene_condition_active.tsv.gz"))
  .rc_write_tsv_gz(extracted$universal,
    file.path(outdir, "pando_tf_peak_gene_universal.tsv.gz"))
  saveRDS(extracted$fit_contracts,
    file.path(outdir, "pando_condition_grn_fits.rds"))
  if (isTRUE(save_pando_objects)) {
    dir.create(file.path(outdir, "pando_objects"), recursive = TRUE,
               showWarnings = FALSE)
    saveRDS(grn, file.path(outdir, "pando_objects", "condition_grn_fit.rds"))
  }
  selected_cells <- unique(unlist(lapply(
    extracted$fit_contracts, function(fit) {
      unlist(fit$condition_cell_ids, use.names = FALSE)
    }
  ), use.names = FALSE))
  answer <- list(
    schema_version = "regcompass_condition_grn_common_dictionary_v1",
    analysis_mode = "condition_grn",
    condition_coefficients_calculated = TRUE,
    pando_fit_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
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
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF across conditions",
      grn_fit =
        "global-plus-condition candidate discovery, exact edge union, fixed-dictionary Gaussian GLM",
      condition_effect =
        "unscaled fixed-dictionary condition coefficient",
      coefficient_contract =
        "same_exact_edge_dictionary_unscaled_gaussian_glm",
      significance = "estimable and BH adjusted P below 0.05",
      penalty_regulatory_evidence =
        "paired-cell TF-by-ATAC projection using penalty_effect without effect-size or model-R2 gates"
    ),
    group_cols = c(condition_col, celltype_col)
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
