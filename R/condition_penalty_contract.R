# API-based Pando compatibility and RegCompass penalty-entry gates.

.RC_PANDO_PENALTY_CORR_THRESHOLD <- 0.05
.RC_PANDO_PENALTY_ESTIMATE_THRESHOLD <- 0.05

.rc_attach_condition_penalty_corr <- function(fit) {
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  if (!"edge_id" %in% colnames(coefficient)) {
    stop("Condition-GRN coefficients lack edge_id.", call. = FALSE)
  }
  had_native_corr <- "corr" %in% colnames(coefficient)
  native_corr <- if (had_native_corr) {
    suppressWarnings(as.numeric(coefficient$corr))
  } else {
    rep(NA_real_, nrow(coefficient))
  }
  corr <- native_corr
  missing_corr <- !is.finite(corr)
  if (any(missing_corr)) {
    dictionary <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
    if (!all(c("edge_id", "max_abs_tf_target_cor") %in%
             colnames(dictionary))) {
      stop(
        "Condition-GRN fit lacks both coefficient corr and dictionary ",
        "max_abs_tf_target_cor required by the RegCompass penalty gate.",
        call. = FALSE
      )
    }
    dictionary_corr <- suppressWarnings(as.numeric(
      dictionary$max_abs_tf_target_cor[
        match(coefficient$edge_id, dictionary$edge_id)
      ]
    ))
    corr[missing_corr] <- dictionary_corr[missing_corr]
  }
  coefficient$corr <- corr
  coefficient$corr_source <- ifelse(
    had_native_corr & is.finite(native_corr),
    "condition_coefficient_corr",
    "edge_dictionary_max_abs_tf_target_cor"
  )
  fit$coefficients <- coefficient
  fit
}

.rc_condition_penalty_gate <- function(coefficient) {
  required <- c("estimate", "padj", "corr", "estimable")
  if (!is.data.frame(coefficient) ||
      !all(required %in% colnames(coefficient))) {
    stop(
      "Condition-GRN coefficients must contain estimate, padj, corr, and estimable ",
      "columns before RegCompass penalty filtering.",
      call. = FALSE
    )
  }
  estimate <- suppressWarnings(as.numeric(coefficient$estimate))
  padj <- suppressWarnings(as.numeric(coefficient$padj))
  corr <- suppressWarnings(as.numeric(coefficient$corr))
  coefficient$estimable %in% TRUE &
    is.finite(padj) & padj < 0.05 &
    is.finite(corr) & abs(corr) >= .RC_PANDO_PENALTY_CORR_THRESHOLD &
    is.finite(estimate) &
    abs(estimate) >= .RC_PANDO_PENALTY_ESTIMATE_THRESHOLD
}

.rc_apply_condition_penalty_gate <- function(fit) {
  fit <- .rc_attach_condition_penalty_corr(fit)
  validation_fit <- fit
  if (!is.null(validation_fit$regcompass_penalty_filter)) {
    validation_coefficient <- as.data.frame(
      validation_fit$coefficients, stringsAsFactors = FALSE
    )
    legacy_gate <- validation_coefficient$estimable %in% TRUE &
      is.finite(suppressWarnings(as.numeric(validation_coefficient$padj))) &
      suppressWarnings(as.numeric(validation_coefficient$padj)) < 0.05
    validation_coefficient$significant <- legacy_gate
    validation_coefficient$penalty_effect <- ifelse(
      legacy_gate,
      suppressWarnings(as.numeric(validation_coefficient$estimate)),
      0
    )
    validation_fit$coefficients <- validation_coefficient
  }
  .rc_require_pando_condition_grn_fit(validation_fit)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  gate <- .rc_condition_penalty_gate(coefficient)
  coefficient$significant <- gate
  coefficient$penalty_effect <- ifelse(
    gate, suppressWarnings(as.numeric(coefficient$estimate)), 0
  )
  fit$coefficients <- coefficient
  fit$regcompass_penalty_filter <-
    "estimable & BH padj < 0.05 & abs(corr) >= 0.05 & abs(estimate) >= 0.05"
  fit$regcompass_corr_threshold <- .RC_PANDO_PENALTY_CORR_THRESHOLD
  fit$regcompass_estimate_threshold <-
    .RC_PANDO_PENALTY_ESTIMATE_THRESHOLD
  fit
}

.rc_filter_standard_pando_edges <- function(table) {
  required <- c("estimate", "padj", "corr")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop(
      "Standard Pando requires estimate, padj, and corr columns for RegCompass ",
      "penalty filtering.",
      call. = FALSE
    )
  }
  estimate <- suppressWarnings(as.numeric(table$estimate))
  padj <- suppressWarnings(as.numeric(table$padj))
  corr <- suppressWarnings(as.numeric(table$corr))
  keep <- is.finite(estimate) &
    abs(estimate) >= .RC_PANDO_PENALTY_ESTIMATE_THRESHOLD &
    is.finite(corr) & abs(corr) >= .RC_PANDO_PENALTY_CORR_THRESHOLD &
    is.finite(padj) & padj < .rc_standard_pando_padj_fixed
  if ("estimable" %in% colnames(table)) {
    keep <- keep & table$estimable %in% TRUE
  }
  answer <- table[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    padj = "< 0.05",
    absolute_correlation = ">= 0.05",
    absolute_estimate = ">= 0.05"
  )
  answer
}

.rc_merge_pando_results_penalty_base <- .rc_merge_pando_results

.rc_merge_pando_results <- function(
    condition_result = NULL, standard_results = list(),
    condition_types = character(), standard_types = character(),
    condition_col, celltype_col, outdir) {
  answer <- .rc_merge_pando_results_penalty_base(
    condition_result = condition_result,
    standard_results = standard_results,
    condition_types = condition_types,
    standard_types = standard_types,
    condition_col = condition_col,
    celltype_col = celltype_col,
    outdir = outdir
  )

  if (length(answer$condition_grn_fits)) {
    answer$condition_grn_fits <- lapply(
      answer$condition_grn_fits, .rc_apply_condition_penalty_gate
    )
  }
  if (is.data.frame(answer$tf_peak_gene_condition_all) &&
      nrow(answer$tf_peak_gene_condition_all)) {
    all_edges <- answer$tf_peak_gene_condition_all
    condition_rows <- rep(TRUE, nrow(all_edges))
    if ("analysis_mode" %in% colnames(all_edges)) {
      mode <- as.character(all_edges$analysis_mode)
      condition_rows <- is.na(mode) | !nzchar(mode) | mode == "condition_grn"
    }
    condition_edges <- all_edges[condition_rows, , drop = FALSE]
    if (nrow(condition_edges)) {
      corr_rows <- .rc_bind_frames_fill(lapply(
        answer$condition_grn_fits, function(fit) {
          coefficient <- as.data.frame(fit$coefficients,
                                       stringsAsFactors = FALSE)
          data.frame(
            .cell_type = as.character(fit$cell_type),
            condition = as.character(coefficient$condition),
            edge_id = as.character(coefficient$edge_id),
            corr = suppressWarnings(as.numeric(coefficient$corr)),
            corr_source = as.character(coefficient$corr_source),
            stringsAsFactors = FALSE
          )
        }
      ))
      condition_key <- paste(
        as.character(condition_edges[[celltype_col]]),
        as.character(condition_edges$condition),
        as.character(condition_edges$edge_id),
        sep = "\001"
      )
      lookup_key <- paste(
        corr_rows$.cell_type, corr_rows$condition, corr_rows$edge_id,
        sep = "\001"
      )
      lookup_index <- match(condition_key, lookup_key)
      if (anyNA(lookup_index)) {
        stop("Condition-GRN correlation lookup is incomplete.", call. = FALSE)
      }
      condition_edges$corr <- corr_rows$corr[lookup_index]
      condition_edges$corr_source <- corr_rows$corr_source[lookup_index]
      gate <- .rc_condition_penalty_gate(condition_edges)
      condition_edges$significant <- gate
      condition_edges$penalty_effect <- ifelse(
        gate, suppressWarnings(as.numeric(condition_edges$estimate)), 0
      )
      condition_edges$penalty_eligible <- gate
      condition_edges$active_in_condition <- gate
      all_edges[condition_rows, colnames(condition_edges)] <- condition_edges
      active_condition <- condition_edges[gate, , drop = FALSE]
    } else {
      active_condition <- data.frame()
    }
    existing_active <- answer$tf_peak_gene_condition
    active_standard <- if (is.data.frame(existing_active) &&
        nrow(existing_active) &&
        "analysis_mode" %in% colnames(existing_active)) {
      existing_active[
        as.character(existing_active$analysis_mode) == "standard_pando",
        , drop = FALSE
      ]
    } else {
      data.frame()
    }
    active <- .rc_bind_frames_fill(list(active_condition, active_standard))
    answer$tf_peak_gene_condition_all <- all_edges
    answer$tf_peak_gene_condition <- active
    answer$tf_peak_gene_condition_effect_all <- all_edges
    answer$tf_peak_gene_condition_effect <- active
  }
  answer$normalization_policy$condition_effect_filter <-
    "estimable, BH adjusted P below 0.05, absolute corr at least 0.05, and absolute estimate at least 0.05"
  answer$normalization_policy$standard_edge_filter <-
    "adjusted P below 0.05, absolute corr at least 0.05, and absolute estimate at least 0.05"

  .rc_write_tsv_gz(
    answer$tf_peak_gene_condition_all,
    file.path(outdir, "pando_tf_peak_gene_all.tsv.gz")
  )
  .rc_write_tsv_gz(
    answer$tf_peak_gene_condition,
    file.path(outdir, "pando_tf_peak_gene_active.tsv.gz")
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
