# API-based Pando compatibility and RegCompass penalty-entry contracts.

# Retained only because older Stage-1 provenance contains these fields. They are
# no longer used as post-fit edge gates; zero records the absence of an extra
# correlation or coefficient threshold beyond Pando candidate screening and BH.
.RC_PANDO_PENALTY_CORR_THRESHOLD <- 0
.RC_PANDO_PENALTY_ESTIMATE_THRESHOLD <- 0

.rc_condition_penalty_gate <- function(coefficient) {
  required <- c("estimate", "padj", "estimable")
  if (!is.data.frame(coefficient) ||
      !all(required %in% colnames(coefficient))) {
    stop(
      "Condition-GRN coefficients must contain estimate, padj, and estimable ",
      "columns before RegCompass penalty filtering.",
      call. = FALSE
    )
  }
  estimate <- suppressWarnings(as.numeric(coefficient$estimate))
  padj <- suppressWarnings(as.numeric(coefficient$padj))
  coefficient$estimable %in% TRUE &
    is.finite(padj) & padj < 0.05 &
    is.finite(estimate)
}

.rc_apply_condition_penalty_gate <- function(fit) {
  .rc_require_pando_condition_grn_fit(fit)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  gate <- .rc_condition_penalty_gate(coefficient)
  coefficient$significant <- gate
  coefficient$penalty_effect <- ifelse(
    gate, suppressWarnings(as.numeric(coefficient$estimate)), 0
  )
  fit$coefficients <- coefficient
  fit$regcompass_penalty_filter <- "estimable & BH padj < 0.05"
  fit
}

.rc_filter_standard_pando_edges <- function(table) {
  required <- c("estimate", "padj")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop(
      "Standard Pando requires estimate and padj columns for RegCompass ",
      "penalty filtering.",
      call. = FALSE
    )
  }
  estimate <- suppressWarnings(as.numeric(table$estimate))
  padj <- suppressWarnings(as.numeric(table$padj))
  keep <- is.finite(estimate) &
    is.finite(padj) & padj < .rc_standard_pando_padj_fixed
  if ("estimable" %in% colnames(table)) {
    keep <- keep & table$estimable %in% TRUE
  }
  answer <- table[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    estimable = if ("estimable" %in% colnames(table)) TRUE else NA,
    padj = "< 0.05"
  )
  answer
}

.rc_merge_pando_results <- function(
    condition_result = NULL, standard_results = list(),
    condition_types = character(), standard_types = character(),
    condition_col, celltype_col, outdir) {
  answer <- .rc_merge_pando_results_with_parallel_objects(
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
    "estimable and BH adjusted P below 0.05"
  answer$normalization_policy$standard_edge_filter <-
    "estimable when available and adjusted P below 0.05"

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
