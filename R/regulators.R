#' Rank candidate regulator-reaction associations
#'
#' This helper performs candidate prioritization, not causal driver discovery. It
#' combines available evidence ranks with a robust mean-rank score and annotates
#' coarse evidence tiers for downstream review.
#'
#' @param evidence Data frame with regulator/reaction IDs and numeric evidence
#' columns such as direct association, adjusted association, motif support,
#' enhancer support, or stability scores.
#' @param regulator_col Column containing regulator IDs.
#' @param reaction_col Column containing reaction IDs.
#' @param evidence_cols Numeric evidence columns to rank. Larger values are
#' treated as stronger support unless their names appear in `smaller_is_better`.
#' @param smaller_is_better Evidence columns where smaller values are stronger.
#'
#' @return Data frame sorted by candidate rank. The ranking is candidate
#' prioritization only and must not be interpreted as causal proof.
#' @export
rc_rank_regulators <- function(evidence,
                               regulator_col = "regulator_id",
                               reaction_col = "reaction_id",
                               evidence_cols = NULL,
                               smaller_is_better = character()) {
  if (!is.data.frame(evidence)) stop("`evidence` must be a data.frame.", call. = FALSE)
  required <- c(regulator_col, reaction_col)
  missing <- setdiff(required, colnames(evidence))
  if (length(missing) > 0L) stop("`evidence` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (is.null(evidence_cols)) {
    evidence_cols <- setdiff(names(evidence)[vapply(evidence, is.numeric, logical(1L))], character())
  }
  evidence_cols <- intersect(evidence_cols, colnames(evidence))
  if (length(evidence_cols) == 0L) stop("At least one numeric evidence column is required.", call. = FALSE)
  non_numeric <- evidence_cols[!vapply(evidence[evidence_cols], is.numeric, logical(1L))]
  if (length(non_numeric) > 0L) stop("Evidence columns must be numeric: ", paste(non_numeric, collapse = ", "), call. = FALSE)

  out <- evidence
  rank_cols <- paste0("rank_", evidence_cols)
  for (i in seq_along(evidence_cols)) {
    col <- evidence_cols[[i]]
    x <- out[[col]]
    if (col %in% smaller_is_better) {
      out[[rank_cols[[i]]]] <- rank(x, ties.method = "average", na.last = "keep")
    } else {
      out[[rank_cols[[i]]]] <- rank(-x, ties.method = "average", na.last = "keep")
    }
  }
  rank_mat <- as.matrix(out[rank_cols])
  out$rank_score <- matrixStats::rowMedians(rank_mat, na.rm = TRUE)
  out$n_evidence_available <- rowSums(!is.na(rank_mat))
  motif_vec <- if ("motif_support" %in% evidence_cols) !is.na(out$motif_support) & out$motif_support > 0 else rep(FALSE, nrow(out))
  enhancer_vec <- if ("enhancer_support" %in% evidence_cols) !is.na(out$enhancer_support) & out$enhancer_support > 0 else rep(FALSE, nrow(out))
  out$evidence_tier <- ifelse(motif_vec & enhancer_vec, "motif-and-enhancer-supported",
                              ifelse(motif_vec, "motif-supported", "correlation-only"))
  out$candidate_rank <- rank(out$rank_score, ties.method = "average", na.last = "keep")
  out <- out[order(out$candidate_rank, out[[reaction_col]], out[[regulator_col]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}
