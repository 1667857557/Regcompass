#' Filter Layer 1 reactions to valid downstream candidates
#'
#' Removes reactions that are all-missing after Q95 calibration, optionally have
#' very-low Q95 power, lack any complete GPR-supporting group, or have no finite
#' reaction confidence.
#' @export
rc_filter_valid_reactions <- function(layer1,
                                      require_complete_gpr = TRUE,
                                      exclude_q95_very_low = TRUE,
                                      require_finite_confidence = FALSE) {
  if (is.null(layer1$C_rel) || is.null(layer1$C_raw)) stop("`layer1` must contain `C_rel` and `C_raw`.", call. = FALSE)
  C_rel <- as.matrix(layer1$C_rel)
  C_raw <- as.matrix(layer1$C_raw)
  invalid <- character(0)

  q95 <- layer1$q95_diagnostics
  if (!is.null(q95) && nrow(q95) > 0L) {
    if ("all_missing_reaction_flag" %in% colnames(q95)) {
      invalid <- union(invalid, as.character(q95$reaction_id[q95$all_missing_reaction_flag %in% TRUE]))
    }
    if (isTRUE(exclude_q95_very_low) && "q95_power_class" %in% colnames(q95)) {
      invalid <- union(invalid, as.character(q95$reaction_id[as.character(q95$q95_power_class) == "very_low"]))
    }
  }

  conf <- layer1$reaction_confidence
  if (!is.null(conf) && nrow(conf) > 0L) {
    unsupported_col <- if ("reaction_unsupported_by_complete_gpr_flag" %in% colnames(conf)) {
      "reaction_unsupported_by_complete_gpr_flag"
    } else if ("no_complete_gpr_group_flag" %in% colnames(conf)) {
      "no_complete_gpr_group_flag"
    } else {
      NA_character_
    }
    if (isTRUE(require_complete_gpr) && !is.na(unsupported_col)) {
      invalid <- union(invalid, as.character(conf$reaction_id[conf[[unsupported_col]] %in% TRUE]))
    }
    if (isTRUE(require_finite_confidence) && "reaction_confidence" %in% colnames(conf)) {
      med_conf <- tapply(conf$reaction_confidence, conf$reaction_id, stats::median, na.rm = TRUE)
      invalid <- union(invalid, names(med_conf)[!is.finite(med_conf)])
    }
  }

  valid <- setdiff(rownames(C_rel), invalid)
  list(
    valid_reactions = valid,
    invalid_reactions = invalid,
    C_rel = C_rel[valid, , drop = FALSE],
    C_raw = C_raw[valid, , drop = FALSE]
  )
}

#' Rank Layer 1 reactions after optional validity filtering
#' @export
rc_rank_reactions <- function(layer1,
                              score = c("median_C_rel", "iqr_C_rel"),
                              filter_valid = TRUE) {
  score <- match.arg(score)
  x <- if (isTRUE(filter_valid)) rc_filter_valid_reactions(layer1)$C_rel else as.matrix(layer1$C_rel)
  stat <- switch(
    score,
    median_C_rel = matrixStats::rowMedians(x, na.rm = TRUE),
    iqr_C_rel = matrixStats::rowIQRs(x, na.rm = TRUE)
  )
  out <- data.frame(reaction_id = rownames(x), score = stat, stringsAsFactors = FALSE)
  out[order(out$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
}
