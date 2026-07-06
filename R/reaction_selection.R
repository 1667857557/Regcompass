#' Select candidate reactions for RegCompassR Layer 2
#' @export
rc_select_layer2_reactions <- function(layer1, gem, selected_reactions = NULL,
                                       selection_method = c("auto", "top", "differential", "pathway", "custom"),
                                       top_n = 300, min_C_rel = 0.15, min_confidence = 0.25,
                                       neighbor_depth = 1, max_subgem_reactions = 1000) {
  selection_method <- match.arg(selection_method)
  valid <- rc_validate_gem(gem)
  rxns <- valid$reactions
  if (!is.null(selected_reactions)) {
    keep <- intersect(unique(as.character(selected_reactions)), rxns)
    reason <- stats::setNames(rep("custom", length(keep)), keep)
  } else {
    C <- as.matrix(layer1$C_rel)
    Conf <- if (!is.null(layer1$reaction_confidence)) rc_layer2_confidence_matrix(layer1$reaction_confidence, C) else matrix(1, nrow(C), ncol(C), dimnames = dimnames(C))
    common <- intersect(intersect(rownames(C), rownames(Conf)), rxns)
    invalid <- rc_layer2_invalid_reactions(layer1)
    common <- setdiff(common, invalid)
    evidence <- pmax(rowMedians_safe(C[common, , drop = FALSE]), rowMedians_safe(Conf[common, , drop = FALSE]), na.rm = TRUE)
    pass <- common[(rowMedians_safe(C[common, , drop = FALSE]) >= min_C_rel) | (rowMedians_safe(Conf[common, , drop = FALSE]) >= min_confidence)]
    ranked <- names(sort(evidence[pass], decreasing = TRUE, na.last = NA))
    keep <- utils::head(ranked, top_n)
    reason <- stats::setNames(rep("Layer1 C_rel/confidence threshold or top evidence", length(keep)), keep)
  }
  keep <- rc_add_reaction_neighbors(valid$S, keep, depth = neighbor_depth, limit = max_subgem_reactions)
  reason[setdiff(keep, names(reason))] <- "shared-metabolite neighbor/support reaction"
  keep <- utils::head(unique(keep), max_subgem_reactions)
  data.frame(reaction_id = keep, reaction_selection_reason = unname(reason[keep]), stringsAsFactors = FALSE)
}

rc_add_reaction_neighbors <- function(S, seeds, depth = 1, limit = 1000) {
  seeds <- intersect(seeds, colnames(S))
  keep <- seeds
  for (d in seq_len(max(0, depth))) {
    mets <- rownames(S)[rowSums(abs(S[, keep, drop = FALSE]) > 0) > 0]
    nbr <- colnames(S)[colSums(abs(S[mets, , drop = FALSE]) > 0) > 0]
    keep <- unique(c(keep, nbr))
    if (length(keep) >= limit) return(utils::head(keep, limit))
  }
  keep
}

rowMedians_safe <- function(x) {
  x <- as.matrix(x)
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowMedians(x, na.rm = TRUE) else apply(x, 1, stats::median, na.rm = TRUE)
}

rc_layer2_invalid_reactions <- function(layer1) {
  invalid <- character()
  q95 <- layer1$q95_diagnostics
  if (is.data.frame(q95) && "reaction_id" %in% colnames(q95)) {
    if ("all_missing_reaction_flag" %in% colnames(q95)) invalid <- union(invalid, as.character(q95$reaction_id[q95$all_missing_reaction_flag %in% TRUE]))
    if ("q95_power_class" %in% colnames(q95)) invalid <- union(invalid, as.character(q95$reaction_id[as.character(q95$q95_power_class) == "very_low"]))
  }
  conf <- layer1$reaction_confidence
  if (is.data.frame(conf) && "reaction_id" %in% colnames(conf)) {
    flag_col <- if ("reaction_unsupported_by_complete_gpr_flag" %in% colnames(conf)) "reaction_unsupported_by_complete_gpr_flag" else if ("no_complete_gpr_group_flag" %in% colnames(conf)) "no_complete_gpr_group_flag" else NA_character_
    if (!is.na(flag_col)) invalid <- union(invalid, as.character(conf$reaction_id[conf[[flag_col]] %in% TRUE]))
  }
  if (!is.null(layer1$C_rel)) {
    C <- as.matrix(layer1$C_rel)
    invalid <- union(invalid, rownames(C)[rowSums(is.finite(C)) == 0])
  }
  invalid
}
