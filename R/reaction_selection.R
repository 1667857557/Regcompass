#' Select candidate reactions for RegCompassR Layer 2
#' @export
rc_select_layer2_reactions <- function(layer1, gem, selected_reactions = NULL,
                                       selection_method = c("auto", "top", "differential", "pathway", "custom"),
                                       top_n = 300, min_C_rel = 0.15, min_confidence = 0.25,
                                       neighbor_depth = 1, max_subgem_reactions = 1000,
                                       override_invalid = FALSE) {
  selection_method <- match.arg(selection_method)
  valid <- rc_validate_gem(gem)
  rxns <- valid$reactions
  invalid <- rc_layer2_invalid_reactions(layer1)
  custom_invalid_reaction_warning <- NULL
  if (!is.null(selected_reactions)) {
    requested <- intersect(unique(as.character(selected_reactions)), rxns)
    invalid_requested <- intersect(requested, invalid)
    if (length(invalid_requested) && !isTRUE(override_invalid)) {
      custom_invalid_reaction_warning <- paste("Filtered invalid custom reactions:", paste(invalid_requested, collapse = ", "))
      warning(custom_invalid_reaction_warning, call. = FALSE)
    }
    keep <- if (isTRUE(override_invalid)) requested else setdiff(requested, invalid)
    reason <- stats::setNames(rep(if (isTRUE(override_invalid)) "custom (invalid filtering overridden)" else "custom", length(keep)), keep)
  } else {
    C <- as.matrix(layer1$C_rel)
    Conf <- if (!is.null(layer1$reaction_confidence)) rc_layer2_confidence_matrix(layer1$reaction_confidence, C) else matrix(1, nrow(C), ncol(C), dimnames = dimnames(C))
    common <- intersect(intersect(rownames(C), rownames(Conf)), rxns)
    common <- setdiff(common, invalid)
    evidence <- pmax(rowMedians_safe(C[common, , drop = FALSE]), rowMedians_safe(Conf[common, , drop = FALSE]), na.rm = TRUE)
    pass <- common[(rowMedians_safe(C[common, , drop = FALSE]) >= min_C_rel) | (rowMedians_safe(Conf[common, , drop = FALSE]) >= min_confidence)]
    ranked <- names(sort(evidence[pass], decreasing = TRUE, na.last = NA))
    keep <- utils::head(ranked, top_n)
    reason <- stats::setNames(rep("Layer1 C_rel/confidence threshold or top evidence", length(keep)), keep)
  }
  keep <- rc_add_reaction_neighbors(valid$S, keep, depth = neighbor_depth, limit = max_subgem_reactions)
  if (!isTRUE(override_invalid)) keep <- setdiff(keep, invalid)
  reason[setdiff(keep, names(reason))] <- "shared-metabolite neighbor/support reaction"
  keep <- utils::head(unique(keep), max_subgem_reactions)
  out <- data.frame(reaction_id = keep, reaction_selection_reason = unname(reason[keep]), stringsAsFactors = FALSE)
  attr(out, "custom_invalid_reaction_warning") <- custom_invalid_reaction_warning
  out
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

#' Select target reactions without expanding the network
#' @export
rc_select_target_reactions <- function(layer1, targets = NULL, pathway = NULL, subsystem = NULL, method = c("balanced_top_capacity", "custom", "top_capacity", "condition_specific", "differential", "pathway"), selection_mode = c("balanced_rank", "and", "or", "ranked_product"), group_cols = c("condition", "cell_type"), top_n = 100, top_n_per_group = 30, min_C_rel = 0.15, min_confidence = 0.25, min_units_per_group = 3, require_complete_gpr = TRUE, exclude_very_low_q95_power = TRUE, exclude_low_q95_power = exclude_very_low_q95_power) {
  method <- match.arg(method)
  selection_mode <- match.arg(selection_mode)
  invalid <- rc_target_selection_invalid_reactions(layer1, require_complete_gpr = require_complete_gpr, exclude_low_q95_power = exclude_low_q95_power)
  if (!is.null(targets) || method == "custom") {
    keep <- setdiff(unique(as.character(targets)), invalid)
    return(data.frame(reaction_id = keep, selection_reason = "custom", require_complete_gpr = require_complete_gpr, stringsAsFactors = FALSE))
  }
  C <- as.matrix(layer1$C_rel); conf <- if (!is.null(layer1$reaction_confidence)) rc_layer2_confidence_matrix(layer1$reaction_confidence, C) else matrix(1, nrow(C), ncol(C), dimnames=dimnames(C))
  meta <- layer1$unit_meta %||% layer1$pool_meta %||% layer1$metacell_meta %||% data.frame(unit_id = colnames(C), stringsAsFactors = FALSE)
  id_col <- .rc_first_existing_col(c("unit_id", "pool_id", "sample_celltype_id", "metacell_id"), meta, fallback = colnames(meta)[1])
  meta <- meta[match(colnames(C), as.character(meta[[id_col]])), , drop = FALSE]
  groups <- if (all(group_cols %in% colnames(meta))) interaction(meta[, group_cols, drop = FALSE], drop = TRUE, sep = "|") else factor(rep("all", ncol(C)))
  rxns <- rownames(C); gs <- matrix(NA_real_, nrow(C), length(levels(groups)), dimnames = list(rxns, levels(groups)))
  n_by <- setNames(integer(length(levels(groups))), levels(groups))
  medC <- medConf <- gs
  for (g in levels(groups)) {
    cols <- which(groups == g); n_by[g] <- length(cols)
    medC[, g] <- rowMedians_safe(C[, cols, drop = FALSE])
    medConf[, g] <- rowMedians_safe(conf[, cols, drop = FALSE])
    gs[, g] <- sqrt(pmax(0, medC[, g]) * pmax(0, medConf[, g]))
    gs[n_by[g] < min_units_per_group, g] <- NA_real_
  }
  balanced_score <- apply(gs, 1, max, na.rm = TRUE); balanced_score[!is.finite(balanced_score)] <- NA_real_
  pass <- switch(selection_mode,
                 and = rowSums(medC >= min_C_rel & medConf >= min_confidence, na.rm = TRUE) > 0,
                 or = rowSums(medC >= min_C_rel | medConf >= min_confidence, na.rm = TRUE) > 0,
                 ranked_product = is.finite(balanced_score),
                 balanced_rank = is.finite(balanced_score))
  keep <- setdiff(rxns[pass], invalid)
  ranked <- names(sort(balanced_score[keep], decreasing=TRUE, na.last=NA))
  selected_groups <- if (length(ranked)) apply(gs[ranked,,drop=FALSE], 1, function(x) paste(names(sort(x, decreasing=TRUE, na.last=NA))[seq_len(min(top_n_per_group, sum(is.finite(x))))], collapse = ";")) else character(0)
  keep <- utils::head(ranked, top_n)
  if (length(keep) == 0L) {
    return(data.frame(reaction_id=character(), selection_reason=character(), selected_in_group=character(), balanced_score=numeric(), median_C_rel_by_group=character(), median_confidence_by_group=character(), n_units_by_group=character(), low_power_selection_flag=logical(), gpr_complexity=character(), require_complete_gpr=logical(), stringsAsFactors=FALSE))
  }
  data.frame(reaction_id=keep, selection_reason=paste(method, selection_mode, sep=":"), selected_in_group=unname(selected_groups[keep]), balanced_score=balanced_score[keep], median_C_rel_by_group=apply(round(medC[keep,,drop=FALSE],4),1,paste,collapse=";"), median_confidence_by_group=apply(round(medConf[keep,,drop=FALSE],4),1,paste,collapse=";"), n_units_by_group=paste(n_by, collapse=";"), low_power_selection_flag=any(n_by < min_units_per_group), gpr_complexity=NA_character_, require_complete_gpr=require_complete_gpr, stringsAsFactors=FALSE)
}

rc_target_selection_invalid_reactions <- function(layer1, require_complete_gpr = TRUE, exclude_low_q95_power = TRUE) {
  invalid <- character(0)
  q95 <- layer1$q95_diagnostics
  if (is.data.frame(q95) && "reaction_id" %in% colnames(q95)) {
    if ("all_missing_reaction_flag" %in% colnames(q95)) invalid <- union(invalid, as.character(q95$reaction_id[q95$all_missing_reaction_flag %in% TRUE]))
    if (isTRUE(exclude_low_q95_power) && "q95_power_class" %in% colnames(q95)) invalid <- union(invalid, as.character(q95$reaction_id[as.character(q95$q95_power_class) == "very_low"]))
  }
  conf <- layer1$reaction_confidence
  if (isTRUE(require_complete_gpr) && is.data.frame(conf) && "reaction_id" %in% colnames(conf)) {
    unsupported_col <- if ("reaction_unsupported_by_complete_gpr_flag" %in% colnames(conf)) "reaction_unsupported_by_complete_gpr_flag" else if ("no_complete_gpr_group_flag" %in% colnames(conf)) "no_complete_gpr_group_flag" else NA_character_
    if (!is.na(unsupported_col)) invalid <- union(invalid, as.character(conf$reaction_id[conf[[unsupported_col]] %in% TRUE]))
  }
  invalid
}
