#' Validate a GEM-like object for RegCompassR Layer 2
#'
#' `gem` may be a list with `S`, `lb`, and `ub`, or those fields plus optional
#' `reactions`, `metabolites`, `reaction_meta`, and `medium_policy`.
#' @export
rc_validate_gem <- function(gem, selected_reactions = NULL, allow_zero_support = TRUE) {
  if (!is.list(gem) || is.null(gem$S)) stop("`gem` must be a list containing `S`.", call. = FALSE)
  S <- if (inherits(gem$S, "sparseMatrix")) methods::as(gem$S, "dgCMatrix") else methods::as(gem$S, "dgCMatrix")
  if (is.null(colnames(S))) stop("`gem$S` must have reaction IDs in colnames.", call. = FALSE)
  if (is.null(rownames(S))) rownames(S) <- paste0("met_", seq_len(nrow(S)))
  if (anyDuplicated(colnames(S))) stop("`gem$S` has duplicated reaction IDs.", call. = FALSE)
  rxns <- colnames(S)
  lb <- rc_align_bound(gem$lb, rxns, default = -1000, name = "lb")
  ub <- rc_align_bound(gem$ub, rxns, default = 1000, name = "ub")
  if (any(!is.finite(lb)) || any(!is.finite(ub))) stop("Bounds must be finite numeric values.", call. = FALSE)
  if (any(lb > ub)) stop("Every lower bound must be <= its upper bound.", call. = FALSE)
  if (!is.null(selected_reactions)) {
    missing <- setdiff(selected_reactions, rxns)
    if (length(missing) > 0L) stop("Selected reactions missing from GEM: ", paste(utils::head(missing, 10), collapse = ", "), call. = FALSE)
  }
  zero_cols <- rxns[Matrix::colSums(abs(S) > 0) == 0]
  list(S = S, lb = lb, ub = ub, reactions = rxns, metabolites = rownames(S), zero_column_reactions = zero_cols,
       n_reactions = length(rxns), n_metabolites = nrow(S), valid = TRUE)
}

rc_align_bound <- function(x, rxns, default, name) {
  if (is.null(x)) return(stats::setNames(rep(default, length(rxns)), rxns))
  if (is.data.frame(x)) {
    if (!all(c("reaction_id", name) %in% colnames(x))) stop("Bound data frames need `reaction_id` and `", name, "` columns.", call. = FALSE)
    out <- stats::setNames(rep(default, length(rxns)), rxns)
    out[as.character(x$reaction_id)] <- as.numeric(x[[name]])
    return(out[rxns])
  }
  x <- as.numeric(x)
  if (is.null(names(x))) {
    if (length(x) != length(rxns)) stop("Unnamed `", name, "` must align with GEM reactions.", call. = FALSE)
    names(x) <- rxns
  }
  out <- stats::setNames(rep(default, length(rxns)), rxns)
  out[intersect(names(x), rxns)] <- x[intersect(names(x), rxns)]
  out[rxns]
}
