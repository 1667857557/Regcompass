#' Validate a GEM-like object for RegCompassR Layer 2
#'
#' `gem` may be a list with `S`, `lb`, and `ub`, or those fields plus optional
#' `reactions`, `metabolites`, `reaction_meta`, and `medium_policy`.
rc_validate_gem <- function(gem, selected_reactions = NULL,
                            allow_zero_support = TRUE) {
  if (!is.list(gem) || is.null(gem$S)) {
    stop("`gem` must be a list containing `S`.", call. = FALSE)
  }
  if (!is.logical(allow_zero_support) || length(allow_zero_support) != 1L ||
      is.na(allow_zero_support)) {
    stop("`allow_zero_support` must be TRUE or FALSE.", call. = FALSE)
  }
  S <- .rc_as_dgCMatrix(gem$S)
  if (is.null(colnames(S))) {
    stop("`gem$S` must have reaction IDs in colnames.", call. = FALSE)
  }
  reaction_ids <- trimws(as.character(colnames(S)))
  if (anyNA(reaction_ids) || any(!nzchar(reaction_ids))) {
    stop("`gem$S` reaction IDs must be non-missing and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(reaction_ids)) {
    stop("`gem$S` has duplicated reaction IDs.", call. = FALSE)
  }
  colnames(S) <- reaction_ids
  if (is.null(rownames(S))) rownames(S) <- paste0("met_", seq_len(nrow(S)))
  metabolite_ids <- trimws(as.character(rownames(S)))
  if (anyNA(metabolite_ids) || any(!nzchar(metabolite_ids))) {
    stop("`gem$S` metabolite IDs must be non-missing and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(metabolite_ids)) {
    stop("`gem$S` has duplicated metabolite IDs.", call. = FALSE)
  }
  rownames(S) <- metabolite_ids
  if (length(S@x) && any(!is.finite(S@x))) {
    stop("`gem$S` contains non-finite stoichiometric coefficients.", call. = FALSE)
  }

  lb <- rc_align_bound(gem$lb, reaction_ids, default = -1000, name = "lb")
  ub <- rc_align_bound(gem$ub, reaction_ids, default = 1000, name = "ub")
  if (any(!is.finite(lb)) || any(!is.finite(ub))) {
    stop("Bounds must be finite numeric values.", call. = FALSE)
  }
  if (any(lb > ub)) {
    bad <- reaction_ids[lb > ub]
    stop(
      "Every lower bound must be <= its upper bound; invalid reactions: ",
      paste(utils::head(bad, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(selected_reactions)) {
    selected_reactions <- unique(trimws(as.character(selected_reactions)))
    selected_reactions <- selected_reactions[
      !is.na(selected_reactions) & nzchar(selected_reactions)
    ]
    missing <- setdiff(selected_reactions, reaction_ids)
    if (length(missing) > 0L) {
      stop(
        "Selected reactions missing from GEM: ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }
  }
  zero_cols <- reaction_ids[Matrix::colSums(abs(S) > 0) == 0]
  if (!allow_zero_support && length(zero_cols)) {
    stop(
      "Zero-stoichiometry reactions are not allowed: ",
      paste(utils::head(zero_cols, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  list(
    S = S, lb = lb, ub = ub,
    reactions = reaction_ids, metabolites = metabolite_ids,
    zero_column_reactions = zero_cols,
    n_reactions = length(reaction_ids), n_metabolites = nrow(S),
    valid = TRUE
  )
}

rc_align_bound <- function(x, rxns, default, name, allow_partial = FALSE) {
  rxns <- as.character(rxns)
  if (!is.logical(allow_partial) || length(allow_partial) != 1L ||
      is.na(allow_partial)) {
    stop("`allow_partial` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(default) || length(default) != 1L || !is.finite(default)) {
    stop("`default` must be one finite numeric value.", call. = FALSE)
  }
  if (is.null(x)) {
    return(stats::setNames(rep(default, length(rxns)), rxns))
  }
  if (is.data.frame(x)) {
    if (!all(c("reaction_id", name) %in% colnames(x))) {
      stop(
        "Bound data frames need `reaction_id` and `", name, "` columns.",
        call. = FALSE
      )
    }
    ids <- trimws(as.character(x$reaction_id))
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("Bound reaction IDs must be unique, non-missing, and non-empty.", call. = FALSE)
    }
    unknown <- setdiff(ids, rxns)
    if (length(unknown)) {
      stop(
        "Unknown reaction IDs in `", name, "`: ",
        paste(utils::head(unknown, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    missing <- setdiff(rxns, ids)
    if (length(missing) && !allow_partial) {
      stop(
        "`", name, "` is missing bounds for reactions: ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    values <- suppressWarnings(as.numeric(x[[name]]))
    if (length(values) != length(ids) || any(!is.finite(values))) {
      stop("`", name, "` values must be finite numeric values.", call. = FALSE)
    }
    out <- stats::setNames(rep(default, length(rxns)), rxns)
    out[ids] <- values
    return(out[rxns])
  }

  original_names <- names(x)
  values <- suppressWarnings(as.numeric(x))
  if (is.null(original_names)) {
    if (length(values) != length(rxns)) {
      stop(
        "Unnamed `", name, "` must align with every GEM reaction.",
        call. = FALSE
      )
    }
    if (any(!is.finite(values))) {
      stop("`", name, "` values must be finite.", call. = FALSE)
    }
    return(stats::setNames(values, rxns))
  }
  ids <- trimws(as.character(original_names))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Named bounds must have unique, non-empty reaction IDs.", call. = FALSE)
  }
  unknown <- setdiff(ids, rxns)
  if (length(unknown)) {
    stop(
      "Unknown reaction IDs in `", name, "`: ",
      paste(utils::head(unknown, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  missing <- setdiff(rxns, ids)
  if (length(missing) && !allow_partial) {
    stop(
      "`", name, "` is missing bounds for reactions: ",
      paste(utils::head(missing, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(!is.finite(values))) {
    stop("`", name, "` values must be finite.", call. = FALSE)
  }
  out <- stats::setNames(rep(default, length(rxns)), rxns)
  out[ids] <- values
  out[rxns]
}
