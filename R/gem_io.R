#' Validate a minimal GEM list for RegCompassR QP workflows
#'
#' @param model A list with `S`, `lb`, `ub`, and optional `reaction_id`.
#'
#' @return The validated model with sparse `S` and reaction IDs.
#' @export
rc_validate_gem <- function(model) {
  if (!is.list(model)) {
    stop("`model` must be a list.", call. = FALSE)
  }
  missing <- setdiff(c("S", "lb", "ub"), names(model))
  if (length(missing) > 0) {
    stop("GEM model is missing required fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  S <- Matrix::Matrix(model$S, sparse = TRUE)
  n_rxn <- ncol(S)
  if (length(model$lb) != n_rxn || length(model$ub) != n_rxn) {
    stop("`lb` and `ub` must have length equal to ncol(S).", call. = FALSE)
  }
  if (any(model$lb > model$ub)) {
    stop("All lower bounds must be <= upper bounds.", call. = FALSE)
  }

  reaction_id <- model$reaction_id
  if (is.null(reaction_id)) {
    reaction_id <- colnames(S)
  }
  if (is.null(reaction_id)) {
    reaction_id <- paste0("R", seq_len(n_rxn))
  }
  if (length(reaction_id) != n_rxn || anyNA(reaction_id) || any(!nzchar(reaction_id))) {
    stop("`reaction_id` must contain one non-empty ID per reaction.", call. = FALSE)
  }

  metabolite_id <- model$metabolite_id
  if (is.null(metabolite_id)) metabolite_id <- rownames(S)
  if (is.null(metabolite_id)) metabolite_id <- paste0("M", seq_len(nrow(S)))

  colnames(S) <- reaction_id
  rownames(S) <- metabolite_id
  list(
    S = S,
    lb = as.numeric(model$lb),
    ub = as.numeric(model$ub),
    reaction_id = as.character(reaction_id),
    metabolite_id = as.character(metabolite_id)
  )
}

#' Create a small toy GEM for v0.5 QP tests and examples
#'
#' The toy model has uptake, ATP maintenance, biomass, and lactate demand-like
#' reactions. It is intentionally tiny and is not a biological Human-GEM model.
#'
#' @return A validated GEM list.
#' @export
rc_toy_gem <- function() {
  S <- matrix(
    c(
      1, -1, -1,  0,
      0,  1,  0, -1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("glc_c", "atp_c"),
      c("EX_glc", "ATPM", "BIOMASS", "DM_lac")
    )
  )
  rc_validate_gem(list(
    S = S,
    lb = c(0, 0, 0, 0),
    ub = c(10, 100, 100, 100),
    reaction_id = colnames(S),
    metabolite_id = rownames(S)
  ))
}
