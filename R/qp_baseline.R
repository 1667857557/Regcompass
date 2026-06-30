#' Build a minimal sparse QP for a toy GEM
#'
#' The QP minimizes `0.5 * v' P v + q' v` subject to steady-state mass balance,
#' reaction bounds, and an optional ATP maintenance lower-bound constraint. This
#' v0.5 implementation uses one flux variable per reaction; reversible reactions
#' are represented by negative lower bounds and are not split into plus/minus
#' variables yet.
#'
#' @param S Stoichiometric matrix, metabolites by reactions.
#' @param lb Reaction lower bounds.
#' @param ub Reaction upper bounds.
#' @param penalty Non-negative per-reaction quadratic penalties.
#' @param lambda Small non-negative ridge added to all diagonal penalties.
#' @param atpm_rxn Optional ATP maintenance reaction ID or numeric index.
#' @param atpm_min Minimum ATP maintenance flux when `atpm_rxn` is supplied.
#' @param reaction_id Optional reaction IDs; defaults to `colnames(S)` or `R1...`.
#'
#' @return A list with OSQP matrices `P`, `q`, `A`, `l`, `u` and metadata.
#' @export
rc_build_qp <- function(S, lb, ub, penalty, lambda = 1e-4, atpm_rxn = NULL, atpm_min = 0, reaction_id = NULL) {
  S <- Matrix::Matrix(S, sparse = TRUE)
  n_rxn <- ncol(S)
  if (length(lb) != n_rxn || length(ub) != n_rxn || length(penalty) != n_rxn) {
    stop("`lb`, `ub`, and `penalty` must have length equal to ncol(S).", call. = FALSE)
  }
  if (any(lb > ub)) stop("All lower bounds must be <= upper bounds.", call. = FALSE)
  if (any(penalty < 0 | !is.finite(penalty))) stop("`penalty` must be finite and non-negative.", call. = FALSE)
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda < 0) {
    stop("`lambda` must be a single non-negative number.", call. = FALSE)
  }
  if (is.null(reaction_id)) reaction_id <- colnames(S)
  if (is.null(reaction_id)) reaction_id <- paste0("R", seq_len(n_rxn))
  if (length(reaction_id) != n_rxn) stop("`reaction_id` length must equal ncol(S).", call. = FALSE)

  P <- Matrix::Diagonal(n_rxn, x = as.numeric(penalty) + lambda)
  q <- rep(0, n_rxn)

  A <- rbind(S, Matrix::Diagonal(n_rxn))
  l <- c(rep(0, nrow(S)), as.numeric(lb))
  u <- c(rep(0, nrow(S)), as.numeric(ub))

  atpm_index <- rc_reaction_index(atpm_rxn, reaction_id, allow_null = TRUE)
  if (!is.null(atpm_index)) {
    atpm_row <- Matrix::sparseMatrix(i = 1L, j = atpm_index, x = 1, dims = c(1L, n_rxn))
    A <- rbind(A, atpm_row)
    l <- c(l, atpm_min)
    u <- c(u, Inf)
  }

  list(
    P = Matrix::Matrix(P, sparse = TRUE),
    q = q,
    A = Matrix::Matrix(A, sparse = TRUE),
    l = l,
    u = u,
    reaction_id = as.character(reaction_id),
    n_mass_balance = nrow(S),
    n_reactions = n_rxn,
    atpm_rxn = atpm_rxn,
    atpm_min = atpm_min
  )
}

#' Build a baseline QP directly from a minimal GEM list
#'
#' @param model A GEM list accepted by [rc_validate_gem()].
#' @inheritParams rc_build_qp
#'
#' @return A QP list from [rc_build_qp()].
#' @export
rc_build_baseline_qp <- function(model, penalty = NULL, lambda = 1e-4, atpm_rxn = NULL, atpm_min = 0) {
  model <- rc_validate_gem(model)
  if (is.null(penalty)) penalty <- rep(1, length(model$reaction_id))
  rc_build_qp(model$S, model$lb, model$ub, penalty, lambda = lambda, atpm_rxn = atpm_rxn, atpm_min = atpm_min, reaction_id = model$reaction_id)
}

#' Solve a minimal QP with rosqp/OSQP
#'
#' @param qp A QP list from [rc_build_qp()].
#' @param settings Named list of OSQP settings, for example `list(verbose = FALSE)`.
#'
#' @return Solver result with reaction IDs attached when a primal solution exists.
#' @export
rc_solve_qp <- function(qp, settings = list(verbose = FALSE)) {
  if (!requireNamespace("rosqp", quietly = TRUE)) {
    stop("The `rosqp` package must be installed to solve QPs.", call. = FALSE)
  }
  settings <- utils::modifyList(list(verbose = FALSE), settings)
  pars <- do.call(rosqp::osqpSettings, settings)
  ans <- rosqp::solve_osqp(P = qp$P, q = qp$q, A = qp$A, l = qp$l, u = qp$u, pars = pars)
  if (!is.null(ans$x)) names(ans$x) <- qp$reaction_id
  ans
}

#' Extract OSQP status string from a solver result
#'
#' @param solution Result returned by [rc_solve_qp()].
#'
#' @return A character OSQP status or `NA_character_`.
#' @export
rc_osqp_status <- function(solution) {
  if (!is.null(solution$info$status)) return(as.character(solution$info$status))
  if (!is.null(solution$status)) return(as.character(solution$status))
  NA_character_
}

rc_reaction_index <- function(reaction, reaction_id, allow_null = FALSE) {
  if (is.null(reaction)) {
    if (allow_null) return(NULL)
    stop("`reaction` must not be NULL.", call. = FALSE)
  }
  if (is.numeric(reaction)) {
    idx <- as.integer(reaction)
    if (length(idx) != 1L || is.na(idx) || idx < 1L || idx > length(reaction_id)) {
      stop("Reaction index is out of bounds.", call. = FALSE)
    }
    return(idx)
  }
  idx <- match(as.character(reaction), reaction_id)
  if (is.na(idx)) stop("Reaction ID not found: ", reaction, call. = FALSE)
  idx
}
