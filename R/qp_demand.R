#' Add a selected reaction demand constraint to a baseline QP
#'
#' @param qp_base QP list from [rc_build_qp()] or [rc_build_baseline_qp()].
#' @param reaction_id Reaction ID or numeric index to constrain.
#' @param delta Demand flux. Positive values impose `v_reaction >= delta`; negative
#' values impose `v_reaction <= delta`.
#'
#' @return A modified QP list with one extra demand constraint row.
#' @export
rc_demand_qp <- function(qp_base, reaction_id, delta) {
  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta)) {
    stop("`delta` must be a single numeric value.", call. = FALSE)
  }
  idx <- rc_reaction_index(reaction_id, qp_base$reaction_id)
  demand_row <- Matrix::sparseMatrix(i = 1L, j = idx, x = 1, dims = c(1L, qp_base$n_reactions))
  qp <- qp_base
  qp$A <- rbind(qp$A, demand_row)
  if (delta >= 0) {
    qp$l <- c(qp$l, delta)
    qp$u <- c(qp$u, Inf)
  } else {
    qp$l <- c(qp$l, -Inf)
    qp$u <- c(qp$u, delta)
  }
  qp$demand_reaction <- qp_base$reaction_id[[idx]]
  qp$demand_delta <- delta
  qp
}

#' Solve demand QPs for selected reactions
#'
#' @param qp_base Baseline QP list.
#' @param reactions Character or integer vector of selected reactions.
#' @param delta Demand flux applied to each reaction, or one value per reaction.
#' @param settings OSQP settings passed to [rc_solve_qp()].
#' @param BPPARAM Optional `BiocParallelParam` for selected reaction parallelism.
#'
#' @return A data.frame with reaction ID, delta, OSQP status, objective, and flux.
#' @export
rc_solve_selected_demand_qp <- function(qp_base,
                                        reactions,
                                        delta,
                                        settings = list(verbose = FALSE),
                                        BPPARAM = NULL) {
  if (length(delta) == 1L) delta <- rep(delta, length(reactions))
  if (length(delta) != length(reactions)) stop("`delta` must have length 1 or match `reactions`.", call. = FALSE)

  pieces <- rc_parallel_lapply(seq_along(reactions), function(i) {
    qp <- rc_demand_qp(qp_base, reactions[[i]], delta[[i]])
    sol <- rc_solve_qp(qp, settings = settings)
    status <- rc_osqp_status(sol)
    idx <- rc_reaction_index(reactions[[i]], qp_base$reaction_id)
    data.frame(
      reaction_id = qp_base$reaction_id[[idx]],
      delta = delta[[i]],
      osqp_status = status,
      objective = if (!is.null(sol$info$obj_val)) sol$info$obj_val else NA_real_,
      flux = if (!is.null(sol$x)) sol$x[[idx]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }, BPPARAM = BPPARAM)
  do.call(rbind, pieces)
}
