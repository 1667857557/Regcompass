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
#' @param checkpoint_file Optional path to an RDS checkpoint. Completed results are
#' loaded before solving and saved during execution. Checkpointing is disabled
#' when `BPPARAM` is supplied because parallel workers cannot safely append to the
#' same checkpoint file.
#' @param checkpoint_every Save checkpoint after this many newly solved demand QPs.
#'
#' @return A data.frame with reaction ID, delta, OSQP status, objective, and flux.
#' @export
rc_solve_selected_demand_qp <- function(qp_base,
                                        reactions,
                                        delta,
                                        settings = list(verbose = FALSE),
                                        BPPARAM = NULL,
                                        checkpoint_file = NULL,
                                        checkpoint_every = 100L) {
  if (length(delta) == 1L) delta <- rep(delta, length(reactions))
  if (length(delta) != length(reactions)) stop("`delta` must have length 1 or match `reactions`.", call. = FALSE)
  reaction_ids <- vapply(reactions, function(r) qp_base$reaction_id[[rc_reaction_index(r, qp_base$reaction_id)]], character(1L))
  if (!is.numeric(checkpoint_every) || length(checkpoint_every) != 1L || is.na(checkpoint_every) || checkpoint_every < 1) {
    stop("`checkpoint_every` must be a positive integer.", call. = FALSE)
  }
  checkpoint_every <- as.integer(checkpoint_every)
  if (!is.null(checkpoint_file) && !is.null(BPPARAM)) {
    stop("`checkpoint_file` cannot be used with parallel `BPPARAM`; checkpoint by pool or run serial demand solves.", call. = FALSE)
  }

  make_one <- function(i) {
    qp <- rc_demand_qp(qp_base, reaction_ids[[i]], delta[[i]])
    sol <- rc_solve_qp(qp, settings = settings)
    status <- rc_osqp_status(sol)
    idx <- rc_reaction_index(reaction_ids[[i]], qp_base$reaction_id)
    data.frame(
      reaction_id = qp_base$reaction_id[[idx]],
      delta = delta[[i]],
      osqp_status = status,
      objective = if (!is.null(sol$info$obj_val)) sol$info$obj_val else NA_real_,
      flux = if (!is.null(sol$x)) sol$x[[idx]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  empty_out <- function() data.frame(
    reaction_id = character(), delta = numeric(), osqp_status = character(),
    objective = numeric(), flux = numeric(), stringsAsFactors = FALSE
  )

  if (is.null(checkpoint_file)) {
    pieces <- rc_parallel_lapply(seq_along(reactions), make_one, BPPARAM = BPPARAM)
    return(if (length(pieces) > 0L) do.call(rbind, pieces) else empty_out())
  }

  done <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else data.frame()
  done_keys <- if (nrow(done) > 0L) paste(done$reaction_id, done$delta, sep = "\r") else character()
  target_keys <- paste(reaction_ids, delta, sep = "\r")
  pieces <- if (nrow(done) > 0L) list(done) else list()
  solved_since_checkpoint <- 0L
  for (i in which(!target_keys %in% done_keys)) {
    pieces[[length(pieces) + 1L]] <- make_one(i)
    solved_since_checkpoint <- solved_since_checkpoint + 1L
    if (solved_since_checkpoint >= checkpoint_every) {
      saveRDS(do.call(rbind, pieces), checkpoint_file)
      solved_since_checkpoint <- 0L
    }
  }
  out <- if (length(pieces) > 0L) do.call(rbind, pieces) else empty_out()
  saveRDS(out, checkpoint_file)
  out
}
