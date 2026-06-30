#' Select reactions for Layer 3 demand QP
#'
#' Selected reactions are the union of high-variance Layer 1 reactions,
#' exchange reactions, transport reactions, and user-specified reactions. The
#' high-variance term is a practical pool-level differential proxy: reactions
#' with larger variance in `C_rel` vary more across annotated pools and are
#' prioritized for exact feasibility testing.
#'
#' @param C_rel Reaction-by-pool Layer 1 relative capacity matrix.
#' @param reaction_meta Data frame with `reaction_id`, `is_exchange`, and
#' `is_transport` columns.
#' @param top_n Number of top variable Layer 1 reactions to include.
#' @param include_exchange Include reactions marked as exchange.
#' @param include_transport Include reactions marked as transport.
#' @param user_reactions Optional character vector of user-specified reactions.
#'
#' @return Character vector of unique selected reaction IDs.
#' @export
rc_select_reactions <- function(C_rel,
                                reaction_meta,
                                top_n = 500,
                                include_exchange = TRUE,
                                include_transport = TRUE,
                                user_reactions = NULL) {
  if (is.null(rownames(C_rel))) stop("`C_rel` must have reaction IDs as row names.", call. = FALSE)
  if (!is.data.frame(reaction_meta)) stop("`reaction_meta` must be a data.frame.", call. = FALSE)
  required <- c("reaction_id", "is_exchange", "is_transport")
  missing <- setdiff(required, colnames(reaction_meta))
  if (length(missing) > 0L) stop("`reaction_meta` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 0) {
    stop("`top_n` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.logical(include_exchange) || length(include_exchange) != 1L || is.na(include_exchange)) {
    stop("`include_exchange` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(include_transport) || length(include_transport) != 1L || is.na(include_transport)) {
    stop("`include_transport` must be TRUE or FALSE.", call. = FALSE)
  }

  C_mat <- as.matrix(C_rel)
  storage.mode(C_mat) <- "numeric"
  var_score <- matrixStats::rowVars(C_mat, na.rm = TRUE)
  names(var_score) <- rownames(C_mat)
  top_n <- as.integer(top_n)
  top <- names(sort(var_score, decreasing = TRUE))[seq_len(min(top_n, length(var_score)))]

  selected <- unique(c(
    top,
    if (include_exchange) as.character(reaction_meta$reaction_id[as.logical(reaction_meta$is_exchange)]),
    if (include_transport) as.character(reaction_meta$reaction_id[as.logical(reaction_meta$is_transport)]),
    as.character(user_reactions)
  ))
  selected[!is.na(selected) & nzchar(selected)]
}

#' Estimate selected demand QP workload and execution plan
#'
#' @param n_pools Number of pools to solve.
#' @param selected_reactions Character vector of selected reactions.
#' @param seconds_per_qp Expected median seconds per QP solve. Used only for
#' wall-time estimation.
#' @param workers Number of parallel workers planned for Linux execution.
#' @param checkpoint_every Checkpoint interval in solved QPs.
#'
#' @return A one-row data.frame with QP count, time estimate, parallel plan, and
#' checkpoint recommendation.
#' @export
rc_estimate_selected_demand_qp <- function(n_pools,
                                           selected_reactions,
                                           seconds_per_qp = NA_real_,
                                           workers = 1L,
                                           checkpoint_every = 100L) {
  if (!is.numeric(n_pools) || length(n_pools) != 1L || is.na(n_pools) || n_pools < 0) {
    stop("`n_pools` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) || workers < 1) {
    stop("`workers` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(checkpoint_every) || length(checkpoint_every) != 1L || is.na(checkpoint_every) || checkpoint_every < 1) {
    stop("`checkpoint_every` must be a positive integer.", call. = FALSE)
  }
  n_pools <- as.integer(n_pools)
  workers <- as.integer(workers)
  checkpoint_every <- as.integer(checkpoint_every)
  n_selected <- length(unique(as.character(selected_reactions)))
  n_qp <- n_pools * (1L + n_selected)
  estimated_seconds_serial <- if (is.na(seconds_per_qp)) NA_real_ else n_qp * seconds_per_qp
  data.frame(
    n_pools = n_pools,
    n_selected_reactions = n_selected,
    estimated_QP_count = n_qp,
    seconds_per_qp = seconds_per_qp,
    estimated_seconds_serial = estimated_seconds_serial,
    estimated_seconds_parallel = if (is.na(estimated_seconds_serial)) NA_real_ else estimated_seconds_serial / workers,
    workers = workers,
    parallel_plan = if (workers > 1L) "parallelize pools and/or selected reactions with BiocParallel on Linux" else "serial execution",
    checkpoint_every = checkpoint_every,
    expected_checkpoints = ceiling(n_qp / checkpoint_every),
    stringsAsFactors = FALSE
  )
}
