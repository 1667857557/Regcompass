#' Safe robust sigma-scale estimate
#' @export
rc_safe_scale <- function(x, min_scale = 0.05) {
  mad_sigma <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  iqr_sigma <- stats::IQR(x, na.rm = TRUE) / 1.349
  max(mad_sigma, iqr_sigma, min_scale, na.rm = TRUE)
}

#' Robust row-wise z score clipped to a finite range
#' @export
rc_gene_zscore <- function(X, min_scale = 0.05, z_clip = 6) {
  X <- as.matrix(X)
  z <- X
  for (i in seq_len(nrow(X))) {
    x <- as.numeric(X[i, ])
    med <- stats::median(x, na.rm = TRUE)
    sc <- rc_safe_scale(x, min_scale = min_scale)
    zi <- (x - med) / sc
    zi <- pmax(pmin(zi, z_clip), -z_clip)
    z[i, ] <- zi
  }
  z
}

#' Robust row-wise z score
#' @export
rc_robust_z <- function(x, eps = 1e-6) {
  rc_gene_zscore(x, min_scale = max(eps, 0.05), z_clip = Inf)
}

#' Sigmoid gene score from normalized pool expression
#' @export
rc_gene_score <- function(X, min_scale = 0.05, z_clip = 6) rc_sigmoid(rc_gene_zscore(X, min_scale = min_scale, z_clip = z_clip))

#' Logistic sigmoid transform
#'
#' @param z Numeric vector, matrix, or array.
#'
#' @return `1 / (1 + exp(-z))` with the same shape as `z`.
#' @export
rc_sigmoid <- function(z) 1 / (1 + exp(-z))

#' Boltzmann-weighted average biased toward the minimum
#'
#' This is not a LogSumExp soft minimum. It is a Boltzmann-weighted average of
#' finite input scores where lower scores receive larger weights.
#'
#' @param scores Numeric scores.
#' @param tau Positive temperature. Smaller values put more weight on the minimum.
#'
#' @return A single numeric value between `min(scores)` and `mean(scores)` for
#' finite non-empty `scores`; `NA_real_` if no finite scores are supplied.
#' @export
rc_boltzmann_minavg <- function(scores, tau = 0.08) {
  if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau <= 0) {
    stop("`tau` must be a single positive number.", call. = FALSE)
  }
  scores <- scores[is.finite(scores)]
  if (length(scores) == 0L) return(NA_real_)
  if (length(scores) == 1L) return(scores)

  z <- -scores / tau
  z <- z - max(z)
  w <- exp(z)
  w <- w / sum(w)
  sum(w * scores)
}

#' Compute capacity for one reaction in one pool
#'
#' @param parsed_gpr Parsed GPR rule returned by [rc_parse_gpr_simple()].
#' @param gene_score_vec Named numeric gene scores for one pool.
#' @param tau Positive Boltzmann temperature for AND complexes.
#'
#' @return A single raw Layer 1 reaction capacity potential.
#' @export
rc_reaction_capacity_one <- function(parsed_gpr, gene_score_vec, tau = 0.08) {
  and_caps <- vapply(parsed_gpr, function(and_group) {
    vals <- gene_score_vec[and_group]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) return(NA_real_)
    rc_boltzmann_minavg(vals, tau = tau)
  }, numeric(1))

  if (all(is.na(and_caps))) {
    return(NA_real_)
  }
  sum(and_caps, na.rm = TRUE)
}

#' Compute Layer 1 reaction capacity for all reactions and pools
#'
#' @param gpr_list A named list of parsed GPR rules.
#' @param gene_score A numeric gene-by-pool matrix with lower-case gene IDs in
#' `rownames(gene_score)`. Values should be expression-derived scores, commonly
#' `rc_sigmoid(rc_robust_z(pool_expression))`.
#' @param promiscuity_mode One of `"sqrt"`, `"linear"`, or `"none"`.
#' @param tau Positive Boltzmann temperature for AND complexes.
#' @param BPPARAM Optional `BiocParallelParam` used to parallelize over reactions.
#'
#' @return A reaction-by-pool numeric matrix of raw Layer 1 capacity potentials.
#' @export
rc_reaction_capacity <- function(gpr_list,
                                 gene_score,
                                 promiscuity_mode = c("sqrt", "linear", "none"),
                                 tau = 0.08,
                                 BPPARAM = NULL) {
  promiscuity_mode <- match.arg(promiscuity_mode)
  if (is.null(names(gpr_list))) {
    stop("`gpr_list` must be named by reaction IDs.", call. = FALSE)
  }
  gene_score <- as.matrix(gene_score)
  if (is.null(rownames(gene_score)) || is.null(colnames(gene_score))) {
    stop("`gene_score` must have rownames as genes and colnames as pools.", call. = FALSE)
  }
  rownames(gene_score) <- tolower(rownames(gene_score))

  weights <- rc_promiscuity_weight(gpr_list, mode = promiscuity_mode)
  weighted_score <- gene_score
  common_genes <- intersect(rownames(weighted_score), names(weights))
  weighted_score[common_genes, ] <- sweep(weighted_score[common_genes, , drop = FALSE], 1, weights[common_genes], "*")

  reaction_ids <- names(gpr_list)
  per_reaction <- rc_parallel_lapply(reaction_ids, function(rid) {
    parsed <- gpr_list[[rid]]
    vapply(seq_len(ncol(weighted_score)), function(j) {
      rc_reaction_capacity_one(parsed, weighted_score[, j], tau = tau)
    }, numeric(1))
  }, BPPARAM = BPPARAM)

  out <- do.call(rbind, per_reaction)
  rownames(out) <- reaction_ids
  colnames(out) <- colnames(gene_score)
  out
}

#' Build GPR diagnostics for reactions
#'
#' @param gpr_list A named list of parsed GPR rules.
#' @param gene_ids Genes available in the scored matrix.
#'
#' @return A data.frame with one row per reaction.
#' @export
rc_gpr_diagnostics <- function(gpr_list, gene_ids) {
  gene_ids <- unique(tolower(gene_ids))
  do.call(rbind, lapply(names(gpr_list), function(rid) {
    rule <- gpr_list[[rid]]
    genes <- unique(unlist(rule, use.names = FALSE))
    genes <- genes[nzchar(genes)]
    n_genes <- length(genes)
    n_missing <- sum(!genes %in% gene_ids)
    data.frame(
      reaction_id = rid,
      n_gpr_genes = n_genes,
      n_and_groups = length(rule),
      has_isoenzyme = length(rule) > 1L,
      has_multisubunit = any(vapply(rule, length, integer(1)) > 1L),
      missing_gene_fraction = if (n_genes == 0L) NA_real_ else n_missing / n_genes,
      capacity_missing_flag = n_genes == 0L || n_missing == n_genes,
      stringsAsFactors = FALSE
    )
  }))
}
