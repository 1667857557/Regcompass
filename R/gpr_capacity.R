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
  med <- matrixStats::rowMedians(X, na.rm = TRUE)
  mad_sigma <- matrixStats::rowMads(X, constant = 1.4826, na.rm = TRUE)
  iqr_sigma <- matrixStats::rowIQRs(X, na.rm = TRUE) / 1.349
  sc <- pmax(mad_sigma, iqr_sigma, min_scale, na.rm = TRUE)
  z <- sweep(X, 1, med, "-")
  z <- sweep(z, 1, sc, "/")
  pmax(pmin(z, z_clip), -z_clip)
}

#' Robust row-wise z score
rc_robust_z <- function(x, eps = 1e-6) {
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a single positive finite number.", call. = FALSE)
  }
  if (is.null(dim(x))) {
    x <- as.numeric(x)
    center <- stats::median(x, na.rm = TRUE)
    scale <- rc_safe_scale(x, min_scale = max(eps, 0.05))
    return((x - center) / scale)
  }
  rc_gene_zscore(x, min_scale = max(eps, 0.05), z_clip = Inf)
}

#' Logistic sigmoid transform
#' @export
rc_sigmoid <- function(z) 1 / (1 + exp(-z))

#' Gene capacity score from meta-cell-level logCPM
#' @export
rc_gene_score <- function(X, min_scale = 0.05, z_clip = 6) {
  rc_sigmoid(rc_gene_zscore(X, min_scale = min_scale, z_clip = z_clip))
}

#' Boltzmann-weighted minimum-biased average for GPR AND
#'
#' Default tau = 0.20 is a biologically moderate bottleneck model for multi-subunit
#' enzymes. Smaller tau, such as 0.08, behaves closer to a hard minimum; larger tau
#' moves toward the arithmetic mean.
#' @export
rc_boltzmann_minavg <- function(scores, tau = 0.20) {
  if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau <= 0) stop("`tau` must be a single positive number.", call. = FALSE)
  scores <- scores[is.finite(scores)]
  if (length(scores) == 0L) return(NA_real_)
  if (length(scores) == 1L) return(scores)
  z <- -scores / tau
  z <- z - max(z)
  w <- exp(z)
  w <- w / sum(w)
  sum(w * scores)
}


#' AND aggregation for one GPR complex
#'
#' Implements the plan-supported sensitivity choices: hard minimum,
#' Boltzmann-weighted minimum-biased average, and arithmetic mean.
#' @export
rc_and_capacity <- function(scores, method = c("boltzmann", "min", "mean"), tau = 0.20) {
  method <- match.arg(method)
  scores <- scores[is.finite(scores)]
  if (length(scores) == 0L) return(NA_real_)
  switch(
    method,
    min = min(scores),
    mean = mean(scores),
    boltzmann = rc_boltzmann_minavg(scores, tau = tau)
  )
}

#' OR aggregation across isoenzyme groups
#'
#' The default `sum` returns an unbounded isoenzyme-summed reaction capacity
#' potential for within-reaction comparisons. `max` and `prob_or` provide bounded
#' sensitivity diagnostics, while `sum_sqrtK` dampens isoenzyme-rich reactions.
#' @export
rc_or_capacity <- function(and_capacities,
                           method = c("sum_sqrtK", "max", "prob_or", "sum")) {
  method <- match.arg(method)
  x <- and_capacities[is.finite(and_capacities)]
  if (length(x) == 0L) return(NA_real_)
  switch(
    method,
    sum = sum(x),
    max = max(x),
    prob_or = 1 - prod(1 - pmin(pmax(x, 0), 1)),
    sum_sqrtK = sum(x) / sqrt(length(x))
  )
}

#' Compute capacity for one reaction in one pool
#' @export
rc_reaction_capacity_one <- function(parsed_gpr, gene_score_vec, tau = 0.20, and_method = c("boltzmann", "min", "mean"),
                                     or_method = c("sum_sqrtK", "max", "prob_or", "sum")) {
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)
  and_caps <- vapply(parsed_gpr, function(and_group) {
    genes <- unique(tolower(and_group))
    vals <- gene_score_vec[genes]
    if (length(genes) == 0L || anyNA(vals) || any(!is.finite(vals))) {
      return(NA_real_)
    }
    rc_and_capacity(vals, method = and_method, tau = tau)
  }, numeric(1))

  rc_or_capacity(and_caps, method = or_method)
}

#' Compute raw Layer 1 reaction capacity
#'
#' Uses fixed sqrt promiscuity correction, Boltzmann AND with tau = 0.20 by default,
#' and square-root dampened OR-group summation. These defaults are the main biological model; alternative
#' tau values should be interpreted only as sensitivity to the multi-subunit
#' bottleneck assumption.
#' @export
rc_reaction_capacity <- function(gpr_list,
                                 gene_score,
                                 promiscuity_mode = c("sqrt", "linear", "none"),
                                 tau = 0.20,
                                 and_method = c("boltzmann", "min", "mean"),
                                 or_method = c("sum_sqrtK", "max", "prob_or", "sum"),
                                 BPPARAM = NULL) {
  promiscuity_mode <- match.arg(promiscuity_mode)
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)
  if (is.null(names(gpr_list))) {
    stop("`gpr_list` must be named by reaction IDs.", call. = FALSE)
  }
  gene_score <- as.matrix(gene_score)
  if (is.null(rownames(gene_score)) || is.null(colnames(gene_score))) {
    stop("`gene_score` must have rownames as genes and colnames as pools.", call. = FALSE)
  }
  normalized_genes <- tolower(trimws(rownames(gene_score)))
  if (anyNA(normalized_genes) || any(!nzchar(normalized_genes))) {
    stop("`gene_score` gene row names must be non-missing and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(normalized_genes)) {
    stop("`gene_score` contains duplicated gene IDs after case normalization.", call. = FALSE)
  }
  rownames(gene_score) <- normalized_genes
  reaction_ids <- names(gpr_list)
  if (!length(reaction_ids)) {
    return(matrix(
      numeric(), nrow = 0L, ncol = ncol(gene_score),
      dimnames = list(character(), colnames(gene_score))
    ))
  }
  if (anyNA(reaction_ids) || any(!nzchar(trimws(reaction_ids))) || anyDuplicated(reaction_ids)) {
    stop("`gpr_list` reaction names must be unique, non-missing, and non-empty.", call. = FALSE)
  }

  weights <- rc_promiscuity_weight(gpr_list, mode = promiscuity_mode)
  common_genes <- intersect(rownames(gene_score), names(weights))
  weighted_score <- gene_score
  weighted_score[common_genes, ] <- sweep(weighted_score[common_genes, , drop = FALSE], 1, weights[common_genes], "*")

  per_reaction <- rc_internal_lapply(reaction_ids, function(rid) {
    parsed <- gpr_list[[rid]]
    vapply(seq_len(ncol(weighted_score)), function(j) {
      score_vector <- weighted_score[, j, drop = TRUE]
      names(score_vector) <- rownames(weighted_score)
      rc_reaction_capacity_one(parsed, score_vector, tau = tau, and_method = and_method, or_method = or_method)
    }, numeric(1))
  }, BPPARAM = BPPARAM)

  out <- do.call(rbind, per_reaction)
  rownames(out) <- reaction_ids
  colnames(out) <- colnames(gene_score)
  out
}

#' Reaction-level GPR diagnostics
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
      missing_subunit_fraction = if (n_genes == 0L) NA_real_ else n_missing / n_genes,
      missing_subunit_flag = n_missing > 0L,
      capacity_missing_flag = n_genes == 0L || n_missing == n_genes,
      stringsAsFactors = FALSE
    )
  }))
}

rc_hard_min_capacity <- function(gpr_list, gene_score, BPPARAM = NULL) {
  gene_score <- as.matrix(gene_score)
  if (is.null(rownames(gene_score)) || is.null(colnames(gene_score))) {
    stop("`gene_score` must have gene rownames and unit colnames.", call. = FALSE)
  }
  normalized_genes <- tolower(trimws(rownames(gene_score)))
  if (anyDuplicated(normalized_genes)) {
    stop("`gene_score` contains duplicated gene IDs after case normalization.", call. = FALSE)
  }
  rownames(gene_score) <- normalized_genes
  reaction_ids <- names(gpr_list)
  if (is.null(reaction_ids)) stop("`gpr_list` must be named by reaction IDs.", call. = FALSE)
  if (!length(reaction_ids)) {
    return(matrix(numeric(), nrow = 0L, ncol = ncol(gene_score),
                  dimnames = list(character(), colnames(gene_score))))
  }
  weights <- rc_promiscuity_weight(gpr_list)
  common_genes <- intersect(rownames(gene_score), names(weights))
  gene_score[common_genes, ] <- sweep(
    gene_score[common_genes, , drop = FALSE],
    1, weights[common_genes], "*"
  )
  per_reaction <- rc_internal_lapply(reaction_ids, function(rid) {
    parsed <- gpr_list[[rid]]
    vapply(seq_len(ncol(gene_score)), function(j) {
      and_caps <- vapply(parsed, function(and_group) {
        genes <- unique(tolower(trimws(as.character(and_group))))
        genes <- genes[!is.na(genes) & nzchar(genes)]
        vals <- gene_score[genes, j]
        if (!length(genes) || length(vals) != length(genes) ||
            anyNA(vals) || any(!is.finite(vals))) {
          return(NA_real_)
        }
        min(vals)
      }, numeric(1))
      rc_or_capacity(and_caps, method = "sum_sqrtK")
    }, numeric(1))
  }, BPPARAM = BPPARAM)
  out <- do.call(rbind, per_reaction)
  rownames(out) <- reaction_ids
  colnames(out) <- colnames(gene_score)
  out
}
