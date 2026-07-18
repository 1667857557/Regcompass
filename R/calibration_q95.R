.rc_q95_weighted_quantile <- function(x, weights, probs = 0.95) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- x[keep]
  weights <- weights[keep]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  order_index <- order(x)
  x <- x[order_index]
  weights <- weights[order_index]
  cumulative <- cumsum(weights) / sum(weights)
  vapply(probs, function(probability) {
    if (!is.finite(probability) || probability < 0 || probability > 1) {
      stop("`probs` values must be between zero and one.", call. = FALSE)
    }
    if (probability <= 0) return(x[[1L]])
    if (probability >= 1) return(x[[length(x)]])
    x[[which(cumulative >= probability)[[1L]]]]
  }, numeric(1))
}

.rc_q95_validate_weights <- function(weights, columns) {
  if (is.null(weights)) {
    output <- rep(1, length(columns))
    names(output) <- columns
    return(output / sum(output))
  }
  if (!is.numeric(weights) || length(weights) != length(columns) ||
      any(!is.finite(weights)) || any(weights <= 0)) {
    stop("`weights` must contain one positive finite value per capacity column.",
         call. = FALSE)
  }
  if (!is.null(names(weights))) {
    missing <- setdiff(columns, names(weights))
    if (length(missing)) {
      stop("Named `weights` are missing capacity columns: ",
           paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
    }
    weights <- weights[columns]
  }
  output <- as.numeric(weights)
  output <- output / sum(output)
  names(output) <- columns
  output
}

.rc_q95_validate_balance_ids <- function(balance_ids, columns) {
  if (is.null(balance_ids)) return(NULL)
  if (length(balance_ids) != length(columns)) {
    stop("`balance_ids` must contain one biological-sample ID per capacity column.",
         call. = FALSE)
  }
  if (!is.null(names(balance_ids))) {
    missing <- setdiff(columns, names(balance_ids))
    if (length(missing)) {
      stop("Named `balance_ids` are missing capacity columns: ",
           paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
    }
    balance_ids <- balance_ids[columns]
  }
  output <- trimws(as.character(balance_ids))
  if (anyNA(output) || any(!nzchar(output))) {
    stop("`balance_ids` must be non-missing and non-empty.", call. = FALSE)
  }
  names(output) <- columns
  output
}

.rc_q95_local_weights <- function(global_weights, balance_ids, columns) {
  if (!is.null(balance_ids)) {
    local <- .rc_equal_sample_weights(balance_ids[columns])
    names(local) <- columns
    return(local)
  }
  local <- global_weights[columns]
  local / sum(local)
}

#' Shrink within-reaction Q95 diagnostics across strata
#'
#' This calibration is diagnostic only. It never replaces the bounded absolute
#' reaction support used by the LP.
rc_q95_shrink <- function(
    C_raw, unit_meta = NULL, stratum_col = NULL,
    q = 0.95, n0 = 80, eps = 1e-6, weights = NULL,
    balance_ids = NULL) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) {
    stop("`C_raw` must have reaction rownames and unit colnames.", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1) {
    stop("`q` must be one number strictly between zero and one.", call. = FALSE)
  }
  if (!is.numeric(n0) || length(n0) != 1L || !is.finite(n0) || n0 < 0) {
    stop("`n0` must be one finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be one positive finite number.", call. = FALSE)
  }
  if (!is.null(weights) && !is.null(balance_ids)) {
    stop("Specify only one of `weights` and `balance_ids`.", call. = FALSE)
  }
  balance_ids <- .rc_q95_validate_balance_ids(
    balance_ids,
    colnames(C_raw)
  )
  weights <- if (!is.null(balance_ids)) {
    .rc_equal_sample_weights(balance_ids)
  } else {
    .rc_q95_validate_weights(weights, colnames(C_raw))
  }
  names(weights) <- colnames(C_raw)
  weighted <- length(unique(round(as.numeric(weights), 14))) > 1L
  sample_balanced <- !is.null(balance_ids) || weighted
  balance_scope <- if (!is.null(balance_ids)) {
    "equal_sample_global_and_within_stratum"
  } else if (weighted) {
    "caller_supplied_global_weights"
  } else {
    "equal_unit"
  }
  quantile_one <- function(values, selected_weights) {
    .rc_q95_weighted_quantile(values, selected_weights, probs = q)
  }
  global_q <- vapply(seq_len(nrow(C_raw)), function(i) {
    quantile_one(C_raw[i, ], weights)
  }, numeric(1))
  global_n <- rowSums(is.finite(C_raw))
  global_n_effective <- vapply(seq_len(nrow(C_raw)), function(i) {
    observed <- is.finite(C_raw[i, ])
    selected <- weights[observed]
    if (!length(selected)) return(0)
    sum(selected)^2 / sum(selected^2)
  }, numeric(1))

  if (!is.null(stratum_col)) {
    if (is.null(unit_meta) || !is.data.frame(unit_meta) ||
        !"pool_id" %in% colnames(unit_meta) ||
        !stratum_col %in% colnames(unit_meta)) {
      stop("`unit_meta` with `pool_id` and `stratum_col` is required.",
           call. = FALSE)
    }
    unit_meta <- unit_meta[
      match(colnames(C_raw), as.character(unit_meta$pool_id)),
      ,
      drop = FALSE
    ]
    if (anyNA(unit_meta$pool_id)) {
      stop("`unit_meta` is missing metadata for some capacity columns.",
           call. = FALSE)
    }
    strata <- trimws(as.character(unit_meta[[stratum_col]]))
    if (anyNA(strata) || any(!nzchar(strata))) {
      stop("Q95 strata must be non-missing and non-empty.", call. = FALSE)
    }
  } else {
    strata <- rep("global", ncol(C_raw))
  }

  relative <- matrix(
    NA_real_,
    nrow = nrow(C_raw),
    ncol = ncol(C_raw),
    dimnames = dimnames(C_raw)
  )
  diagnostics <- list()
  index <- 1L
  for (stratum in unique(strata)) {
    columns <- which(strata == stratum)
    stratum_weights <- .rc_q95_local_weights(
      weights,
      balance_ids,
      columns
    )
    n <- rowSums(is.finite(C_raw[, columns, drop = FALSE]))
    n_effective <- vapply(seq_len(nrow(C_raw)), function(i) {
      observed <- is.finite(C_raw[i, columns])
      selected <- stratum_weights[observed]
      if (!length(selected)) return(0)
      sum(selected)^2 / sum(selected^2)
    }, numeric(1))
    rho <- n_effective / (n_effective + n0)
    q_stratum <- vapply(seq_len(nrow(C_raw)), function(i) {
      quantile_one(C_raw[i, columns], stratum_weights)
    }, numeric(1))
    q_used <- ifelse(is.finite(q_stratum), q_stratum, global_q)
    q_shrink <- rho * q_used + (1 - rho) * global_q
    relative[, columns] <- sweep(
      C_raw[, columns, drop = FALSE],
      1L,
      q_shrink + eps,
      "/"
    )
    diagnostics[[index]] <- data.frame(
      reaction_id = rownames(C_raw),
      stratum = stratum,
      n = as.integer(n),
      n_effective = as.numeric(n_effective),
      n_global = as.integer(global_n),
      n_global_effective = as.numeric(global_n_effective),
      q_stratum = as.numeric(q_stratum),
      q_stratum_used = as.numeric(q_used),
      q_global = as.numeric(global_q),
      rho_n = as.numeric(rho),
      q_value = as.numeric(q_shrink),
      sample_balanced = sample_balanced,
      balance_scope = balance_scope,
      n_balancing_samples = if (is.null(balance_ids)) {
        NA_integer_
      } else {
        length(unique(balance_ids[columns]))
      },
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
  relative <- pmin(pmax(relative, 0), 1)
  finite_range <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    diff(range(x))
  })
  noninformative <- !is.finite(finite_range) | finite_range <= eps
  if (any(noninformative)) relative[noninformative, ] <- NA_real_
  Q <- do.call(rbind, diagnostics)
  Q$all_missing_reaction_flag <- Q$n_global == 0L
  Q$stratum_missing_reaction_flag <- Q$n == 0L
  # Compatibility aliases retained for downstream diagnostics and older tests.
  Q$q_shrink <- Q$q_value
  Q$q95_power_class <- factor(
    ifelse(
      Q$n < 5L,
      "very_low",
      ifelse(
        Q$n < 20L,
        "low",
        ifelse(Q$n < 100L, "moderate",
               ifelse(Q$n < 400L, "adequate", "high"))
      )
    ),
    levels = c("very_low", "low", "moderate", "adequate", "high"),
    ordered = TRUE
  )
  list(
    C_rel = relative,
    Q = Q,
    weights = weights,
    balance_ids = balance_ids
  )
}

.rc_q95_bootstrap_one <- function(x, weights, B = 500, q = 0.95) {
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- as.numeric(x[keep])
  weights <- as.numeric(weights[keep])
  if (length(x) < 20L) {
    return(c(q95 = NA_real_, ci_low = NA_real_, ci_high = NA_real_, width = NA_real_))
  }
  probabilities <- weights / sum(weights)
  estimates <- replicate(B, {
    sampled <- sample.int(length(x), size = length(x), replace = TRUE,
                          prob = probabilities)
    .rc_q95_weighted_quantile(
      x[sampled],
      rep(1, length(sampled)),
      probs = q
    )
  })
  interval <- stats::quantile(
    estimates,
    c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  c(
    q95 = .rc_q95_weighted_quantile(x, weights, probs = q),
    ci_low = interval[[1L]],
    ci_high = interval[[2L]],
    width = diff(interval)
  )
}

#' Calibrate reaction-capacity diagnostics
#'
#' `C_abs` and the compatibility field `C_rel` contain zero-preserving bounded
#' absolute support used by downstream penalties. Q95-normalized values are
#' returned only in `C_within_reaction_relative`.
#'
#' @param C_raw Reaction-by-unit raw bounded support matrix.
#' @param eps Positive numerical tolerance.
#' @param bootstrap Compute weighted bootstrap Q95 intervals.
#' @param B Number of bootstrap replicates.
#' @param BPPARAM Reserved parallel backend parameter.
#' @param n0 Empirical-Bayes shrinkage strength toward the global reaction Q95.
#' @param unit_meta Unit metadata with a `pool_id` column.
#' @param stratum_col Optional diagnostic stratum column.
#' @param weights Optional positive global unit weights. Retained for callers
#'   that provide a custom weighting estimand.
#' @param balance_ids Optional biological-sample IDs aligned to capacity
#'   columns. When supplied, samples receive equal total mass globally and the
#'   weights are recomputed inside every diagnostic stratum so each represented
#'   sample again receives equal total mass.
#' @return A list containing absolute LP support, within-reaction diagnostics,
#'   Q95 diagnostics, and normalized weights.
rc_q95_calibrate <- function(
    C_raw, eps = 1e-6, bootstrap = TRUE,
    B = 500, BPPARAM = NULL, n0 = 80,
    unit_meta = NULL, stratum_col = NULL, weights = NULL,
    balance_ids = NULL) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) {
    stop("`C_raw` must have reaction rownames and unit colnames.", call. = FALSE)
  }
  if (!is.logical(bootstrap) || length(bootstrap) != 1L || is.na(bootstrap)) {
    stop("`bootstrap` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(bootstrap) &&
      (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1)) {
    stop("`B` must be one positive finite number.", call. = FALSE)
  }
  diagnostic <- rc_q95_shrink(
    C_raw,
    unit_meta = unit_meta,
    stratum_col = stratum_col,
    q = 0.95,
    n0 = n0,
    eps = eps,
    weights = weights,
    balance_ids = balance_ids
  )
  absolute <- pmin(pmax(C_raw, 0), 1)
  all_missing <- rowSums(is.finite(C_raw)) == 0L
  if (any(all_missing)) absolute[all_missing, ] <- NA_real_
  all_zero <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) > 0L && max(x) <= eps
  })
  if (any(all_zero)) absolute[all_zero, ] <- 0

  Q <- diagnostic$Q
  reaction_index <- match(Q$reaction_id, rownames(C_raw))
  minimum <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  })
  maximum <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  })
  informative <- rowSums(is.finite(diagnostic$C_rel)) > 0L
  Q$quantile_used <- 0.95
  Q$n_finite <- Q$n
  Q$n_finite_global <- Q$n_global
  Q$low_n_flag <- Q$n_effective < 20
  Q$all_zero_reaction_flag <-
    is.finite(maximum[reaction_index]) & maximum[reaction_index] <= eps
  Q$constant_reaction_flag <-
    is.finite(minimum[reaction_index]) &
    is.finite(maximum[reaction_index]) &
    abs(maximum[reaction_index] - minimum[reaction_index]) <= eps
  Q$raw_out_of_unit_interval_flag <- vapply(
    reaction_index,
    function(i) any(is.finite(C_raw[i, ]) &
                      (C_raw[i, ] < -eps | C_raw[i, ] > 1 + eps)),
    logical(1)
  )
  Q$relative_capacity_informative <- informative[reaction_index]
  Q$calibration_role <- "diagnostic_only_not_lp_capacity"

  if (isTRUE(bootstrap)) {
    normalized_weights <- diagnostic$weights
    bootstrap_rows <- lapply(seq_len(nrow(Q)), function(i) {
      reaction <- Q$reaction_id[[i]]
      stratum <- Q$stratum[[i]]
      columns <- seq_len(ncol(C_raw))
      if (!is.null(stratum_col)) {
        aligned_meta <- unit_meta[
          match(colnames(C_raw), as.character(unit_meta$pool_id)),
          ,
          drop = FALSE
        ]
        columns <- which(as.character(aligned_meta[[stratum_col]]) == stratum)
      }
      bootstrap_weights <- .rc_q95_local_weights(
        normalized_weights,
        diagnostic$balance_ids,
        columns
      )
      .rc_q95_bootstrap_one(
        C_raw[reaction, columns],
        bootstrap_weights,
        B = as.integer(B),
        q = 0.95
      )
    })
    boot <- do.call(rbind, bootstrap_rows)
    Q$q95_bootstrap <- boot[, "q95"]
    Q$q95_ci_low <- boot[, "ci_low"]
    Q$q95_ci_high <- boot[, "ci_high"]
    Q$q95_ci_width <- boot[, "width"]
    Q$q95_unstable_flag <- is.finite(Q$q95_ci_width) &
      (Q$q95_ci_width / pmax(Q$q_value, eps)) > 0.5
  }
  list(
    C_rel = absolute,
    C_abs = absolute,
    C_within_reaction_relative = diagnostic$C_rel,
    Q = Q,
    weights = diagnostic$weights,
    balance_ids = diagnostic$balance_ids,
    BPPARAM = BPPARAM
  )
}

rc_q95_bootstrap <- function(x, B = 500) {
  .rc_q95_bootstrap_one(x, rep(1, length(x)), B = B, q = 0.95)
}

rc_empty_gpr_evidence_matrix <- function(gpr_list, unit_ids) {
  genes <- unique(tolower(unlist(gpr_list, use.names = FALSE)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) genes <- character(0)
  matrix(NA_real_, nrow = length(genes), ncol = length(unit_ids), dimnames = list(genes, unit_ids))
}

rc_set_confidence_source <- function(confidence, source, detection_available = NULL) {
  confidence$confidence_source <- source
  confidence$evidence_source <- source
  if (!is.null(detection_available)) confidence$detection_available <- detection_available
  confidence
}

rc_softmin_conf <- function(v, tau = 0.20) {
  v <- as.numeric(v)
  if (any(!is.finite(v))) return(NA_real_)
  w <- exp(-v / tau)
  sum(w * v) / sum(w)
}

rc_apply_low_confidence_quantile <- function(df, low_confidence_quantile) {
  df$low_confidence_reaction_flag <- NA
  if (is.null(low_confidence_quantile)) return(df)
  if (!is.numeric(low_confidence_quantile) || length(low_confidence_quantile) != 1L || is.na(low_confidence_quantile) || low_confidence_quantile <= 0 || low_confidence_quantile >= 1) {
    stop("`low_confidence_quantile` must be NULL or a single number between 0 and 1.", call. = FALSE)
  }
  split_src <- split(seq_len(nrow(df)), df$confidence_source)
  for (idx in split_src) {
    x <- df$reaction_confidence[idx]
    if (!any(is.finite(x))) next
    threshold <- stats::quantile(x, probs = low_confidence_quantile, na.rm = TRUE, names = FALSE)
    df$low_confidence_reaction_flag[idx] <- is.finite(x) & x < threshold
  }
  df
}

rc_reaction_confidence_gpr_aware <- function(
    gpr_list,
    gene_confidence,
    unit_detection = NULL,
    tau_conf = 0.20,
    and_method = c("softmin", "min", "mean"),
    or_method = c("max", "prob_or", "sum_sqrtK"),
    missing_group_policy = c("complete_group", "partial_group"),
    low_confidence_quantile = NULL,
    low_confidence_threshold = NULL) {
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)
  missing_group_policy <- match.arg(missing_group_policy)
  if (!is.null(low_confidence_threshold)) {
    if (!is.numeric(low_confidence_threshold) ||
        length(low_confidence_threshold) != 1L ||
        !is.finite(low_confidence_threshold) ||
        low_confidence_threshold < 0 ||
        low_confidence_threshold > 1) {
      stop(
        "`low_confidence_threshold` must be NULL or one number between 0 and 1.",
        call. = FALSE
      )
    }
    if (!is.null(low_confidence_quantile)) {
      stop(
        "Specify only one of `low_confidence_threshold` and `low_confidence_quantile`.",
        call. = FALSE
      )
    }
  }
  gene_confidence <- as.matrix(gene_confidence)
  if (is.null(rownames(gene_confidence)) ||
      is.null(colnames(gene_confidence))) {
    stop(
      "`gene_confidence` must have gene rownames and unit colnames.",
      call. = FALSE
    )
  }
  normalized <- tolower(trimws(rownames(gene_confidence)))
  if (anyNA(normalized) || any(!nzchar(normalized)) ||
      anyDuplicated(normalized)) {
    stop(
      "`gene_confidence` must have unique non-empty genes after normalization.",
      call. = FALSE
    )
  }
  rownames(gene_confidence) <- normalized
  unit_ids <- colnames(gene_confidence)

  rows <- lapply(names(gpr_list), function(rid) {
    groups <- gpr_list[[rid]]
    group_data <- lapply(groups, function(group) {
      genes <- unique(tolower(trimws(as.character(group))))
      genes <- genes[!is.na(genes) & nzchar(genes)]
      hit <- intersect(genes, rownames(gene_confidence))
      observed <- if (!length(genes)) NA_real_ else length(hit) / length(genes)
      complete <- length(genes) > 0L && length(hit) == length(genes)
      eligible <- complete ||
        (identical(missing_group_policy, "partial_group") && length(hit) > 0L)
      if (!eligible) {
        return(list(
          values = rep(NA_real_, length(unit_ids)),
          complete = complete,
          eligible = FALSE,
          observed = observed
        ))
      }
      matrix_in <- gene_confidence[hit, , drop = FALSE]
      values <- if (nrow(matrix_in) == 1L) {
        as.numeric(matrix_in[1L, ])
      } else if (identical(and_method, "min")) {
        apply(matrix_in, 2L, function(x) {
          if (any(!is.finite(x))) NA_real_ else min(x)
        })
      } else if (identical(and_method, "mean")) {
        colMeans(matrix_in, na.rm = FALSE)
      } else {
        apply(matrix_in, 2L, rc_softmin_conf, tau = tau_conf)
      }
      list(
        values = values,
        complete = complete,
        eligible = TRUE,
        observed = observed
      )
    })
    if (!length(group_data)) {
      group_data <- list(list(
        values = rep(NA_real_, length(unit_ids)),
        complete = FALSE,
        eligible = FALSE,
        observed = NA_real_
      ))
    }
    group_matrix <- do.call(rbind, lapply(group_data, `[[`, "values"))
    complete <- vapply(group_data, `[[`, logical(1), "complete")
    eligible <- vapply(group_data, `[[`, logical(1), "eligible")
    observed <- vapply(group_data, `[[`, numeric(1), "observed")
    values <- rep(NA_real_, length(unit_ids))
    if (any(eligible)) {
      selected <- group_matrix[eligible, , drop = FALSE]
      values <- if (identical(or_method, "max")) {
        apply(selected, 2L, function(x) {
          x <- x[is.finite(x)]
          if (!length(x)) NA_real_ else max(x)
        })
      } else if (identical(or_method, "prob_or")) {
        apply(selected, 2L, function(x) {
          x <- x[is.finite(x)]
          if (!length(x)) return(NA_real_)
          1 - prod(1 - pmin(pmax(x, 0), 1))
        })
      } else {
        apply(selected, 2L, function(x) {
          x <- x[is.finite(x)]
          if (!length(x)) return(NA_real_)
          sum(x) / sqrt(length(x))
        })
      }
    }
    all_genes <- unique(tolower(unlist(groups, use.names = FALSE)))
    all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
    observed_genes <- intersect(all_genes, rownames(gene_confidence))
    coverage <- if (!length(all_genes)) {
      NA_real_
    } else {
      length(observed_genes) / length(all_genes)
    }
    data.frame(
      reaction_id = rid,
      pool_id = unit_ids,
      reaction_confidence = values,
      reaction_evidence_score = values,
      reaction_confidence_unpenalized = values,
      reaction_evidence_score_unpenalized = values,
      confidence_source = "gpr_aware_gene_confidence",
      evidence_source = "gpr_aware_gene_confidence",
      n_gpr_genes_total = length(all_genes),
      n_gpr_genes_multiome = length(observed_genes),
      multiome_coverage_fraction = coverage,
      missing_gpr_gene_fraction = if (is.na(coverage)) NA_real_ else 1 - coverage,
      missing_gene_fraction = if (is.na(coverage)) NA_real_ else 1 - coverage,
      missing_subunit_confidence_penalty = NA_real_,
      detection_available = !is.null(unit_detection),
      mean_gpr_detection_rate = NA_real_,
      low_confidence_threshold = if (is.null(low_confidence_threshold)) {
        NA_real_
      } else {
        low_confidence_threshold
      },
      n_and_groups_total = length(groups),
      n_and_groups_complete = sum(complete),
      n_and_groups_eligible = sum(eligible),
      complete_and_group_fraction = if (!length(groups)) NA_real_ else mean(complete),
      best_and_group_observed_fraction = if (all(is.na(observed))) {
        NA_real_
      } else {
        max(observed, na.rm = TRUE)
      },
      any_incomplete_gpr_group_flag = any(!complete),
      reaction_unsupported_by_complete_gpr_flag = !any(complete),
      missing_required_subunit_flag = any(!complete),
      no_complete_gpr_group_flag = !any(complete),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (!is.null(low_confidence_threshold)) {
    out$low_confidence_reaction_flag <-
      is.finite(out$reaction_confidence) &
      out$reaction_confidence < low_confidence_threshold
    return(out)
  }
  rc_apply_low_confidence_quantile(out, low_confidence_quantile)
}

rc_gpr_confidence_one <- function(parsed_gpr, evidence_vec) {
  and_vals <- vapply(parsed_gpr, function(and_group) {
    genes <- unique(tolower(and_group))
    vals <- evidence_vec[genes]
    vals <- vals[is.finite(vals)]
    if (length(genes) == 0L || length(vals) == 0L) {
      return(NA_real_)
    }
    min(vals)
  }, numeric(1))
  finite <- and_vals[is.finite(and_vals)]
  if (length(finite) == 0L) return(NA_real_)
  max(finite)
}

rc_gpr_confidence_matrix <- function(gpr, evidence) {
  vapply(seq_len(ncol(evidence)), function(j) {
    rc_gpr_confidence_one(gpr, evidence[, j])
  }, numeric(1))
}

rc_reaction_confidence <- function(gpr_list,
                                   gene_confidence = NULL,
                                   unit_detection = NULL,
                                   method = c("gpr_aware"),
                                   tau_conf = 0.20,
                                   and_method = c("softmin", "min", "mean"),
                                   or_method = c("max", "prob_or", "sum_sqrtK"),
                                   low_confidence_quantile = NULL,
                                   low_confidence_threshold = NULL,
                                   unit_ids = NULL) {
  method <- match.arg(method)
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)

  if (identical(method, "gpr_aware")) {
    evidence <- if (!is.null(gene_confidence)) gene_confidence else unit_detection
    source <- if (!is.null(gene_confidence)) "gpr_aware_gene_confidence" else "gpr_aware_rna_detection"
    detection_available <- !is.null(unit_detection)
    if (is.null(evidence)) {
      if (is.null(unit_ids)) stop("Provide `gene_confidence` or `unit_detection` for gpr_aware confidence.", call. = FALSE)
      evidence <- rc_empty_gpr_evidence_matrix(gpr_list, unit_ids)
      source <- "gpr_aware_no_evidence"
      detection_available <- FALSE
    }
    out <- rc_reaction_confidence_gpr_aware(
      gpr_list = gpr_list,
      gene_confidence = evidence,
      unit_detection = unit_detection,
      tau_conf = tau_conf,
      and_method = and_method,
      or_method = or_method,
      low_confidence_quantile = low_confidence_quantile,
      low_confidence_threshold = low_confidence_threshold
    )
    return(rc_set_confidence_source(out, source, detection_available = detection_available))
  }
}

rc_gpr_gene_ids <- function(gpr_list) {
  genes <- unique(tolower(unlist(gpr_list, use.names = FALSE)))
  genes[!is.na(genes) & nzchar(genes)]
}

rc_metabolic_gpr_genes <- function(gpr_table) {
  if (is.data.frame(gpr_table) && all(c("reaction_id", "and_group_id", "gene") %in% colnames(gpr_table))) {
    genes <- unique(as.character(gpr_table$gene))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    return(genes)
  }
  if (is.data.frame(gpr_table) && all(c("reaction_id", "gpr") %in% colnames(gpr_table))) {
    rules <- as.character(gpr_table$gpr)
    rules <- gsub("[()]", " ", rules)
    rules <- gsub("\\s+(and|or)\\s+", "|", rules, perl = TRUE, ignore.case = TRUE)
    tokens <- unlist(strsplit(rules, "|", fixed = TRUE), use.names = FALSE)
    genes <- unique(trimws(tokens))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    return(genes)
  }
  gpr_list <- if (is.data.frame(gpr_table)) rc_parse_gpr_table(gpr_table) else gpr_table
  genes <- unique(as.character(unlist(gpr_list, use.names = FALSE)))
  genes[!is.na(genes) & nzchar(genes)]
}

rc_safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}

rc_q95_bootstrap_diagnostics <- function(C_raw, Q, unit_meta = NULL, stratum_col = NULL, B = 500, BPPARAM = NULL) {
  C_raw <- as.matrix(C_raw)
  strata <- if (!is.null(stratum_col)) {
    unit_meta <- unit_meta[match(colnames(C_raw), unit_meta$pool_id), , drop = FALSE]
    as.character(unit_meta[[stratum_col]])
  } else rep("global", ncol(C_raw))
  rows <- seq_len(nrow(Q))
  boot <- rc_parallel_lapply(rows, function(i) {
    rid <- Q$reaction_id[[i]]; st <- Q$stratum[[i]]
    cols <- which(strata == st)
    rc_q95_bootstrap(C_raw[rid, cols], B = B)
  }, BPPARAM = BPPARAM)
  do.call(rbind, boot)
}
