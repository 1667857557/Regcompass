.rc_multitask_grn_defaults <- function() {
  list(
    alpha = 0.5,
    global_penalty_factor = 1,
    deviation_penalty_factor = 2,
    lambda_rule = "lambda.1se",
    nfolds = 5L,
    n_stability = 25L,
    stability_fraction = 0.8,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    min_abs_effect = 0,
    min_cv_rsq = 0,
    zero_tolerance = 1e-8,
    candidate_screen_threshold = 0,
    max_edges_per_target = 500L,
    seed = 12345L
  )
}

.rc_validate_multitask_grn_args <- function(args = list()) {
  if (!is.list(args)) {
    stop("`multitask_args` must be a list.", call. = FALSE)
  }
  defaults <- .rc_multitask_grn_defaults()
  unknown <- setdiff(names(args), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown `multitask_args`: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  out <- modifyList(defaults, args)
  fraction <- c(
    "alpha", "stability_fraction", "min_selection_frequency",
    "min_sign_stability"
  )
  for (name in fraction) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 1) {
      stop("`multitask_args$", name, "` must be in [0, 1].", call. = FALSE)
    }
  }
  positive <- c("global_penalty_factor", "deviation_penalty_factor")
  for (name in positive) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop("`multitask_args$", name, "` must be positive.", call. = FALSE)
    }
  }
  nonnegative <- c(
    "min_abs_effect", "min_cv_rsq", "zero_tolerance",
    "candidate_screen_threshold"
  )
  for (name in nonnegative) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0) {
      stop("`multitask_args$", name, "` must be non-negative.", call. = FALSE)
    }
  }
  integer_fields <- c("nfolds", "n_stability", "max_edges_per_target", "seed")
  for (name in integer_fields) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < if (name == "n_stability") 0 else 1 ||
        abs(value - round(value)) > sqrt(.Machine$double.eps)) {
      stop("`multitask_args$", name, "` must be an integer in its valid range.",
           call. = FALSE)
    }
    out[[name]] <- as.integer(value)
  }
  out$lambda_rule <- match.arg(
    as.character(out$lambda_rule), c("lambda.1se", "lambda.min")
  )
  if (out$nfolds < 3L) {
    stop("`multitask_args$nfolds` must be at least 3.", call. = FALSE)
  }
  out
}

.rc_multitask_contrast <- function(condition) {
  levels <- sort(unique(as.character(condition)))
  if (!length(levels)) stop("No conditions were available.", call. = FALSE)
  if (length(levels) == 1L) {
    contrast <- matrix(numeric(), nrow = 1L, ncol = 0L)
    rownames(contrast) <- levels
    return(contrast)
  }
  contrast <- stats::contr.sum(length(levels))
  rownames(contrast) <- levels
  colnames(contrast) <- paste0("condition_deviation_", seq_len(ncol(contrast)))
  contrast
}

.rc_residualize_matrix <- function(x, block) {
  x <- as.matrix(x)
  block <- as.character(block)
  out <- x
  rows <- split(seq_len(nrow(x)), block)
  for (index in rows) {
    out[index, ] <- sweep(
      x[index, , drop = FALSE], 2L,
      colMeans(x[index, , drop = FALSE]), "-"
    )
  }
  out
}

.rc_residualize_vector <- function(x, block) {
  x <- as.numeric(x)
  block <- as.character(block)
  out <- x
  rows <- split(seq_along(x), block)
  for (index in rows) out[index] <- x[index] - mean(x[index])
  out
}

.rc_condition_balanced_weights <- function(condition) {
  condition <- as.character(condition)
  n <- table(condition)
  weight <- 1 / as.numeric(n[condition])
  weight / mean(weight)
}

.rc_multitask_foldid <- function(condition, sample = NULL, nfolds = 5L, seed = 1L) {
  condition <- as.character(condition)
  set.seed(seed)
  if (!is.null(sample)) {
    sample <- paste(condition, as.character(sample), sep = "\001")
    samples_by_condition <- split(unique(sample), condition[match(unique(sample), sample)])
    min_samples <- min(vapply(samples_by_condition, length, integer(1)))
    if (is.finite(min_samples) && min_samples >= 3L) {
      k <- min(nfolds, min_samples)
      sample_fold <- integer()
      for (values in samples_by_condition) {
        values <- sample(values, length(values))
        fold <- rep(seq_len(k), length.out = length(values))
        sample_fold[values] <- fold
      }
      return(as.integer(sample_fold[sample]))
    }
  }
  count <- table(condition)
  k <- min(nfolds, min(count))
  if (k < 3L) {
    stop("At least three cells per condition are required for cross-validation.",
         call. = FALSE)
  }
  foldid <- integer(length(condition))
  for (level in names(count)) {
    index <- which(condition == level)
    index <- sample(index, length(index))
    foldid[index] <- rep(seq_len(k), length.out = length(index))
  }
  foldid
}

.rc_edge_screen_score <- function(x, y, condition) {
  condition <- as.character(condition)
  scores <- rep(0, ncol(x))
  for (level in unique(condition)) {
    index <- which(condition == level)
    if (length(index) < 3L) next
    y_one <- y[index]
    if (!is.finite(stats::sd(y_one)) || stats::sd(y_one) <= 0) next
    value <- suppressWarnings(stats::cor(x[index, , drop = FALSE], y_one))
    value[!is.finite(value)] <- 0
    scores <- pmax(scores, abs(as.numeric(value)))
  }
  scores
}

.rc_decode_multitask_coefficients <- function(coef, n_edges, contrast) {
  beta <- coef[seq_len(n_edges)]
  if (!ncol(contrast)) {
    delta <- matrix(0, nrow = nrow(contrast), ncol = n_edges)
  } else {
    gamma <- coef[n_edges + seq_len(n_edges * ncol(contrast))]
    gamma <- matrix(gamma, nrow = n_edges, ncol = ncol(contrast))
    delta <- contrast %*% t(gamma)
  }
  theta <- sweep(delta, 2L, beta, "+")
  rownames(delta) <- rownames(contrast)
  rownames(theta) <- rownames(contrast)
  list(beta = beta, delta = delta, theta = theta)
}

.rc_extract_glmnet_vector <- function(fit, s, expected) {
  value <- as.matrix(stats::coef(fit, s = s))[, 1L]
  if ("(Intercept)" %in% names(value)) value <- value[names(value) != "(Intercept)"]
  value <- as.numeric(value)
  if (length(value) != expected) {
    stop("glmnet returned an unexpected coefficient vector length.", call. = FALSE)
  }
  value
}

.rc_weighted_rsq <- function(observed, predicted, weight) {
  keep <- is.finite(observed) & is.finite(predicted) & is.finite(weight) & weight > 0
  if (sum(keep) < 3L) return(NA_real_)
  observed <- observed[keep]
  predicted <- predicted[keep]
  weight <- weight[keep]
  center <- sum(weight * observed) / sum(weight)
  denominator <- sum(weight * (observed - center)^2)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  1 - sum(weight * (observed - predicted)^2) / denominator
}

.rc_fit_multitask_target <- function(
    edges, target, rna, atac, meta, condition_col, sample_col, args) {
  empty <- list(
    global = data.frame(), condition = data.frame(), diagnostics = data.frame(
      target = target, status = "no_candidates", stringsAsFactors = FALSE
    )
  )
  if (!nrow(edges)) return(empty)
  cells <- colnames(rna)
  condition <- as.character(meta[[condition_col]])
  sample <- if (!is.null(sample_col)) as.character(meta[[sample_col]]) else NULL
  y <- as.numeric(rna[target, cells, drop = TRUE])
  tf <- as.matrix(Matrix::t(rna[edges$tf_feature_id, cells, drop = FALSE]))
  peak <- as.matrix(Matrix::t(atac[edges$atac_feature_id, cells, drop = FALSE]))
  x_raw <- tf * peak
  x_raw[!is.finite(x_raw)] <- 0
  y[!is.finite(y)] <- 0

  screen_score <- .rc_edge_screen_score(x_raw, y, condition)
  keep <- is.finite(screen_score) &
    screen_score >= args$candidate_screen_threshold
  if (!any(keep)) {
    empty$diagnostics$status <- "no_edges_after_shared_screen"
    return(empty)
  }
  edges <- edges[keep, , drop = FALSE]
  x_raw <- x_raw[, keep, drop = FALSE]
  tf <- tf[, keep, drop = FALSE]
  screen_score <- screen_score[keep]
  order_index <- order(-screen_score, edges$edge_id)
  if (length(order_index) > args$max_edges_per_target) {
    order_index <- order_index[seq_len(args$max_edges_per_target)]
  }
  edges <- edges[order_index, , drop = FALSE]
  x_raw <- x_raw[, order_index, drop = FALSE]
  tf <- tf[, order_index, drop = FALSE]
  screen_score <- screen_score[order_index]

  residual_block <- if (!is.null(sample)) sample else condition
  x_residual <- .rc_residualize_matrix(x_raw, residual_block)
  y_residual <- .rc_residualize_vector(y, residual_block)
  condition_rows <- split(seq_along(condition), condition)
  edge_scale <- sqrt(Reduce(`+`, lapply(condition_rows, function(index) {
    colMeans(x_residual[index, , drop = FALSE]^2)
  })) / length(condition_rows))
  variable <- is.finite(edge_scale) & edge_scale > args$zero_tolerance
  if (!any(variable)) {
    empty$diagnostics$status <- "no_within_group_edge_variation"
    return(empty)
  }
  edges <- edges[variable, , drop = FALSE]
  x_residual <- x_residual[, variable, drop = FALSE]
  tf <- tf[, variable, drop = FALSE]
  screen_score <- screen_score[variable]
  edge_scale <- edge_scale[variable]
  x_scaled <- sweep(x_residual, 2L, edge_scale, "/")

  tf_reference <- Reduce(`+`, lapply(condition_rows, function(index) {
    colMeans(tf[index, , drop = FALSE])
  })) / length(condition_rows)
  contrast <- .rc_multitask_contrast(condition)
  contrast_rows <- contrast[condition, , drop = FALSE]
  design_blocks <- list(x_scaled)
  if (ncol(contrast_rows)) {
    for (j in seq_len(ncol(contrast_rows))) {
      design_blocks[[length(design_blocks) + 1L]] <-
        x_scaled * as.numeric(contrast_rows[, j])
    }
  }
  design <- do.call(cbind, design_blocks)
  n_edges <- ncol(x_scaled)
  penalty_factor <- c(
    rep(args$global_penalty_factor, n_edges),
    rep(args$deviation_penalty_factor, n_edges * ncol(contrast))
  )
  weight <- .rc_condition_balanced_weights(condition)
  foldid <- .rc_multitask_foldid(
    condition = condition,
    sample = sample,
    nfolds = args$nfolds,
    seed = args$seed + sum(utf8ToInt(target))
  )
  cvfit <- tryCatch(
    glmnet::cv.glmnet(
      x = design,
      y = y_residual,
      weights = weight,
      foldid = foldid,
      family = "gaussian",
      alpha = args$alpha,
      intercept = FALSE,
      standardize = FALSE,
      penalty.factor = penalty_factor,
      keep = TRUE,
      parallel = FALSE
    ),
    error = function(error) error
  )
  if (inherits(cvfit, "error")) {
    empty$diagnostics$status <- "glmnet_failed"
    empty$diagnostics$error_message <- conditionMessage(cvfit)
    return(empty)
  }
  lambda <- cvfit[[args$lambda_rule]]
  expected <- ncol(design)
  coef <- .rc_extract_glmnet_vector(cvfit, args$lambda_rule, expected)
  decoded <- .rc_decode_multitask_coefficients(coef, n_edges, contrast)
  lambda_index <- which.min(abs(cvfit$lambda - lambda))
  prediction <- cvfit$fit.preval[, lambda_index]
  cv_rsq <- .rc_weighted_rsq(y_residual, prediction, weight)

  condition_levels <- rownames(contrast)
  selection_count <- matrix(0, nrow = length(condition_levels), ncol = n_edges)
  sign_sum <- matrix(0, nrow = length(condition_levels), ncol = n_edges)
  n_stability_success <- 0L
  if (args$n_stability > 0L) {
    set.seed(args$seed + 31L + sum(utf8ToInt(target)))
    for (b in seq_len(args$n_stability)) {
      selected_cells <- unlist(lapply(condition_rows, function(index) {
        size <- max(3L, floor(length(index) * args$stability_fraction))
        sample(index, min(length(index), size), replace = FALSE)
      }), use.names = FALSE)
      fit <- tryCatch(
        glmnet::glmnet(
          x = design[selected_cells, , drop = FALSE],
          y = y_residual[selected_cells],
          weights = weight[selected_cells],
          family = "gaussian",
          alpha = args$alpha,
          lambda = lambda,
          intercept = FALSE,
          standardize = FALSE,
          penalty.factor = penalty_factor
        ),
        error = function(error) NULL
      )
      if (is.null(fit)) next
      value <- tryCatch(
        .rc_extract_glmnet_vector(fit, lambda, expected),
        error = function(error) NULL
      )
      if (is.null(value)) next
      theta <- .rc_decode_multitask_coefficients(
        value, n_edges, contrast
      )$theta
      selected <- abs(theta) > args$zero_tolerance
      selection_count <- selection_count + selected
      sign_sum <- sign_sum + sign(theta) * selected
      n_stability_success <- n_stability_success + 1L
    }
  }
  if (n_stability_success > 0L) {
    selection_frequency <- selection_count / n_stability_success
    sign_stability <- matrix(0, nrow = nrow(selection_count), ncol = ncol(selection_count))
    nonzero <- selection_count > 0
    sign_stability[nonzero] <- abs(sign_sum[nonzero]) / selection_count[nonzero]
  } else {
    selected <- abs(decoded$theta) > args$zero_tolerance
    selection_frequency <- selected * 1
    sign_stability <- selected * 1
  }
  stable_estimate <- decoded$theta * selection_frequency * sign_stability
  active <- selection_frequency >= args$min_selection_frequency &
    sign_stability >= args$min_sign_stability &
    abs(decoded$theta) >= args$min_abs_effect &
    is.finite(cv_rsq) & cv_rsq >= args$min_cv_rsq

  global <- cbind(
    edges,
    data.frame(
      global_estimate = as.numeric(decoded$beta),
      candidate_screen_score = screen_score,
      edge_scale = edge_scale,
      tf_reference = tf_reference,
      cv_rsq = cv_rsq,
      lambda = lambda,
      stringsAsFactors = FALSE
    )
  )
  condition_output <- do.call(rbind, lapply(seq_along(condition_levels), function(i) {
    data.frame(
      edges,
      condition = condition_levels[[i]],
      global_estimate = as.numeric(decoded$beta),
      condition_deviation = as.numeric(decoded$delta[i, ]),
      effective_estimate = as.numeric(decoded$theta[i, ]),
      selection_frequency = as.numeric(selection_frequency[i, ]),
      sign_stability = as.numeric(sign_stability[i, ]),
      stability_weight = as.numeric(selection_frequency[i, ] * sign_stability[i, ]),
      stable_estimate = as.numeric(stable_estimate[i, ]),
      estimate = as.numeric(stable_estimate[i, ]),
      active_edge = as.logical(active[i, ]),
      candidate_screen_score = screen_score,
      edge_scale = edge_scale,
      tf_reference = tf_reference,
      atac_projection_weight = as.numeric(
        stable_estimate[i, ] * tf_reference / pmax(edge_scale, args$zero_tolerance)
      ),
      cv_rsq = cv_rsq,
      rsq = cv_rsq,
      lambda = lambda,
      padj = NA_real_,
      evidence_type = "multitask_stability_selected",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  sign_by_edge <- split(
    condition_output,
    condition_output$edge_id
  )
  flip <- vapply(sign_by_edge, function(one) {
    one <- one[one$active_edge %in% TRUE, , drop = FALSE]
    length(unique(sign(one$effective_estimate))) > 1L
  }, logical(1))
  condition_output$sign_flip_flag <- unname(flip[condition_output$edge_id])

  diagnostics <- data.frame(
    target = target,
    status = "ok",
    n_structural_candidates = nrow(edges),
    n_model_edges = n_edges,
    n_active_condition_edges = sum(active),
    n_active_conditions = sum(rowSums(active) > 0),
    cv_rsq = cv_rsq,
    lambda = lambda,
    lambda_rule = args$lambda_rule,
    n_stability_requested = args$n_stability,
    n_stability_success = n_stability_success,
    residualization_block = if (is.null(sample)) "condition" else "sample",
    stringsAsFactors = FALSE
  )
  list(global = global, condition = condition_output, diagnostics = diagnostics)
}

.rc_fit_multitask_celltype_grn <- function(
    design, rna, atac, meta, condition_col, sample_col = NULL,
    multitask_args = list()) {
  args <- .rc_validate_multitask_grn_args(multitask_args)
  Pando::validate_grn_design(design)
  cells <- as.character(design$feature_contract$cell_ids)
  if (!setequal(cells, colnames(rna)) || !setequal(cells, colnames(atac))) {
    stop("Pando design cells and normalized assays do not match.", call. = FALSE)
  }
  rna <- rna[, cells, drop = FALSE]
  atac <- atac[, cells, drop = FALSE]
  meta <- meta[match(cells, rownames(meta)), , drop = FALSE]
  if (anyNA(rownames(meta))) stop("Cell metadata do not align to the Pando design.", call. = FALSE)
  candidates <- design$candidate_edges
  if (!nrow(candidates)) {
    stop("Pando produced no structural TF-peak-target candidates.", call. = FALSE)
  }
  missing_rna <- setdiff(
    unique(c(candidates$tf_feature_id, candidates$target_feature_id)),
    rownames(rna)
  )
  missing_atac <- setdiff(unique(candidates$atac_feature_id), rownames(atac))
  if (length(missing_rna) || length(missing_atac)) {
    stop("Pando design feature IDs are absent from normalized assay matrices.",
         call. = FALSE)
  }
  targets <- unique(as.character(candidates$target))
  fits <- lapply(targets, function(target) {
    .rc_fit_multitask_target(
      edges = candidates[candidates$target == target, , drop = FALSE],
      target = target,
      rna = rna,
      atac = atac,
      meta = meta,
      condition_col = condition_col,
      sample_col = sample_col,
      args = args
    )
  })
  bind <- function(name) {
    value <- do.call(rbind, lapply(fits, `[[`, name))
    if (is.null(value)) data.frame() else value
  }
  list(
    global = bind("global"),
    condition = bind("condition"),
    diagnostics = bind("diagnostics"),
    params = args
  )
}
