# Leakage-resistant direct condition-theta multitask regression.
#
# The canonical fitter applies the elastic-net penalty directly to the
# condition-specific coefficients theta[e, c]. The shared backbone and zero-sum
# deviations are derived after fitting:
#   beta[e] = mean_c theta[e, c]
#   delta[e, c] = theta[e, c] - beta[e]
# This permits exact condition-specific zeros while preserving one candidate
# universe, one edge scale and one lambda-selection rule across conditions.

.rc_condition_matrix_centers <- function(x, condition) {
  x <- as.matrix(x)
  condition <- as.character(condition)
  levels <- sort(unique(condition))
  centers <- do.call(rbind, lapply(levels, function(level) {
    colMeans(x[condition == level, , drop = FALSE])
  }))
  rownames(centers) <- levels
  colnames(centers) <- colnames(x)
  centers
}

.rc_condition_vector_centers <- function(x, condition) {
  x <- as.numeric(x)
  condition <- as.character(condition)
  vapply(sort(unique(condition)), function(level) {
    mean(x[condition == level])
  }, numeric(1))
}

.rc_apply_matrix_centers <- function(x, condition, centers) {
  x <- as.matrix(x)
  condition <- as.character(condition)
  missing <- setdiff(unique(condition), rownames(centers))
  if (length(missing)) {
    stop(
      "Validation conditions are absent from the training fold: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  x - centers[condition, , drop = FALSE]
}

.rc_apply_vector_centers <- function(x, condition, centers) {
  x <- as.numeric(x)
  condition <- as.character(condition)
  missing <- setdiff(unique(condition), names(centers))
  if (length(missing)) {
    stop(
      "Validation conditions are absent from the training fold: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  x - unname(centers[condition])
}

.rc_equal_condition_edge_scale <- function(x_centered, condition) {
  condition <- as.character(condition)
  rows <- split(seq_along(condition), condition)
  sqrt(Reduce(`+`, lapply(rows, function(index) {
    colMeans(x_centered[index, , drop = FALSE]^2)
  })) / length(rows))
}

.rc_order_target_edges <- function(edges) {
  if (!is.data.frame(edges) || !nrow(edges)) return(edges)
  if (!"edge_id" %in% colnames(edges) || anyNA(edges$edge_id) ||
      any(!nzchar(trimws(as.character(edges$edge_id)))) ||
      anyDuplicated(as.character(edges$edge_id))) {
    stop("Target candidate edges require unique non-empty `edge_id` values.",
         call. = FALSE)
  }
  order_index <- order(as.character(edges$edge_id))
  out <- edges[order_index, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.rc_condition_theta_design_matrix <- function(
    x_scaled, condition, condition_levels = NULL) {
  x_scaled <- as.matrix(x_scaled)
  condition <- as.character(condition)
  if (is.null(condition_levels)) {
    condition_levels <- sort(unique(condition))
  } else {
    condition_levels <- as.character(condition_levels)
  }
  if (length(setdiff(unique(condition), condition_levels))) {
    stop("Condition-specific design levels do not cover all observations.",
         call. = FALSE)
  }
  edge_names <- colnames(x_scaled)
  if (is.null(edge_names)) edge_names <- paste0("edge_", seq_len(ncol(x_scaled)))
  blocks <- lapply(condition_levels, function(level) {
    block <- x_scaled * as.numeric(condition == level)
    colnames(block) <- paste(level, edge_names, sep = "::")
    block
  })
  do.call(cbind, blocks)
}

.rc_decode_condition_theta <- function(
    coefficient, n_edges, condition_levels, edge_ids = NULL) {
  condition_levels <- as.character(condition_levels)
  expected <- n_edges * length(condition_levels)
  coefficient <- as.numeric(coefficient)
  if (length(coefficient) != expected) {
    stop("Condition-theta coefficient vector has an unexpected length.",
         call. = FALSE)
  }
  theta <- t(matrix(
    coefficient,
    nrow = n_edges,
    ncol = length(condition_levels)
  ))
  rownames(theta) <- condition_levels
  if (!is.null(edge_ids)) colnames(theta) <- as.character(edge_ids)
  beta <- colMeans(theta)
  delta <- sweep(theta, 2L, beta, "-")
  list(beta = beta, delta = delta, theta = theta)
}

.rc_multitask_manual_cv <- function(
    x_raw, y, condition, penalty_factor, args, full_edge_scale) {
  condition <- as.character(condition)
  condition_levels <- sort(unique(condition))
  full_x_centered <- .rc_residualize_matrix(x_raw, condition)
  full_y_centered <- .rc_residualize_vector(y, condition)
  full_x_scaled <- sweep(full_x_centered, 2L, full_edge_scale, "/")
  full_design <- .rc_condition_theta_design_matrix(
    full_x_scaled, condition, condition_levels
  )
  full_weight <- .rc_condition_balanced_weights(condition)
  full_path <- glmnet::glmnet(
    x = full_design,
    y = full_y_centered,
    weights = full_weight,
    family = "gaussian",
    alpha = args$alpha,
    intercept = FALSE,
    standardize = FALSE,
    penalty.factor = penalty_factor
  )
  lambda <- as.numeric(full_path$lambda)
  if (!length(lambda) || any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("glmnet did not return a valid common lambda path.", call. = FALSE)
  }

  foldid <- .rc_multitask_foldid(
    condition = condition,
    nfolds = args$nfolds,
    seed = args$seed
  )
  folds <- sort(unique(foldid))
  n_lambda <- length(lambda)
  oof_prediction <- matrix(NA_real_, nrow = length(y), ncol = n_lambda)
  oof_observed <- rep(NA_real_, length(y))
  fold_mse <- matrix(NA_real_, nrow = length(folds), ncol = n_lambda)

  for (fold_index in seq_along(folds)) {
    fold <- folds[[fold_index]]
    validation <- which(foldid == fold)
    training <- which(foldid != fold)
    train_condition <- condition[training]
    validation_condition <- condition[validation]

    x_centers <- .rc_condition_matrix_centers(
      x_raw[training, , drop = FALSE], train_condition
    )
    y_centers <- .rc_condition_vector_centers(y[training], train_condition)
    x_train <- .rc_apply_matrix_centers(
      x_raw[training, , drop = FALSE], train_condition, x_centers
    )
    y_train <- .rc_apply_vector_centers(
      y[training], train_condition, y_centers
    )
    x_validation <- .rc_apply_matrix_centers(
      x_raw[validation, , drop = FALSE], validation_condition, x_centers
    )
    y_validation <- .rc_apply_vector_centers(
      y[validation], validation_condition, y_centers
    )

    train_scale <- .rc_equal_condition_edge_scale(x_train, train_condition)
    invalid_scale <- !is.finite(train_scale) |
      train_scale <= args$zero_tolerance
    train_scale[invalid_scale] <- full_edge_scale[invalid_scale]
    train_scale[!is.finite(train_scale) |
                  train_scale <= args$zero_tolerance] <- 1
    x_train <- sweep(x_train, 2L, train_scale, "/")
    x_validation <- sweep(x_validation, 2L, train_scale, "/")

    train_design <- .rc_condition_theta_design_matrix(
      x_train, train_condition, condition_levels
    )
    validation_design <- .rc_condition_theta_design_matrix(
      x_validation, validation_condition, condition_levels
    )
    train_weight <- .rc_condition_balanced_weights(train_condition)
    fit <- glmnet::glmnet(
      x = train_design,
      y = y_train,
      weights = train_weight,
      family = "gaussian",
      alpha = args$alpha,
      lambda = lambda,
      intercept = FALSE,
      standardize = FALSE,
      penalty.factor = penalty_factor
    )
    prediction <- as.matrix(stats::predict(
      fit,
      newx = validation_design,
      s = lambda,
      type = "response"
    ))
    if (!identical(dim(prediction), c(length(validation), n_lambda))) {
      stop("glmnet returned an unexpected validation prediction matrix.",
           call. = FALSE)
    }
    oof_prediction[validation, ] <- prediction
    oof_observed[validation] <- y_validation
    validation_weight <- full_weight[validation]
    squared <- sweep(prediction, 1L, y_validation, "-")^2
    fold_mse[fold_index, ] <- colSums(
      squared * validation_weight
    ) / sum(validation_weight)
  }

  if (anyNA(oof_prediction) || anyNA(oof_observed)) {
    stop("Manual cross-validation did not produce complete out-of-fold values.",
         call. = FALSE)
  }
  squared <- sweep(oof_prediction, 1L, oof_observed, "-")^2
  cvm <- colSums(squared * full_weight) / sum(full_weight)
  cvsd <- apply(fold_mse, 2L, stats::sd) / sqrt(nrow(fold_mse))
  cvsd[!is.finite(cvsd)] <- 0
  lambda_min_index <- which.min(cvm)
  eligible <- which(cvm <= cvm[[lambda_min_index]] + cvsd[[lambda_min_index]])
  lambda_1se_index <- eligible[[which.max(lambda[eligible])]]
  selected_index <- if (identical(args$lambda_rule, "lambda.min")) {
    lambda_min_index
  } else {
    lambda_1se_index
  }
  selected_lambda <- lambda[[selected_index]]
  coefficient <- .rc_extract_glmnet_vector(
    full_path, selected_lambda, ncol(full_design)
  )

  list(
    fit = full_path,
    coefficient = coefficient,
    lambda = selected_lambda,
    lambda_index = selected_index,
    lambda_grid = lambda,
    cvm = cvm,
    cvsd = cvsd,
    foldid = foldid,
    oof_observed = oof_observed,
    oof_prediction = oof_prediction[, selected_index],
    full_x_centered = full_x_centered,
    full_y_centered = full_y_centered,
    full_x_scaled = full_x_scaled,
    full_design = full_design,
    full_weight = full_weight,
    condition_levels = condition_levels,
    preprocessing =
      "fold_specific_training_condition_centering_and_edge_scaling",
    coefficient_parameterization = "direct_condition_theta"
  )
}

.rc_fit_multitask_target_direct <- function(
    edges, target, rna, atac, meta, condition_col, args) {
  empty <- list(
    global = data.frame(),
    condition = data.frame(),
    diagnostics = data.frame(
      target = target,
      status = "no_candidates",
      stringsAsFactors = FALSE
    )
  )
  if (!nrow(edges)) return(empty)
  edges <- .rc_order_target_edges(edges)
  cells <- colnames(rna)
  condition <- trimws(as.character(meta[[condition_col]]))
  y <- as.numeric(rna[target, cells, drop = TRUE])
  tf <- as.matrix(Matrix::t(
    rna[edges$tf_feature_id, cells, drop = FALSE]
  ))
  peak <- as.matrix(Matrix::t(
    atac[edges$atac_feature_id, cells, drop = FALSE]
  ))
  colnames(tf) <- edges$edge_id
  colnames(peak) <- edges$edge_id
  x_raw <- tf * peak
  colnames(x_raw) <- edges$edge_id
  x_raw[!is.finite(x_raw)] <- 0
  y[!is.finite(y)] <- 0

  screen_score <- .rc_edge_screen_score(x_raw, y, condition)
  if (nrow(edges) > args$max_edges_per_target) {
    keep <- seq_len(args$max_edges_per_target)
    edges <- edges[keep, , drop = FALSE]
    x_raw <- x_raw[, keep, drop = FALSE]
    tf <- tf[, keep, drop = FALSE]
    screen_score <- screen_score[keep]
  }

  x_residual <- .rc_residualize_matrix(x_raw, condition)
  y_residual <- .rc_residualize_vector(y, condition)
  edge_scale <- .rc_equal_condition_edge_scale(x_residual, condition)
  variable <- is.finite(edge_scale) & edge_scale > args$zero_tolerance
  if (!any(variable)) {
    empty$diagnostics$status <- "no_within_condition_edge_variation"
    return(empty)
  }
  edges <- edges[variable, , drop = FALSE]
  x_raw <- x_raw[, variable, drop = FALSE]
  x_residual <- x_residual[, variable, drop = FALSE]
  tf <- tf[, variable, drop = FALSE]
  screen_score <- screen_score[variable]
  edge_scale <- edge_scale[variable]

  condition_rows <- split(seq_along(condition), condition)
  tf_reference <- Reduce(`+`, lapply(condition_rows, function(index) {
    colMeans(tf[index, , drop = FALSE])
  })) / length(condition_rows)
  condition_levels <- sort(unique(condition))
  n_edges <- ncol(x_raw)
  if (!isTRUE(all.equal(
    as.numeric(args$global_penalty_factor),
    as.numeric(args$deviation_penalty_factor),
    tolerance = sqrt(.Machine$double.eps)
  ))) {
    stop(
      paste(
        "Direct condition-theta sparsity requires equal global and deviation",
        "penalty factors. These legacy fields now act as one common theta",
        "penalty; set both to the same value."
      ),
      call. = FALSE
    )
  }
  theta_penalty_factor <- as.numeric(args$global_penalty_factor)
  penalty_factor <- rep(
    theta_penalty_factor,
    n_edges * length(condition_levels)
  )
  cvfit <- tryCatch(
    .rc_multitask_manual_cv(
      x_raw = x_raw,
      y = y,
      condition = condition,
      penalty_factor = penalty_factor,
      args = args,
      full_edge_scale = edge_scale
    ),
    error = function(error) error
  )
  if (inherits(cvfit, "error")) {
    empty$diagnostics$status <- "glmnet_cv_failed"
    empty$diagnostics$error_message <- conditionMessage(cvfit)
    return(empty)
  }
  decoded <- .rc_decode_condition_theta(
    cvfit$coefficient,
    n_edges = n_edges,
    condition_levels = condition_levels,
    edge_ids = edges$edge_id
  )
  cv_rsq <- .rc_weighted_rsq(
    cvfit$oof_observed,
    cvfit$oof_prediction,
    cvfit$full_weight
  )

  selection_count <- matrix(
    0,
    nrow = length(condition_levels),
    ncol = n_edges,
    dimnames = list(condition_levels, edges$edge_id)
  )
  sign_sum <- selection_count
  n_bootstrap_success <- 0L
  set.seed(args$seed + 31L + sum(utf8ToInt(target)))
  for (b in seq_len(args$n_bootstrap)) {
    bootstrap_index <- .rc_condition_stratified_bootstrap_indices(condition)
    bootstrap_condition <- condition[bootstrap_index]
    bootstrap_x <- x_raw[bootstrap_index, , drop = FALSE]
    bootstrap_y <- y[bootstrap_index]
    bootstrap_x <- .rc_residualize_matrix(
      bootstrap_x, bootstrap_condition
    )
    bootstrap_y <- .rc_residualize_vector(
      bootstrap_y, bootstrap_condition
    )
    bootstrap_x <- sweep(bootstrap_x, 2L, edge_scale, "/")
    bootstrap_design <- .rc_condition_theta_design_matrix(
      bootstrap_x, bootstrap_condition, condition_levels
    )
    bootstrap_weight <- .rc_condition_balanced_weights(bootstrap_condition)
    fit <- tryCatch(
      glmnet::glmnet(
        x = bootstrap_design,
        y = bootstrap_y,
        weights = bootstrap_weight,
        family = "gaussian",
        alpha = args$alpha,
        lambda = cvfit$lambda,
        intercept = FALSE,
        standardize = FALSE,
        penalty.factor = penalty_factor
      ),
      error = function(error) NULL
    )
    if (is.null(fit)) next
    value <- tryCatch(
      .rc_extract_glmnet_vector(
        fit, cvfit$lambda, ncol(bootstrap_design)
      ),
      error = function(error) NULL
    )
    if (is.null(value)) next
    theta <- .rc_decode_condition_theta(
      value,
      n_edges = n_edges,
      condition_levels = condition_levels,
      edge_ids = edges$edge_id
    )$theta
    selected <- abs(theta) > args$zero_tolerance
    selection_count <- selection_count + selected
    sign_sum <- sign_sum + sign(theta) * selected
    n_bootstrap_success <- n_bootstrap_success + 1L
  }
  if (n_bootstrap_success > 0L) {
    selection_frequency <- selection_count / n_bootstrap_success
    sign_stability <- matrix(
      0,
      nrow = nrow(selection_count),
      ncol = ncol(selection_count),
      dimnames = dimnames(selection_count)
    )
    nonzero <- selection_count > 0
    sign_stability[nonzero] <-
      abs(sign_sum[nonzero]) / selection_count[nonzero]
  } else {
    selection_frequency <- selection_count
    sign_stability <- selection_count
  }
  stable_estimate <- decoded$theta * selection_frequency * sign_stability
  bootstrap_fraction <- n_bootstrap_success / args$n_bootstrap
  bootstrap_adequate <- is.finite(bootstrap_fraction) &&
    bootstrap_fraction >= args$min_bootstrap_success_fraction
  effect_threshold <- max(args$min_abs_effect, args$zero_tolerance)

  global <- cbind(
    edges,
    data.frame(
      global_estimate = as.numeric(decoded$beta),
      candidate_screen_score = screen_score,
      edge_scale = edge_scale,
      tf_reference = tf_reference,
      cv_rsq = cv_rsq,
      lambda = cvfit$lambda,
      coefficient_parameterization = "direct_condition_theta",
      theta_penalty_factor = theta_penalty_factor,
      cv_preprocessing = cvfit$preprocessing,
      bootstrap_method = "condition_stratified_full_size_nonparametric",
      n_bootstrap_requested = args$n_bootstrap,
      n_bootstrap_success = n_bootstrap_success,
      bootstrap_success_fraction = bootstrap_fraction,
      bootstrap_completion_adequate = bootstrap_adequate,
      stringsAsFactors = FALSE
    )
  )
  condition_output <- do.call(rbind, lapply(
    seq_along(condition_levels), function(i) {
      data.frame(
        edges,
        condition = condition_levels[[i]],
        global_estimate = as.numeric(decoded$beta),
        condition_deviation = as.numeric(decoded$delta[i, ]),
        effective_estimate = as.numeric(decoded$theta[i, ]),
        selection_frequency = as.numeric(selection_frequency[i, ]),
        sign_stability = as.numeric(sign_stability[i, ]),
        stability_weight = as.numeric(
          selection_frequency[i, ] * sign_stability[i, ]
        ),
        stable_estimate = as.numeric(stable_estimate[i, ]),
        estimate = as.numeric(stable_estimate[i, ]),
        active_edge = bootstrap_adequate &
          as.numeric(selection_frequency[i, ]) >= args$min_selection_frequency &
          as.numeric(sign_stability[i, ]) >= args$min_sign_stability &
          abs(as.numeric(decoded$theta[i, ])) > effect_threshold &
          is.finite(cv_rsq) & cv_rsq >= args$min_cv_rsq,
        candidate_screen_score = screen_score,
        edge_scale = edge_scale,
        tf_reference = tf_reference,
        atac_projection_weight = as.numeric(
          stable_estimate[i, ] * tf_reference /
            pmax(edge_scale, args$zero_tolerance)
        ),
        cv_rsq = cv_rsq,
        rsq = cv_rsq,
        lambda = cvfit$lambda,
        coefficient_parameterization = "direct_condition_theta",
        theta_penalty_factor = theta_penalty_factor,
        cv_preprocessing = cvfit$preprocessing,
        bootstrap_method = "condition_stratified_full_size_nonparametric",
        n_bootstrap_requested = args$n_bootstrap,
        n_bootstrap_success = n_bootstrap_success,
        bootstrap_success_fraction = bootstrap_fraction,
        bootstrap_completion_adequate = bootstrap_adequate,
        padj = NA_real_,
        evidence_type = "direct_theta_bootstrap_stability_selected",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  ))
  sign_by_edge <- split(condition_output, condition_output$edge_id)
  flip <- vapply(sign_by_edge, function(one) {
    one <- one[one$active_edge %in% TRUE, , drop = FALSE]
    length(unique(sign(one$effective_estimate))) > 1L
  }, logical(1))
  condition_output$sign_flip_flag <- unname(
    flip[condition_output$edge_id]
  )

  diagnostics <- data.frame(
    target = target,
    status = "ok",
    n_structural_candidates = nrow(edges),
    n_model_edges = n_edges,
    n_active_condition_edges = sum(condition_output$active_edge %in% TRUE),
    n_active_conditions = length(unique(
      condition_output$condition[condition_output$active_edge %in% TRUE]
    )),
    cv_rsq = cv_rsq,
    lambda = cvfit$lambda,
    lambda_rule = args$lambda_rule,
    coefficient_parameterization = "direct_condition_theta",
    theta_penalty_factor = theta_penalty_factor,
    cv_preprocessing = cvfit$preprocessing,
    residualization_block = "condition",
    bootstrap_method = "condition_stratified_full_size_nonparametric",
    n_bootstrap_requested = args$n_bootstrap,
    n_bootstrap_success = n_bootstrap_success,
    bootstrap_success_fraction = bootstrap_fraction,
    min_bootstrap_success_fraction = args$min_bootstrap_success_fraction,
    bootstrap_completion_adequate = bootstrap_adequate,
    active_effect_threshold = effect_threshold,
    stringsAsFactors = FALSE
  )
  list(global = global, condition = condition_output, diagnostics = diagnostics)
}

# One explicit estimator binding replaces the former order/active wrapper chain.
.rc_fit_multitask_target <- .rc_fit_multitask_target_direct
