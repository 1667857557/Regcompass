# Shared numerical utilities for the canonical multitask GRN estimator.
#
# Authoritative parameter defaults and validation are defined once in
# multitask_grn_parameter_policy.R. The direct condition-theta estimator is
# defined once in multitask_grn_cv_contract.R. This file deliberately contains
# no estimator/default shadow definitions; Collate order therefore cannot change
# the mathematical model by rebinding the same function names.

.rc_residualize_matrix <- function(x, condition) {
  x <- as.matrix(x)
  condition <- as.character(condition)
  if (nrow(x) != length(condition)) {
    stop("Condition labels must match matrix rows.", call. = FALSE)
  }
  out <- x
  for (index in split(seq_len(nrow(x)), condition)) {
    out[index, ] <- sweep(
      x[index, , drop = FALSE],
      2L,
      colMeans(x[index, , drop = FALSE]),
      "-"
    )
  }
  out
}

.rc_residualize_vector <- function(x, condition) {
  x <- as.numeric(x)
  condition <- as.character(condition)
  if (length(x) != length(condition)) {
    stop("Condition labels must match vector length.", call. = FALSE)
  }
  out <- x
  for (index in split(seq_along(x), condition)) {
    out[index] <- x[index] - mean(x[index])
  }
  out
}

.rc_condition_balanced_weights <- function(condition) {
  condition <- as.character(condition)
  if (!length(condition) || anyNA(condition) || any(!nzchar(condition))) {
    stop("Condition labels must be complete and non-empty.", call. = FALSE)
  }
  count <- table(condition)
  weight <- 1 / as.numeric(count[condition])
  weight / mean(weight)
}

.rc_multitask_foldid <- function(condition, nfolds = 5L, seed = 1L) {
  condition <- as.character(condition)
  count <- table(condition)
  k <- min(as.integer(nfolds), min(count))
  if (k < 3L) {
    stop(
      "At least three cells per condition are required for cross-validation.",
      call. = FALSE
    )
  }
  set.seed(seed)
  foldid <- integer(length(condition))
  for (level in names(count)) {
    index <- which(condition == level)
    index <- base::sample(index, length(index), replace = FALSE)
    foldid[index] <- rep(seq_len(k), length.out = length(index))
  }
  foldid
}

.rc_condition_stratified_bootstrap_indices <- function(condition) {
  condition <- as.character(condition)
  unlist(lapply(split(seq_along(condition), condition), function(index) {
    base::sample(index, length(index), replace = TRUE)
  }), use.names = FALSE)
}

.rc_edge_screen_score <- function(x, y, condition) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  condition <- as.character(condition)
  if (nrow(x) != length(y) || length(y) != length(condition)) {
    stop("Predictors, target, and condition labels must align.", call. = FALSE)
  }
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

.rc_extract_glmnet_vector <- function(fit, s, expected) {
  value <- as.matrix(stats::coef(fit, s = s))[, 1L]
  if ("(Intercept)" %in% names(value)) {
    value <- value[names(value) != "(Intercept)"]
  }
  value <- as.numeric(value)
  if (length(value) != expected) {
    stop("glmnet returned an unexpected coefficient vector length.",
         call. = FALSE)
  }
  value
}

.rc_weighted_rsq <- function(observed, predicted, weight) {
  keep <- is.finite(observed) & is.finite(predicted) &
    is.finite(weight) & weight > 0
  if (sum(keep) < 3L) return(NA_real_)
  observed <- observed[keep]
  predicted <- predicted[keep]
  weight <- weight[keep]
  center <- sum(weight * observed) / sum(weight)
  denominator <- sum(weight * (observed - center)^2)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  1 - sum(weight * (observed - predicted)^2) / denominator
}
