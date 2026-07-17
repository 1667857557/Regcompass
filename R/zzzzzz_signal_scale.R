# Final signal-scale correction.
# RNA input is log1p(CPM); ATAC input is a non-negative normalized accessibility
# signal. Applying expm1 to both assays over-saturates RNA and is not meaningful
# for TF-IDF-like ATAC values. The canonical transform therefore operates on the
# supplied normalized signal directly.

.rc_absolute_activity_score <- function(X, half_saturation = 1) {
  X <- as.matrix(X)
  if (!is.numeric(half_saturation) || length(half_saturation) != 1L ||
      !is.finite(half_saturation) || half_saturation <= 0) {
    stop("`half_saturation` must be one positive finite number.", call. = FALSE)
  }
  observed <- is.finite(X)
  signal <- pmax(X, 0)
  score <- signal / (signal + half_saturation)
  score[observed & signal <= 0] <- 0
  score[!observed] <- NA_real_
  dimnames(score) <- dimnames(X)
  attr(score, "score_semantics") <- paste(
    "zero-preserving bounded support from non-negative normalized signal;",
    "not a probability or enzyme capacity"
  )
  score
}
