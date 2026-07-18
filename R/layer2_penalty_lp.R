rc_compass_score_from_penalty <- function(P, feasible,
                                          method = c("ecdf", "mad_sigmoid"),
                                          variation_tolerance = 1e-8) {
  method <- match.arg(method)
  P <- as.matrix(P)
  feasible <- as.matrix(feasible)
  valid_dimnames <- function(x) {
    !is.null(rownames(x)) && !is.null(colnames(x)) &&
      !anyNA(rownames(x)) && !anyNA(colnames(x)) &&
      all(nzchar(rownames(x))) && all(nzchar(colnames(x))) &&
      !anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x))
  }
  if (!is.numeric(P) || !is.logical(feasible) ||
      !valid_dimnames(P) || !valid_dimnames(feasible)) {
    stop(
      "`P` must be numeric and `feasible` logical; both require unique ",
      "non-empty target and unit IDs.",
      call. = FALSE
    )
  }
  if (!setequal(rownames(P), rownames(feasible)) ||
      !setequal(colnames(P), colnames(feasible))) {
    stop("`P` and `feasible` must contain identical target and unit IDs.",
         call. = FALSE)
  }
  feasible <- feasible[rownames(P), colnames(P), drop = FALSE]
  if (anyNA(feasible)) {
    stop("`feasible` cannot contain missing values.", call. = FALSE)
  }
  if (!is.numeric(variation_tolerance) ||
      length(variation_tolerance) != 1L ||
      !is.finite(variation_tolerance) || variation_tolerance < 0) {
    stop("`variation_tolerance` must be one finite non-negative number.",
         call. = FALSE)
  }
  score <- matrix(NA_real_, nrow(P), ncol(P), dimnames = dimnames(P))
  noninformative <- logical(nrow(P))
  for (i in seq_len(nrow(P))) {
    index <- feasible[i, ] & is.finite(P[i, ])
    x <- P[i, index]
    if (length(x) < 2L || diff(range(x)) <= variation_tolerance) {
      noninformative[[i]] <- TRUE
      next
    }
    if (identical(method, "ecdf")) {
      ranks <- rank(x, ties.method = "average")
      score[i, index] <- 1 - (ranks - 1) / (length(x) - 1)
    } else {
      center <- stats::median(x)
      scale <- max(stats::mad(x, constant = 1.4826),
                   stats::IQR(x) / 1.349, variation_tolerance)
      score[i, index] <- rc_sigmoid((center - x) / scale)
    }
  }
  attr(score, "score_semantics") <- if (identical(method, "ecdf")) {
    "within_target_relative_penalty_rank_not_probability"
  } else {
    "within_target_robust_penalty_transform_not_probability"
  }
  attr(score, "noninformative_target") <- stats::setNames(
    noninformative, rownames(P)
  )
  score
}

rc_layer2_unit_matrices <- function(
    layer1, unit, sample_col, celltype_col, condition_col) {
  unit <- match.arg(unit, c("sample_celltype", "metacell"))
  C <- as.matrix(layer1$C_rel)
  Conf <- rc_layer2_confidence_matrix(layer1$reaction_confidence, C)
  valid_dimnames <- function(x) {
    !is.null(rownames(x)) && !is.null(colnames(x)) &&
      !anyNA(rownames(x)) && !anyNA(colnames(x)) &&
      all(nzchar(rownames(x))) && all(nzchar(colnames(x))) &&
      !anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x))
  }
  if (!valid_dimnames(C) || !valid_dimnames(Conf)) {
    stop(
      "Layer 1 capacity and confidence matrices require unique non-empty ",
      "reaction and unit IDs.",
      call. = FALSE
    )
  }
  if (!setequal(colnames(C), colnames(Conf))) {
    stop("Layer 1 capacity and confidence matrices contain different units.",
         call. = FALSE)
  }
  Conf <- Conf[, colnames(C), drop = FALSE]
  common <- intersect(rownames(C), rownames(Conf))
  if (!length(common)) {
    stop("Layer 1 capacity and confidence matrices share no reactions.",
         call. = FALSE)
  }
  C <- C[common, , drop = FALSE]
  Conf <- Conf[common, , drop = FALSE]
  if (identical(unit, "metacell")) {
    return(list(
      C_rel = C,
      reaction_confidence = Conf,
      unit_meta = layer1$unit_meta,
      summary = "metacell-level Layer1 matrices"
    ))
  }
  if (is.null(layer1$unit_meta) || !is.data.frame(layer1$unit_meta)) {
    stop("sample_celltype units require data-frame `layer1$unit_meta`.",
         call. = FALSE)
  }
  pm <- layer1$unit_meta
  required <- c("pool_id", sample_col, celltype_col)
  missing <- setdiff(required, colnames(pm))
  if (length(missing)) {
    stop("unit_meta missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  pool_ids <- trimws(as.character(pm$pool_id))
  if (anyNA(pool_ids) || any(!nzchar(pool_ids)) || anyDuplicated(pool_ids)) {
    stop("`unit_meta$pool_id` must contain unique non-empty unit IDs.",
         call. = FALSE)
  }
  if (!setequal(pool_ids, colnames(C))) {
    stop("`unit_meta$pool_id` does not exactly match Layer 1 matrix units.",
         call. = FALSE)
  }
  pm$pool_id <- pool_ids
  pm <- pm[match(colnames(C), pm$pool_id), , drop = FALSE]
  group_cols <- c(
    sample_col,
    if (condition_col %in% colnames(pm)) condition_col else NULL,
    celltype_col
  )
  invalid_group <- vapply(pm[, group_cols, drop = FALSE], function(x) {
    value <- trimws(as.character(x))
    anyNA(value) || any(!nzchar(value))
  }, logical(1))
  if (any(invalid_group)) {
    stop("Layer 2 grouping metadata contain missing or empty values: ",
         paste(group_cols[invalid_group], collapse = ", "), call. = FALSE)
  }
  gid <- interaction(
    pm[, group_cols, drop = FALSE],
    sep = "|", drop = TRUE, lex.order = TRUE
  )
  agg <- function(M) {
    do.call(cbind, lapply(
      split(pm$pool_id, gid),
      function(cols) matrixStats::rowMedians(M[, cols, drop = FALSE], na.rm = TRUE)
    ))
  }
  outC <- agg(C)
  outF <- agg(Conf)
  unit_meta <- unique(data.frame(
    unit_id = as.character(gid),
    pm[, group_cols, drop = FALSE],
    stringsAsFactors = FALSE
  ))
  unit_meta <- unit_meta[match(colnames(outC), unit_meta$unit_id), , drop = FALSE]
  rownames(unit_meta) <- NULL
  list(
    C_rel = outC,
    reaction_confidence = outF,
    unit_meta = unit_meta,
    summary = "Layer1 aggregated to sample x celltype medians"
  )
}

rc_align_layer2_evidence <- function(M, rxns, fill = NA_real_) {
  M <- as.matrix(M)
  out <- matrix(fill, nrow = length(rxns), ncol = ncol(M), dimnames = list(rxns, colnames(M)))
  common <- intersect(rxns, rownames(M))
  out[common, ] <- M[common, , drop = FALSE]
  out
}

rc_layer2_confidence_matrix <- function(confidence, C_rel) {
  C_rel <- as.matrix(C_rel)
  if (is.data.frame(confidence) && all(c("reaction_id", "pool_id", "reaction_confidence") %in% colnames(confidence))) {
    out <- matrix(NA_real_, nrow = nrow(C_rel), ncol = ncol(C_rel), dimnames = dimnames(C_rel))
    rid <- as.character(confidence$reaction_id)
    pid <- as.character(confidence$pool_id)
    ok <- rid %in% rownames(out) & pid %in% colnames(out)
    if (any(ok)) {
      idx <- cbind(match(rid[ok], rownames(out)), match(pid[ok], colnames(out)))
      out[idx] <- as.numeric(confidence$reaction_confidence[ok])
    }
    return(out)
  }
  M <- as.matrix(confidence)
  if (is.null(rownames(M)) || is.null(colnames(M))) stop("`reaction_confidence` must be a matrix-like object with dimnames or a long data frame.", call. = FALSE)
  M
}
