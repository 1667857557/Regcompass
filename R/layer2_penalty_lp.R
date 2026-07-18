.rc_layer2_penalty_engine <- function(C_rel, Conf, epsilon = 1e-6, epsilon_C = 1e-3, epsilon_Conf = 1e-3,
                              penalty_cap = 20, support_reactions = character(), support_penalty = 0.05) {
  C_raw <- as.matrix(C_rel)
  Conf_raw <- as.matrix(Conf)
  valid_evidence_matrix <- function(x) {
    is.numeric(x) && !is.null(rownames(x)) && !is.null(colnames(x)) &&
      !anyNA(rownames(x)) && !anyNA(colnames(x)) &&
      all(nzchar(rownames(x))) && all(nzchar(colnames(x))) &&
      !anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x))
  }
  if (!valid_evidence_matrix(C_raw) || !valid_evidence_matrix(Conf_raw)) {
    stop(
      "Capacity and confidence must be numeric matrices with unique, ",
      "non-empty reaction and unit IDs.",
      call. = FALSE
    )
  }
  if (!setequal(rownames(C_raw), rownames(Conf_raw)) ||
      !setequal(colnames(C_raw), colnames(Conf_raw))) {
    stop(
      "Capacity and confidence matrices must contain identical reaction and unit IDs.",
      call. = FALSE
    )
  }
  Conf_raw <- Conf_raw[rownames(C_raw), colnames(C_raw), drop = FALSE]
  constants <- c(
    epsilon = epsilon,
    epsilon_C = epsilon_C,
    epsilon_Conf = epsilon_Conf,
    penalty_cap = penalty_cap
  )
  scalar_constants <- vapply(
    list(epsilon, epsilon_C, epsilon_Conf, penalty_cap),
    function(value) is.numeric(value) && length(value) == 1L,
    logical(1)
  )
  if (!all(scalar_constants) || any(!is.finite(constants)) ||
      any(constants <= 0) ||
      !is.numeric(support_penalty) || !length(support_penalty) ||
      any(!is.finite(support_penalty)) || any(support_penalty < 0)) {
    stop(
      "Epsilon and cap values must be positive finite numbers; support ",
      "penalties must be finite and non-negative.",
      call. = FALSE
    )
  }
  outside_unit_interval <- function(x) {
    observed <- x[is.finite(x)]
    any(observed < 0 | observed > 1)
  }
  if (outside_unit_interval(C_raw) || outside_unit_interval(Conf_raw)) {
    stop(
      "Observed capacity and confidence values must lie in [0, 1].",
      call. = FALSE
    )
  }
  Cc <- ifelse(is.na(C_raw), epsilon_C, pmax(C_raw, epsilon_C))
  Fc <- ifelse(is.na(Conf_raw), epsilon_Conf, pmax(Conf_raw, epsilon_Conf))
  E <- Cc * Fc
  P <- pmax(0, pmin(-log(E + epsilon), penalty_cap))
  P[!is.finite(P)] <- penalty_cap
  support_reaction_ids <- character()
  support_penalty_by_reaction <- numeric()
  if (length(support_reactions)) {
    if (is.null(names(support_reactions))) {
      support_reaction_ids <- intersect(as.character(support_reactions), rownames(P))
      support_penalty_by_reaction <- stats::setNames(rep(rc_layer2_support_penalty_for_type("support", support_penalty), length(support_reaction_ids)), support_reaction_ids)
    } else {
      support_reaction_ids <- intersect(names(support_reactions), rownames(P))
      support_penalty_by_reaction <- stats::setNames(as.numeric(support_reactions[support_reaction_ids]), support_reaction_ids)
      if (any(!is.finite(support_penalty_by_reaction)) ||
          any(support_penalty_by_reaction < 0)) {
        stop("Named support penalties must be finite and non-negative.",
             call. = FALSE)
      }
    }
    if (length(support_penalty_by_reaction)) P[names(support_penalty_by_reaction), ] <- support_penalty_by_reaction
  }
  missing_evidence <- is.na(C_raw) | is.na(Conf_raw)
  list(penalty = P, components = list(C_rel = C_raw, reaction_confidence = Conf_raw, evidence_product = E,
                                      missing_evidence = missing_evidence, support_reaction = rownames(P) %in% support_reaction_ids,
                                      support_penalty_by_reaction = support_penalty_by_reaction,
                                      penalty = P, epsilon_C = epsilon_C, epsilon_Conf = epsilon_Conf,
                                      penalty_cap = penalty_cap, support_penalty = support_penalty))
}
rc_layer2_penalty <- function(C_rel, Conf, epsilon = 1e-6,
                              epsilon_C = 1e-3, epsilon_Conf = 1e-3,
                              penalty_cap = 20,
                              support_reactions = character(),
                              support_penalty = 0.05,
                              allow_structural_support_override = FALSE) {
  if (length(support_reactions) && !isTRUE(allow_structural_support_override)) {
    stop(
      "FASTCORE/support membership is structural, not biological evidence. Set `allow_structural_support_override = TRUE` only for sensitivity analysis.",
      call. = FALSE
    )
  }
  .rc_layer2_penalty_engine(
    C_rel, Conf,
    epsilon = epsilon,
    epsilon_C = epsilon_C,
    epsilon_Conf = epsilon_Conf,
    penalty_cap = penalty_cap,
    support_reactions = support_reactions,
    support_penalty = support_penalty
  )
}

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

rc_layer2_support_penalties <- function(gem, rxns, C_rel = NULL, Conf = NULL,
                                        support_classes = c("exchange", "demand", "sink", "support"),
                                        transport_penalty_mode = c("normal", "support", "reduced"),
                                        support_penalty = c(exchange = 0.05, demand = 0.1, sink = 0.1, support = 0.05),
                                        transport_reduced_penalty = 1) {
  transport_penalty_mode <- match.arg(transport_penalty_mode)
  meta <- gem$reaction_meta
  if (is.null(meta) || !is.data.frame(meta) || !"reaction_id" %in% colnames(meta)) return(numeric())
  type <- rc_layer2_reaction_type(meta)
  names(type) <- as.character(meta$reaction_id)
  type <- type[intersect(rxns, names(type))]
  out <- numeric()
  support_type <- type[type %in% support_classes]
  if (length(support_type)) {
    out <- c(out, stats::setNames(vapply(support_type, rc_layer2_support_penalty_for_type, numeric(1), support_penalty), names(support_type)))
  }
  transport <- type[type == "transport"]
  if (length(transport) && identical(transport_penalty_mode, "support")) {
    out <- c(out, stats::setNames(rep(rc_layer2_support_penalty_for_type("support", support_penalty), length(transport)), names(transport)))
  } else if (length(transport) && identical(transport_penalty_mode, "reduced")) {
    no_gpr <- !rc_layer2_has_gpr(meta[match(names(transport), as.character(meta$reaction_id)), , drop = FALSE], C_rel, Conf)
    if (any(no_gpr)) out <- c(out, stats::setNames(rep(transport_reduced_penalty, sum(no_gpr)), names(transport)[no_gpr]))
  }
  out[!duplicated(names(out))]
}

rc_layer2_reaction_type <- function(meta) {
  text_cols <- intersect(c("role", "type", "reaction_type", "category", "subsystem", "name", "reaction_name"), colnames(meta))
  txt <- if (length(text_cols)) apply(meta[, text_cols, drop = FALSE], 1, paste, collapse = " ") else rep("", nrow(meta))
  txt <- tolower(txt)
  out <- rep("other", length(txt))
  out[out == "other" & grepl("exchange", txt)] <- "exchange"
  out[out == "other" & grepl("demand", txt)] <- "demand"
  out[out == "other" & grepl("sink", txt)] <- "sink"
  out[out == "other" & grepl("support|artificial", txt)] <- "support"
  out[out == "other" & grepl("transport", txt)] <- "transport"
  out
}

rc_layer2_support_penalty_for_type <- function(type, support_penalty) {
  if (length(support_penalty) == 1L || is.null(names(support_penalty))) return(as.numeric(support_penalty)[1])
  if (type %in% names(support_penalty)) return(as.numeric(support_penalty[[type]]))
  if ("support" %in% names(support_penalty)) return(as.numeric(support_penalty[["support"]]))
  as.numeric(support_penalty[[1]])
}

rc_layer2_has_gpr <- function(meta, C_rel = NULL, Conf = NULL) {
  if (!is.data.frame(meta) || !"reaction_id" %in% colnames(meta)) {
    stop("`meta` must contain `reaction_id`.", call. = FALSE)
  }
  gpr_columns <- intersect(
    c("gpr", "grRule", "gene_reaction_rule",
      "gene_reaction_rule_string", "genes"),
    colnames(meta)
  )
  has_metadata_gpr <- if (length(gpr_columns)) {
    apply(meta[, gpr_columns, drop = FALSE], 1, function(value) {
      any(!is.na(value) & nzchar(trimws(as.character(value))))
    })
  } else {
    rep(FALSE, nrow(meta))
  }

  has_evidence <- rep(FALSE, nrow(meta))
  if (!is.null(C_rel) && !is.null(Conf)) {
    C <- as.matrix(C_rel)
    F <- as.matrix(Conf)
    if (!is.null(rownames(C)) && !is.null(rownames(F))) {
      reaction_id <- as.character(meta$reaction_id)
      present <- reaction_id %in% rownames(C) & reaction_id %in% rownames(F)
      if (any(present)) {
        selected <- reaction_id[present]
        has_evidence[present] <- rowSums(
          is.finite(C[selected, , drop = FALSE]) |
            is.finite(F[selected, , drop = FALSE])
        ) > 0
      }
    }
  }
  has_metadata_gpr | has_evidence
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
