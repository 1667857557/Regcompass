#' Compute decomposed multiome penalties for microCOMPASS
.rc_compute_multiome_penalty_core <- function(
    C_rel, reaction_confidence, gpr_diagnostics = NULL,
    reaction_roles = NULL,
    weights = c(expr = 1.0, confidence = 0.5, missing = 1.0,
                gpr_missing = 0),
    eps = 1e-6, penalty_cap = 20,
    support_penalty = c(
      exchange = 1.0, demand = 20, sink = 20,
      artificial_support = 20, cofactor_recycle = 0.50, transport = 1.00
    ),
    missing_penalty = 1) {
  if (!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps <= 0 ||
      !is.numeric(penalty_cap) || length(penalty_cap) != 1L ||
      !is.finite(penalty_cap) || penalty_cap <= 0 ||
      !is.numeric(missing_penalty) || length(missing_penalty) != 1L ||
      !is.finite(missing_penalty) || missing_penalty < 0) {
    stop(
      "Penalty constants must be finite and satisfy eps > 0, cap > 0, missing >= 0.",
      call. = FALSE
    )
  }

  C_input <- as.matrix(C_rel)
  F_input <- rc_layer2_confidence_matrix(
    reaction_confidence,
    C_input
  )
  if (is.null(rownames(C_input)) || is.null(colnames(C_input)) ||
      is.null(rownames(F_input)) || is.null(colnames(F_input)) ||
      anyDuplicated(rownames(C_input)) ||
      anyDuplicated(colnames(C_input)) ||
      anyDuplicated(rownames(F_input)) ||
      anyDuplicated(colnames(F_input))) {
    stop(
      "Capacity and confidence matrices require unique reaction and unit IDs.",
      call. = FALSE
    )
  }

  reactions <- union(rownames(C_input), rownames(F_input))
  units <- union(colnames(C_input), colnames(F_input))
  align <- function(matrix_in, fill = NA_real_) {
    output <- matrix(
      fill,
      nrow = length(reactions),
      ncol = length(units),
      dimnames = list(reactions, units)
    )
    common_r <- intersect(reactions, rownames(matrix_in))
    common_u <- intersect(units, colnames(matrix_in))
    output[common_r, common_u] <-
      matrix_in[common_r, common_u, drop = FALSE]
    output
  }

  C <- align(C_input)
  F_original <- align(F_input)
  C[is.finite(C)] <- .rc_clamp01(C[is.finite(C)])
  F_original[is.finite(F_original)] <-
    .rc_clamp01(F_original[is.finite(F_original)])

  missing_expression_flag <- !is.finite(C)
  P_expr <- 1 - C
  P_expr[missing_expression_flag] <- missing_penalty

  finite_regulation <- is.finite(F_original)
  F_effective <- matrix(
    1,
    nrow = nrow(F_original),
    ncol = ncol(F_original),
    dimnames = dimnames(F_original)
  )
  F_effective[finite_regulation] <- pmin(
    1,
    pmax(2 * F_original[finite_regulation], eps)
  )
  P_conf <- -log(pmax(F_effective, eps))
  P_conf[!is.finite(P_conf)] <- 0

  # Missing expression already receives the same maximum inverse-expression
  # penalty as observed zero; do not add a second missingness penalty.
  P_missing <- matrix(
    0,
    nrow = nrow(C),
    ncol = ncol(C),
    dimnames = dimnames(C)
  )

  P_gpr <- matrix(
    0,
    nrow = nrow(C),
    ncol = ncol(C),
    dimnames = dimnames(C)
  )
  gpr_missing_fraction <- stats::setNames(
    rep(0, nrow(C)),
    rownames(C)
  )
  if (!is.null(gpr_diagnostics)) {
    if (!is.data.frame(gpr_diagnostics) ||
        !all(c("reaction_id", "missing_gene_fraction") %in%
             colnames(gpr_diagnostics))) {
      stop(
        "`gpr_diagnostics` must contain reaction_id and missing_gene_fraction.",
        call. = FALSE
      )
    }
    if (anyDuplicated(as.character(gpr_diagnostics$reaction_id))) {
      stop("`gpr_diagnostics$reaction_id` must be unique.",
           call. = FALSE)
    }
    hit <- intersect(
      rownames(C),
      as.character(gpr_diagnostics$reaction_id)
    )
    values <- as.numeric(
      gpr_diagnostics$missing_gene_fraction[
        match(hit, as.character(gpr_diagnostics$reaction_id))
      ]
    )
    values[!is.finite(values)] <- 0
    gpr_missing_fraction[hit] <- pmin(pmax(values, 0), 1)
    P_gpr <- matrix(
      gpr_missing_fraction[rownames(C)],
      nrow = nrow(C),
      ncol = ncol(C),
      dimnames = dimnames(C)
    )
  }

  default_weights <- c(
    expr = 1,
    confidence = 0.5,
    missing = 1,
    gpr_missing = 0
  )
  if (is.null(names(weights)) ||
      any(!names(weights) %in% names(default_weights)) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop(
      "`weights` must be a named non-negative vector using expr, confidence, missing, or gpr_missing.",
      call. = FALSE
    )
  }
  W <- default_weights
  W[names(weights)] <- weights
  P_base <- W[["expr"]] * P_expr +
    W[["confidence"]] * P_conf +
    W[["missing"]] * P_missing +
    W[["gpr_missing"]] * P_gpr

  role <- stats::setNames(rep("internal", nrow(C)), rownames(C))
  role_source <- stats::setNames(
    rep("unknown", nrow(C)),
    rownames(C)
  )
  role_confidence <- stats::setNames(
    rep(NA_character_, nrow(C)),
    rownames(C)
  )
  if (!is.null(reaction_roles)) {
    roles <- if (is.data.frame(reaction_roles)) {
      reaction_roles
    } else {
      as.data.frame(reaction_roles)
    }
    if (!all(c("reaction_id", "role") %in% colnames(roles))) {
      stop(
        "`reaction_roles` must contain reaction_id and role.",
        call. = FALSE
      )
    }
    if (anyDuplicated(as.character(roles$reaction_id))) {
      stop("`reaction_roles$reaction_id` must be unique.",
           call. = FALSE)
    }
    hit <- intersect(rownames(C), as.character(roles$reaction_id))
    index <- match(hit, as.character(roles$reaction_id))
    role[hit] <- as.character(roles$role[index])
    if ("role_source" %in% colnames(roles)) {
      role_source[hit] <- as.character(roles$role_source[index])
    }
    if ("role_confidence" %in% colnames(roles)) {
      role_confidence[hit] <- as.character(roles$role_confidence[index])
    }
  }

  if (is.null(names(support_penalty)) ||
      any(!is.finite(support_penalty)) ||
      any(support_penalty < 0)) {
    stop(
      "`support_penalty` must be a named finite non-negative vector.",
      call. = FALSE
    )
  }

  structural_role <- role %in%
    c("exchange", "demand", "sink", "artificial_support")
  curated_role <- role %in% names(support_penalty) &
    role_source %in% c("curated", "model_high_confidence")
  override <- (structural_role | curated_role) &
    role %in% names(support_penalty)
  support_penalty_used <- stats::setNames(
    rep(NA_real_, nrow(C)),
    rownames(C)
  )
  support_penalty_used[override] <- as.numeric(
    support_penalty[role[override]]
  )

  P <- P_base
  if (any(override)) {
    P[override, ] <- support_penalty_used[override]
  }
  P <- pmin(pmax(P, 0), penalty_cap)
  P[!is.finite(P)] <- penalty_cap

  P_role <- matrix(
    0,
    nrow = nrow(C),
    ncol = ncol(C),
    dimnames = dimnames(C)
  )
  P_role[override, ] <- P[override, , drop = FALSE]

  list(
    penalty = P,
    components = list(
      P_expr = P_expr,
      P_conf = P_conf,
      P_missing = P_missing,
      P_gpr = P_gpr,
      P_role = P_role,
      P_base = P_base,
      C_abs = C,
      C_rel = C,
      reaction_regulatory_support = F_original,
      reaction_confidence = F_original,
      reaction_confidence_effective = F_effective,
      missing_expression_flag = missing_expression_flag,
      missing_regulatory_support_flag = !finite_regulation,
      missing_evidence_flag =
        missing_expression_flag | !finite_regulation,
      gpr_missing_fraction = gpr_missing_fraction,
      role = role,
      role_source = role_source,
      role_confidence = role_confidence,
      role_override_flag = override,
      support_penalty_used = support_penalty_used
    ),
    evidence_policy = paste(
      "COMPASS-like inverse-expression penalty from bounded absolute RNA",
      "support; observed zero, missing expression and no-GPR reactions share",
      "the maximum expression penalty; exchange flux is controlled by one",
      "shared model-bound medium"
    ),
    penalty_formula = paste(
      "w_expr*(1-C_abs) + repression_modifier +",
      "w_gpr*gpr_missing; missing C_abs uses the same maximum",
      "inverse-expression penalty as C_abs=0"
    )
  )
}


rc_compute_multiome_penalty <- function(...) {
  answer <- .rc_compute_multiome_penalty_core(...)
  reaction_ids <- rownames(answer$penalty)
  named_components <- c(
    "gpr_missing_fraction", "role", "role_source", "role_confidence",
    "role_override_flag", "support_penalty_used"
  )
  for (name in named_components) {
    value <- answer$components[[name]]
    if (!is.null(value) && is.null(dim(value)) &&
        length(value) == length(reaction_ids)) {
      names(value) <- reaction_ids
      answer$components[[name]] <- value
    }
  }
  answer$evidence_policy <- "penalty_only"
  answer$evidence_description <- paste(
    "This is not the original COMPASS expression-neighbourhood penalty.",
    "RegCompass uses a COMPASS-like inverse-support expression term;",
    "multiome evidence modifies the LP objective penalty only and does not",
    "directly change stoichiometry or internal reaction bounds."
  )
  answer
}
