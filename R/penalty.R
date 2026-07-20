.rc_condition_role_vectors <- function(reaction_ids, reaction_roles = NULL) {
  role <- stats::setNames(rep("internal", length(reaction_ids)), reaction_ids)
  role_source <- stats::setNames(rep("unknown", length(reaction_ids)), reaction_ids)
  if (!is.null(reaction_roles)) {
    roles <- as.data.frame(reaction_roles)
    if (!all(c("reaction_id", "role") %in% colnames(roles))) {
      stop("`reaction_roles` must contain reaction_id and role.", call. = FALSE)
    }
    hit <- intersect(reaction_ids, as.character(roles$reaction_id))
    idx <- match(hit, as.character(roles$reaction_id))
    role[hit] <- as.character(roles$role[idx])
    if ("role_source" %in% colnames(roles)) {
      role_source[hit] <- as.character(roles$role_source[idx])
    }
  }
  list(role = role, role_source = role_source)
}

# Compute the v1.7.0 COMPASS-like penalty from multiome reaction expression.
# Regulatory evidence is integrated into gene support before GPR aggregation;
# there is no independent reaction-confidence term in the canonical model.
rc_compute_multiome_penalty <- function(
    reaction_expression,
    reaction_roles = NULL,
    eps = 1e-6,
    penalty_cap = 20,
    support_penalty = c(
      exchange = 1.0,
      demand = 20,
      sink = 20,
      artificial_support = 20
    ),
    missing_penalty = 1) {
  E <- as.matrix(reaction_expression)
  if (!is.numeric(E) || is.null(rownames(E)) || is.null(colnames(E)) ||
      anyDuplicated(rownames(E)) || anyDuplicated(colnames(E))) {
    stop(
      "Reaction expression requires a numeric matrix with unique dimnames.",
      call. = FALSE
    )
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0 ||
      !is.numeric(penalty_cap) || length(penalty_cap) != 1L ||
      !is.finite(penalty_cap) || penalty_cap <= 0 ||
      !is.numeric(missing_penalty) || length(missing_penalty) != 1L ||
      !is.finite(missing_penalty) || missing_penalty < 0) {
    stop("Penalty constants are invalid.", call. = FALSE)
  }
  required_structural_roles <- c(
    "exchange", "demand", "sink", "artificial_support"
  )
  if (!is.numeric(support_penalty) || is.null(names(support_penalty)) ||
      anyDuplicated(names(support_penalty)) ||
      any(!required_structural_roles %in% names(support_penalty)) ||
      any(!is.finite(support_penalty)) || any(support_penalty < 0)) {
    stop(
      "`support_penalty` must provide finite non-negative costs for all structural roles.",
      call. = FALSE
    )
  }

  finite <- is.finite(E)
  E_nonnegative <- pmax(E, 0)
  P_expr <- matrix(
    missing_penalty,
    nrow = nrow(E),
    ncol = ncol(E),
    dimnames = dimnames(E)
  )
  P_expr[finite] <- 1 / (1 + log2(1 + E_nonnegative[finite]))

  roles <- .rc_condition_role_vectors(rownames(E), reaction_roles)
  role <- roles$role
  role_source <- roles$role_source
  override <- stats::setNames(
    as.logical(role %in% required_structural_roles),
    rownames(E)
  )

  penalty <- P_expr
  if (any(override)) {
    penalty[override, ] <- as.numeric(support_penalty[role[override]])
  }
  penalty <- pmin(pmax(penalty, eps), penalty_cap)
  penalty[!is.finite(penalty)] <- penalty_cap

  list(
    penalty = penalty,
    components = list(
      reaction_expression = E,
      P_expr = P_expr,
      role = role,
      role_source = role_source,
      role_override_flag = override,
      missing_expression_flag = !finite
    ),
    evidence_policy = "penalty_only",
    evidence_policy_detail = paste(
      "single integrated expression penalty with fixed costs only for",
      "exchange/demand/sink/artificial-support reactions"
    ),
    penalty_version = "v1.7.0_gene_integrated_multiome_penalty",
    evidence_description = paste(
      "Condition-specific Pando coefficients learned from RNA+ATAC weight",
      "accessibility-only regulatory deviations integrated into gene support",
      "before GPR aggregation; expression-linked reactions use",
      "1/(1+log2(1+reaction_expression))."
    ),
    penalty_formula = "1 / (1 + log2(1 + E_multiome))"
  )
}
