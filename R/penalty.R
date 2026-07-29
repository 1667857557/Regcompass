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

# Compute the canonical COMPASS-like cost from multiome reaction expression.
# Regulatory evidence is integrated before GPR aggregation; no independent
# reaction-confidence term is added.
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
    )) {
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
      !is.finite(penalty_cap) || penalty_cap <= 0) {
    stop(
      "`eps` and `penalty_cap` must be finite positive constants.",
      call. = FALSE
    )
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

  observed <- is.finite(E)
  E_effective <- E
  E_effective[observed] <- pmax(E_effective[observed], 0)
  P_expr <- 1 / (1 + log2(1 + E_effective))
  dimnames(P_expr) <- dimnames(E)

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
  finite_penalty <- is.finite(penalty)
  penalty[finite_penalty] <- pmin(
    pmax(penalty[finite_penalty], eps),
    penalty_cap
  )
  availability <- observed
  if (any(override)) {
    availability[override, ] <- TRUE
  }

  list(
    penalty = penalty,
    components = list(
      reaction_expression = E,
      effective_reaction_expression = E_effective,
      P_expr = P_expr,
      role = role,
      role_source = role_source,
      role_override_flag = override,
      penalty_available = availability,
      missing_expression_flag = !observed,
      observed_zero_expression_flag = observed & E_effective <= 0
    ),
    evidence_policy = "penalty_only",
    evidence_policy_detail = paste(
      "unmeasured reaction expression remains unavailable (NA), whereas an",
      "observed zero receives the strictest expression-linked penalty;",
      "fixed costs are used only for exchange/demand/sink/artificial-support reactions"
    ),
    penalty_version = "gene_integrated_multiome_penalty_v2",
    evidence_description = paste(
      "Condition-specific Pando coefficients learned from RNA+ATAC weight",
      "accessibility-only regulatory deviations integrated into gene support",
      "before GPR aggregation; expression-linked reactions use",
      "1/(1+log2(1+reaction_expression)); missing expression remains NA."
    ),
    penalty_formula = "1 / (1 + log2(1 + pmax(E_multiome, 0))); missing E_multiome := NA"
  )
}
