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
      !is.finite(missing_penalty) || missing_penalty != 1) {
    stop(
      "Penalty constants are invalid; `missing_penalty` must remain 1 so unmeasured expression is treated as zero expression.",
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
  E_effective[!observed] <- 0
  E_effective <- pmax(E_effective, 0)
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
  penalty <- pmin(pmax(penalty, eps), penalty_cap)
  penalty[!is.finite(penalty)] <- penalty_cap

  list(
    penalty = penalty,
    components = list(
      reaction_expression = E,
      effective_reaction_expression = E_effective,
      P_expr = P_expr,
      role = role,
      role_source = role_source,
      role_override_flag = override,
      missing_expression_flag = !observed,
      zero_or_missing_expression_flag = !observed | E_effective <= 0
    ),
    evidence_policy = "penalty_only",
    evidence_policy_detail = paste(
      "unmeasured and explicit zero reaction expression are both treated as",
      "zero support and receive the strictest expression-linked penalty;",
      "fixed costs are used only for exchange/demand/sink/artificial-support reactions"
    ),
    penalty_version = "gene_integrated_multiome_penalty_v1",
    evidence_description = paste(
      "Condition-specific Pando coefficients learned from RNA+ATAC weight",
      "accessibility-only regulatory deviations integrated into gene support",
      "before GPR aggregation; expression-linked reactions use",
      "1/(1+log2(1+reaction_expression)), with missing expression zero-filled."
    ),
    penalty_formula = "1 / (1 + log2(1 + pmax(E_multiome, 0))); missing E_multiome := 0"
  )
}
