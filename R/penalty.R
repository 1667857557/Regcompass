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

rc_compute_multiome_penalty <- function(
    reaction_expression,
    reaction_roles = NULL,
    eps = 1e-6,
    penalty_cap = 1,
    support_penalty = c(
      exchange = 1.0, demand = 1.0, sink = 1.0, artificial_support = 1.0
    )) {
  E <- as.matrix(reaction_expression)
  if (!is.numeric(E) || is.null(rownames(E)) || is.null(colnames(E)) ||
      anyDuplicated(rownames(E)) || anyDuplicated(colnames(E))) {
    stop("Reaction expression requires a numeric matrix with unique dimnames.",
         call. = FALSE)
  }
  required_roles <- c("exchange", "demand", "sink", "artificial_support")
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0 ||
      !is.numeric(penalty_cap) || length(penalty_cap) != 1L ||
      !is.finite(penalty_cap) || penalty_cap <= 0 ||
      !is.numeric(support_penalty) || is.null(names(support_penalty)) ||
      any(!required_roles %in% names(support_penalty)) ||
      any(!is.finite(support_penalty)) || any(support_penalty < 0)) {
    stop("Penalty controls are invalid.", call. = FALSE)
  }
  expression_observed <- is.finite(E)
  E_effective <- E
  E_effective[!expression_observed] <- 0
  E_effective <- pmax(E_effective, 0)
  P_expr <- 1 / (1 + log2(1 + E_effective))
  dimnames(P_expr) <- dimnames(E)
  roles <- .rc_condition_role_vectors(rownames(E), reaction_roles)
  override <- stats::setNames(roles$role %in% required_roles, rownames(E))
  penalty <- P_expr
  if (any(override)) {
    penalty[override, ] <- as.numeric(support_penalty[roles$role[override]])
  }
  penalty <- pmin(pmax(penalty, eps), penalty_cap)
  dimnames(penalty) <- dimnames(E)
  missing <- !expression_observed
  list(
    penalty = penalty,
    components = list(
      reaction_expression = E,
      effective_reaction_expression = E_effective,
      P_expr = P_expr,
      role = roles$role,
      role_source = roles$role_source,
      role_override_flag = override,
      penalty_available = matrix(TRUE, nrow(E), ncol(E), dimnames = dimnames(E)),
      expression_observed = expression_observed,
      missing_expression_flag = missing,
      missing_expression_imputed_zero = missing,
      maximum_expression_penalty_flag = missing &
        !matrix(override, nrow(E), ncol(E)),
      observed_zero_expression_flag = expression_observed & E_effective <= 0
    ),
    evidence_policy = "penalty_only",
    evidence_policy_detail = paste(
      "GPR expression uses complete AND branches and additive OR branches;",
      "unavailable final reaction expression is set to zero before penalty conversion"
    ),
    missing_expression_policy = "compass_missing_expression_max_penalty",
    structural_reaction_policy =
      "compass_maximum_expression_penalty_for_all_structural_roles",
    penalty_version = "compass_gpr_missing_zero_penalty_v4",
    evidence_description = paste(
      "RNA and Pando regulatory evidence are integrated into gene support before",
      "GPR aggregation. All reactions use the COMPASS cost scale",
      "1/(1+log2(1+E)); missing E and structural reactions receive cost 1."
    ),
    penalty_formula = paste(
      "P=1/(1+log2(1+pmax(E,0))); nonfinite E:=0, hence P:=1;",
      "exchange, demand, sink and artificial_support use P:=1"
    )
  )
}
