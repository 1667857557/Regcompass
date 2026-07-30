# Exact COMPASS-style GPR aggregation and missing-expression penalties.

rc_and_capacity <- function(scores, method = c("min", "median", "mean")) {
  method <- match.arg(method)
  scores <- as.numeric(scores)
  if (!length(scores)) return(NA_real_)
  if (identical(method, "min")) {
    if (any(!is.finite(scores))) return(NA_real_)
    return(min(scores))
  }
  scores[!is.finite(scores)] <- 0
  if (identical(method, "median")) stats::median(scores) else mean(scores)
}

rc_or_capacity <- function(
    and_capacities,
    method = c("sum", "max", "sum_sqrtK", "prob_or")) {
  method <- match.arg(method)
  x <- as.numeric(and_capacities)
  finite <- is.finite(x)
  if (!any(finite)) return(NA_real_)
  x <- x[finite]
  switch(
    method,
    sum = sum(x),
    max = max(x),
    sum_sqrtK = sum(x) / sqrt(length(x)),
    prob_or = 1 - prod(1 - pmin(pmax(x, 0), 1))
  )
}

rc_reaction_capacity_one <- function(
    parsed_gpr, gene_score_vec,
    and_method = c("min", "median", "mean"),
    or_method = c("sum", "max", "sum_sqrtK", "prob_or")) {
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)
  and_caps <- vapply(parsed_gpr, function(and_group) {
    genes <- unique(tolower(and_group))
    if (!length(genes)) return(NA_real_)
    vals <- gene_score_vec[genes]
    rc_and_capacity(vals, method = and_method)
  }, numeric(1))
  rc_or_capacity(and_caps, method = or_method)
}

.rc_reaction_capacity_previous <- rc_reaction_capacity
rc_reaction_capacity <- function(
    gpr_list,
    gene_score,
    promiscuity_mode = c("none", "sqrt", "linear"),
    and_method = c("min", "median", "mean"),
    or_method = c("sum", "max", "sum_sqrtK", "prob_or"),
    BPPARAM = NULL) {
  promiscuity_mode <- match.arg(promiscuity_mode)
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)
  .rc_reaction_capacity_previous(
    gpr_list = gpr_list,
    gene_score = gene_score,
    promiscuity_mode = promiscuity_mode,
    and_method = and_method,
    or_method = or_method,
    BPPARAM = BPPARAM
  )
}

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
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) ||
      eps <= 0 || !is.numeric(penalty_cap) || length(penalty_cap) != 1L ||
      !is.finite(penalty_cap) || penalty_cap <= 0) {
    stop("`eps` and `penalty_cap` must be finite positive constants.",
         call. = FALSE)
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

  expression_observed <- is.finite(E)
  E_effective <- E
  E_effective[!expression_observed] <- 0
  E_effective <- pmax(E_effective, 0)
  P_expr <- 1 / (1 + log2(1 + E_effective))
  dimnames(P_expr) <- dimnames(E)

  roles <- .rc_condition_role_vectors(rownames(E), reaction_roles)
  role <- roles$role
  role_source <- roles$role_source
  override <- stats::setNames(
    as.logical(role %in% required_structural_roles), rownames(E)
  )
  penalty <- P_expr
  if (any(override)) {
    penalty[override, ] <- as.numeric(support_penalty[role[override]])
  }
  penalty <- pmin(pmax(penalty, eps), penalty_cap)
  dimnames(penalty) <- dimnames(E)
  penalty_available <- matrix(
    TRUE, nrow(E), ncol(E), dimnames = dimnames(E)
  )
  missing_expression_flag <- !expression_observed
  maximum_expression_penalty_flag <-
    missing_expression_flag & !matrix(override, nrow(E), ncol(E))

  list(
    penalty = penalty,
    components = list(
      reaction_expression = E,
      effective_reaction_expression = E_effective,
      P_expr = P_expr,
      role = role,
      role_source = role_source,
      role_override_flag = override,
      penalty_available = penalty_available,
      expression_observed = expression_observed,
      missing_expression_flag = missing_expression_flag,
      missing_expression_imputed_zero = missing_expression_flag,
      maximum_expression_penalty_flag = maximum_expression_penalty_flag,
      observed_zero_expression_flag = expression_observed & E_effective <= 0
    ),
    evidence_policy = "compass_missing_expression_max_penalty",
    evidence_policy_detail = paste(
      "COMPASS-compatible reaction expression: complete AND branches are",
      "combined by the selected AND function, OR branches are summed while",
      "ignoring unavailable branches, and an unavailable final reaction",
      "expression is set to zero before penalty conversion"
    ),
    penalty_version = "compass_gpr_missing_zero_penalty_v3",
    evidence_description = paste(
      "Expression-linked reactions use 1/(1+log2(1+E)); missing E is assigned",
      "E=0 and therefore the maximum expression-linked penalty of 1. Fixed",
      "costs remain restricted to structural support reactions."
    ),
    penalty_formula = paste(
      "P = 1 / (1 + log2(1 + pmax(E, 0)));",
      "nonfinite E := 0, hence P := 1"
    )
  )
}
