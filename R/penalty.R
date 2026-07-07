#' Compute decomposed multiome penalties for microCOMPASS
#' @export
rc_compute_multiome_penalty <- function(C_rel, reaction_confidence, gpr_diagnostics = NULL, reaction_roles = NULL,
                                        weights = c(expr = 1.0, confidence = 0.5, missing = 1.0), eps = 1e-6,
                                        penalty_cap = 20,
                                        support_penalty = c(exchange = 0.05, demand = 0.10, sink = 0.10, artificial_support = 0.05, cofactor_recycle = 0.50, transport = 1.00),
                                        missing_penalty = 5) {
  C <- as.matrix(C_rel); F <- rc_layer2_confidence_matrix(reaction_confidence, C)
  rx <- union(rownames(C), rownames(F)); units <- union(colnames(C), colnames(F))
  align <- function(M, fill) { out <- matrix(fill, length(rx), length(units), dimnames = list(rx, units)); out[intersect(rx, rownames(M)), intersect(units, colnames(M))] <- M[intersect(rx, rownames(M)), intersect(units, colnames(M)), drop = FALSE]; out }
  C <- align(C, NA_real_); F <- align(F, NA_real_)
  P_expr <- -log(pmax(C, eps)); P_expr[!is.finite(P_expr)] <- penalty_cap
  P_conf <- -log(pmax(F, eps)); P_conf[!is.finite(P_conf)] <- 0
  missing_flag <- is.na(C) | is.na(F)
  P_missing <- matrix(0, nrow(C), ncol(C), dimnames = dimnames(C)); P_missing[missing_flag] <- missing_penalty
  role <- rep("internal", nrow(C)); names(role) <- rownames(C)
  role_source <- rep("unknown", nrow(C)); names(role_source) <- rownames(C)
  role_confidence <- rep(NA_character_, nrow(C)); names(role_confidence) <- rownames(C)
  if (!is.null(reaction_roles)) {
    rr <- if (is.data.frame(reaction_roles)) reaction_roles else as.data.frame(reaction_roles)
    if (all(c("reaction_id", "role") %in% colnames(rr))) {
      hit <- intersect(rownames(C), as.character(rr$reaction_id))
      m <- match(hit, as.character(rr$reaction_id))
      role[hit] <- as.character(rr$role[m])
      if ("role_source" %in% colnames(rr)) role_source[hit] <- as.character(rr$role_source[m])
      if ("role_confidence" %in% colnames(rr)) role_confidence[hit] <- as.character(rr$role_confidence[m])
    }
  }
  P_role <- matrix(0, nrow(C), ncol(C), dimnames = dimnames(C))
  role_override_flag <- role %in% c("exchange", "demand", "sink", "artificial_support") &
    role_source %in% c("curated", "model_high_confidence")
  support_penalty_used <- rep(NA_real_, nrow(C)); names(support_penalty_used) <- rownames(C)
  for (nm in intersect(names(support_penalty), unique(role[role_override_flag]))) {
    support_penalty_used[role_override_flag & role == nm] <- as.numeric(support_penalty[[nm]])
  }
  transport_evidence_flag <- role == "transport" & is.finite(rowMeans(C, na.rm = TRUE))
  W <- c(expr = 1, confidence = 0.5, missing = 1); W[names(weights)] <- weights
  P_base <- W["expr"] * P_expr + W["confidence"] * P_conf + W["missing"] * P_missing
  P <- P_base
  if (any(role_override_flag)) {
    P[role_override_flag, ] <- matrix(support_penalty_used[role_override_flag],
                                      nrow = sum(role_override_flag), ncol = ncol(P),
                                      dimnames = list(names(role)[role_override_flag], colnames(P)))
  }
  P_role[role_override_flag, ] <- P[role_override_flag, , drop = FALSE]
  P <- pmin(pmax(P, 0), penalty_cap); P[!is.finite(P)] <- penalty_cap
  list(penalty = P, components = list(P_expr = P_expr, P_conf = P_conf, P_missing = P_missing, P_role = P_role, P_base = P_base, C_rel = C, reaction_confidence = F, missing_evidence_flag = missing_flag, role = role, role_source = role_source, role_confidence = role_confidence, role_override_flag = role_override_flag, support_penalty_used = support_penalty_used, transport_evidence_flag = transport_evidence_flag))
}
