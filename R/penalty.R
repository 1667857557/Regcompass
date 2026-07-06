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
  if (!is.null(reaction_roles)) {
    rr <- if (is.data.frame(reaction_roles)) reaction_roles else as.data.frame(reaction_roles)
    if (all(c("reaction_id", "role") %in% colnames(rr))) role[intersect(rownames(C), as.character(rr$reaction_id))] <- as.character(rr$role[match(intersect(rownames(C), as.character(rr$reaction_id)), as.character(rr$reaction_id))])
  }
  P_role <- matrix(0, nrow(C), ncol(C), dimnames = dimnames(C))
  for (nm in intersect(names(support_penalty), unique(role))) P_role[role == nm, ] <- as.numeric(support_penalty[[nm]])
  W <- c(expr = 1, confidence = 0.5, missing = 1); W[names(weights)] <- weights
  P <- W["expr"] * P_expr + W["confidence"] * P_conf + W["missing"] * P_missing + P_role
  P <- pmin(pmax(P, 0), penalty_cap); P[!is.finite(P)] <- penalty_cap
  list(penalty = P, components = list(P_expr = P_expr, P_conf = P_conf, P_missing = P_missing, P_role = P_role, C_rel = C, reaction_confidence = F, missing_evidence_flag = missing_flag, role = role))
}
