.rc_currency_ids <- function(gem, currency_metabolites = NULL) {
  if (!is.null(currency_metabolites)) return(unique(as.character(currency_metabolites)))
  if (!is.null(gem$metabolite_meta) && "is_currency" %in% colnames(gem$metabolite_meta)) return(as.character(gem$metabolite_meta$metabolite_id[gem$metabolite_meta$is_currency %in% TRUE]))
  c("h", "h2o", "atp", "adp", "amp", "pi", "ppi", "nad", "nadh", "nadp", "nadph", "fad", "fadh2", "coa", "co2", "o2")
}
#' Build a target reaction-specific micro-GEM
#' @export
rc_build_target_microgem <- function(gem, target_reaction, medium_table = NULL, condition = NULL, k_hop = 2,
                                     include_same_subsystem = TRUE, include_transport = TRUE, include_exchange = TRUE,
                                     include_demand_sink = TRUE, include_cofactor_modules = TRUE,
                                     currency_metabolites = NULL, max_reactions = 500, strict_closure = FALSE) {
  if (!target_reaction %in% colnames(gem$S)) stop("`target_reaction` is not present in `gem$S`.", call. = FALSE)
  gem <- rc_annotate_reaction_roles(gem)
  gv <- rc_validate_gem(gem); S <- gv$S; meta <- gem$reaction_meta[match(gv$reactions, gem$reaction_meta$reaction_id), , drop = FALSE]
  keep <- target_reaction
  if (include_same_subsystem && "subsystem" %in% colnames(meta)) {
    ss <- meta$subsystem[match(target_reaction, meta$reaction_id)]; if (!is.na(ss) && nzchar(ss)) keep <- union(keep, as.character(meta$reaction_id[meta$subsystem == ss]))
  }
  currency <- .rc_currency_ids(gem, currency_metabolites)
  for (d in seq_len(k_hop)) {
    mets <- rownames(S)[Matrix::rowSums(abs(S[, keep, drop = FALSE]) > 0) > 0]
    mets <- setdiff(mets, currency)
    nbr <- colnames(S)[Matrix::colSums(abs(S[mets, , drop = FALSE]) > 0) > 0]
    keep <- union(keep, nbr); if (length(keep) >= max_reactions) { keep <- utils::head(keep, max_reactions); break }
  }
  role <- stats::setNames(as.character(meta$role), meta$reaction_id)
  boundary_mets <- rownames(S)[Matrix::rowSums(abs(S[, keep, drop = FALSE]) > 0) == 1]
  support_roles <- c(if (include_transport) "transport", if (include_exchange) "exchange", if (include_demand_sink) c("demand", "sink"), if (include_cofactor_modules) "cofactor_recycle")
  support <- names(role)[role %in% support_roles & Matrix::colSums(abs(S[boundary_mets, , drop = FALSE]) > 0) > 0]
  keep <- utils::head(unique(c(keep, support)), max_reactions)
  sub <- gem; sub$S <- S[, keep, drop = FALSE]; sub$lb <- gv$lb[keep]; sub$ub <- gv$ub[keep]; sub$reaction_meta <- meta[match(keep, meta$reaction_id), , drop = FALSE]
  mets_used <- rownames(sub$S)[Matrix::rowSums(abs(sub$S) > 0) > 0]; sub$S <- sub$S[mets_used, , drop = FALSE]
  if (!is.null(gem$metabolite_meta)) sub$metabolite_meta <- gem$metabolite_meta[match(mets_used, as.character(gem$metabolite_meta$metabolite_id)), , drop = FALSE]
  if (!is.null(gem$gpr_table)) sub$gpr_table <- gem$gpr_table[as.character(gem$gpr_table$reaction_id) %in% keep, , drop = FALSE]
  med_diag <- data.frame(); if (!is.null(medium_table)) { app <- rc_apply_medium_constraints(sub, medium_table, condition = condition, strict = FALSE); sub <- app$gem; med_diag <- app$medium_diagnostics }
  sub$target_reaction <- target_reaction; sub$reaction_roles <- sub$reaction_meta[, intersect(c("reaction_id", "role", "role_source", "role_confidence"), colnames(sub$reaction_meta)), drop = FALSE]
  sub$medium_diagnostics <- med_diag
  sub$closure_diagnostics <- rc_check_microgem_closure(sub, target_reaction)
  sub$build_params <- list(k_hop = k_hop, max_reactions = max_reactions, strict_closure = strict_closure)
  sub
}
#' Check strict closure diagnostics for a micro-GEM
#' @export
rc_check_microgem_closure <- function(microgem, target_reaction, solver = "highs", flux_threshold = 1e-8) {
  gv <- rc_validate_gem(microgem)
  f <- rc_compass_vmax_directional(gv$S, gv$lb, gv$ub, target_reaction, "forward", solver, 60)
  role <- if (!is.null(microgem$reaction_meta) && "role" %in% colnames(microgem$reaction_meta)) as.character(microgem$reaction_meta$role) else rep("unknown", length(gv$reactions))
  nnz_m <- Matrix::rowSums(abs(gv$S) > 0)
  data.frame(target_reaction = target_reaction, strict_target_feasible = isTRUE(f$feasible), strict_vmax = f$vmax, n_boundary_metabolites = sum(nnz_m == 1), n_deadend_metabolites = sum(nnz_m <= 1), n_exchange_reactions = sum(role == "exchange"), n_transport_reactions = sum(role == "transport"), n_support_reactions = sum(role %in% c("exchange", "demand", "sink", "artificial_support", "cofactor_recycle")), top_unbalanced_boundary_metabolites = paste(utils::head(rownames(gv$S)[nnz_m == 1], 5), collapse = ";"), closure_warning_flag = !isTRUE(f$feasible), stringsAsFactors = FALSE)
}
