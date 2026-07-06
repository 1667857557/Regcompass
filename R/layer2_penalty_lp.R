#' Run COMPASS-like two-step penalty LP for RegCompassR Layer 2
#' @export
rc_run_layer2_compass_lp <- function(layer1, gem, unit = c("sample_celltype", "pool"),
                                     sample_col = "sample_id", celltype_col = "cell_type", condition_col = "condition",
                                     selected_reactions = NULL,
                                     selection_method = c("auto", "top", "differential", "pathway", "custom"),
                                     top_n = 300, min_C_rel = 0.15, min_confidence = 0.25,
                                     neighbor_depth = 1, max_subgem_reactions = 1000,
                                     omega = 0.95, penalty_epsilon = 1e-6, penalty_cap = 20,
                                     solver = c("gurobi", "highs", "glpk"), time_limit = 60, BPPARAM = NULL) {
  unit <- match.arg(unit); solver <- match.arg(solver); selection_method <- match.arg(selection_method)
  if (is.null(layer1$C_rel)) stop("`layer1` must contain `C_rel`.", call. = FALSE)
  if (is.null(layer1$reaction_confidence)) stop("`layer1` must contain `reaction_confidence`.", call. = FALSE)
  gv <- rc_validate_gem(gem)
  mats <- rc_layer2_unit_matrices(layer1, unit, sample_col, celltype_col, condition_col)
  sel <- rc_select_layer2_reactions(layer1, gem, selected_reactions, selection_method, top_n, min_C_rel, min_confidence, neighbor_depth, max_subgem_reactions)
  sub_rxns <- sel$reaction_id
  S <- gv$S[, sub_rxns, drop = FALSE]; lb <- gv$lb[sub_rxns]; ub <- gv$ub[sub_rxns]
  C <- rc_align_layer2_evidence(mats$C_rel, sub_rxns, fill = 1)
  Conf <- rc_align_layer2_evidence(mats$reaction_confidence, sub_rxns, fill = 1)
  pen <- rc_layer2_penalty(C, Conf, epsilon = penalty_epsilon, penalty_cap = penalty_cap)
  nr <- length(sub_rxns); nu <- ncol(C)
  penalty_mat <- vmax_mat <- matrix(NA_real_, nr, nu, dimnames = list(sub_rxns, colnames(C)))
  feasible_mat <- matrix(FALSE, nr, nu, dimnames = list(sub_rxns, colnames(C)))
  status_mat <- matrix("not_run", nr, nu, dimnames = list(sub_rxns, colnames(C)))
  diagnostics <- list()
  for (u in colnames(C)) for (r in sub_rxns) {
    p <- pen$penalty[, u]
    ans <- rc_compass_two_step_lp(S, lb, ub, target_reaction = r, penalties = p, omega = omega, solver = solver, time_limit = time_limit)
    vmax_mat[r, u] <- ans$vmax; penalty_mat[r, u] <- ans$penalty; feasible_mat[r, u] <- isTRUE(ans$feasible)
    status_mat[r, u] <- ans$solver_status
  }
  score_mat <- rc_compass_score_from_penalty(penalty_mat, feasible_mat, epsilon = penalty_epsilon)
  list(L2_compass_like_score = score_mat, L2_compass_like_penalty = penalty_mat, L2_vmax_internal = vmax_mat,
       L2_feasible_flag = feasible_mat, L2_solver_status = status_mat, penalty_components = pen$components,
       selected_reactions = sel, subgem_diagnostics = data.frame(n_reactions_subgem = nr, n_metabolites_subgem = nrow(S), blocked_reaction_flag = rowSums(feasible_mat) == 0, reaction_id = sub_rxns),
       medium_policy = if (is.null(gem$medium_policy)) "base_bounds" else gem$medium_policy, unit_meta = mats$unit_meta,
       layer1_summary_used = mats$summary, method = "COMPASS-like two-step penalty LP")
}

rc_layer2_penalty <- function(C_rel, Conf, epsilon = 1e-6, epsilon_C = 1e-3, epsilon_Conf = 1e-3, penalty_cap = 20) {
  Cc <- pmax(as.matrix(C_rel), epsilon_C, na.rm = TRUE); Fc <- pmax(as.matrix(Conf), epsilon_Conf, na.rm = TRUE)
  E <- Cc * Fc; P <- pmin(-log(E + epsilon), penalty_cap); P[!is.finite(P)] <- penalty_cap
  list(penalty = P, components = list(C_rel = C_rel, reaction_confidence = Conf, evidence_product = E, penalty = P, epsilon_C = epsilon_C, epsilon_Conf = epsilon_Conf, penalty_cap = penalty_cap))
}

rc_compass_score_from_penalty <- function(P, feasible, epsilon = 1e-6) {
  score <- P
  for (i in seq_len(nrow(P))) {
    x <- P[i, ]; med <- stats::median(x[feasible[i, ]], na.rm = TRUE); sc <- stats::mad(x[feasible[i, ]], constant = 1.4826, na.rm = TRUE)
    if (!is.finite(sc) || sc <= 0) sc <- epsilon
    z <- (med - x) / (sc + epsilon); score[i, ] <- rc_sigmoid(z)
  }
  score[!feasible] <- 0; score[!is.finite(score)] <- 0; score
}

rc_layer2_unit_matrices <- function(layer1, unit, sample_col, celltype_col, condition_col) {
  C <- as.matrix(layer1$C_rel); Conf <- as.matrix(layer1$reaction_confidence)
  common <- intersect(rownames(C), rownames(Conf)); C <- C[common,,drop=FALSE]; Conf <- Conf[common,,drop=FALSE]
  if (unit == "pool") return(list(C_rel = C, reaction_confidence = Conf, unit_meta = layer1$pool_meta, summary = "pool-level Layer1 matrices"))
  if (is.null(layer1$pool_meta)) stop("sample_celltype units require `layer1$pool_meta`.", call. = FALSE)
  pm <- layer1$pool_meta; required <- c("pool_id", sample_col, celltype_col); missing <- setdiff(required, colnames(pm)); if (length(missing)) stop("pool_meta missing: ", paste(missing, collapse=", "), call.=FALSE)
  pm <- pm[pm$pool_id %in% colnames(C), , drop=FALSE]
  group_cols <- c(sample_col, if (condition_col %in% colnames(pm)) condition_col else NULL, celltype_col)
  gid <- interaction(pm[, group_cols, drop=FALSE], sep="|", drop=TRUE)
  agg <- function(M) do.call(cbind, lapply(split(pm$pool_id, gid), function(cols) rowMedians_safe(M[, cols, drop=FALSE])))
  outC <- agg(C); outF <- agg(Conf); unit_meta <- unique(data.frame(unit_id = as.character(gid), pm[, group_cols, drop=FALSE], stringsAsFactors=FALSE))
  rownames(unit_meta) <- NULL; list(C_rel = outC, reaction_confidence = outF, unit_meta = unit_meta, summary = "Layer1 aggregated to sample × celltype medians")
}

rc_align_layer2_evidence <- function(M, rxns, fill = 1) {
  M <- as.matrix(M)
  out <- matrix(fill, nrow = length(rxns), ncol = ncol(M), dimnames = list(rxns, colnames(M)))
  common <- intersect(rxns, rownames(M))
  out[common, ] <- M[common, , drop = FALSE]
  out
}
