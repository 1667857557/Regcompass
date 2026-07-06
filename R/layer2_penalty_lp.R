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
  C <- rc_align_layer2_evidence(mats$C_rel, sub_rxns, fill = NA_real_)
  Conf <- rc_align_layer2_evidence(mats$reaction_confidence, sub_rxns, fill = NA_real_)
  support_rxns <- rc_layer2_support_reactions(gem, sub_rxns)
  pen <- rc_layer2_penalty(C, Conf, epsilon = penalty_epsilon, penalty_cap = penalty_cap, support_reactions = support_rxns)
  nr <- length(sub_rxns); nu <- ncol(C)
  penalty_mat <- vmax_mat <- matrix(NA_real_, nr, nu, dimnames = list(sub_rxns, colnames(C)))
  feasible_mat <- matrix(FALSE, nr, nu, dimnames = list(sub_rxns, colnames(C)))
  status_mat <- matrix("not_run", nr, nu, dimnames = list(sub_rxns, colnames(C)))
  diagnostics <- vector("list", nr * nu)
  diag_i <- 0L
  for (u in colnames(C)) for (r in sub_rxns) {
    p <- pen$penalty[, u]
    ans <- rc_compass_two_step_lp(S, lb, ub, target_reaction = r, penalties = p, omega = omega, solver = solver, time_limit = time_limit)
    vmax_mat[r, u] <- ans$vmax; penalty_mat[r, u] <- ans$penalty; feasible_mat[r, u] <- isTRUE(ans$feasible)
    status_mat[r, u] <- ans$solver_status
    diag_i <- diag_i + 1L
    diagnostics[[diag_i]] <- data.frame(reaction_id = r, unit_id = u, solver_status = ans$solver_status,
                                        step1_status = ans$step1_status, step2_status = ans$step2_status,
                                        objective_value = ans$penalty, vmax = ans$vmax,
                                        number_constraints = ans$number_constraints, number_variables = ans$number_variables,
                                        runtime = ans$runtime, stringsAsFactors = FALSE)
  }
  score_mat <- rc_compass_score_from_penalty(penalty_mat, feasible_mat, epsilon = penalty_epsilon)
  list(L2_compass_like_score = score_mat, L2_compass_like_penalty = penalty_mat, L2_vmax_internal = vmax_mat,
       L2_feasible_flag = feasible_mat, L2_solver_status = status_mat, penalty_components = pen$components,
       selected_reactions = sel, subgem_diagnostics = data.frame(n_reactions_subgem = nr, n_metabolites_subgem = nrow(S), blocked_reaction_flag = rowSums(feasible_mat) == 0, reaction_id = sub_rxns),
       lp_diagnostics = do.call(rbind, diagnostics),
       medium_policy = if (is.null(gem$medium_policy)) "base_bounds" else gem$medium_policy, unit_meta = mats$unit_meta,
       layer1_summary_used = mats$summary, method = "COMPASS-like two-step penalty LP")
}

rc_layer2_penalty <- function(C_rel, Conf, epsilon = 1e-6, epsilon_C = 1e-3, epsilon_Conf = 1e-3,
                              penalty_cap = 20, support_reactions = character(), support_penalty = 0) {
  C_raw <- as.matrix(C_rel); Conf_raw <- as.matrix(Conf)
  Cc <- ifelse(is.na(C_raw), epsilon_C, pmax(C_raw, epsilon_C))
  Fc <- ifelse(is.na(Conf_raw), epsilon_Conf, pmax(Conf_raw, epsilon_Conf))
  E <- Cc * Fc
  P <- pmax(0, pmin(-log(E + epsilon), penalty_cap))
  P[!is.finite(P)] <- penalty_cap
  support_reactions <- intersect(support_reactions, rownames(P))
  if (length(support_reactions)) P[support_reactions, ] <- support_penalty
  missing_evidence <- is.na(C_raw) | is.na(Conf_raw)
  list(penalty = P, components = list(C_rel = C_rel, reaction_confidence = Conf, evidence_product = E,
                                      missing_evidence = missing_evidence, support_reaction = rownames(P) %in% support_reactions,
                                      penalty = P, epsilon_C = epsilon_C, epsilon_Conf = epsilon_Conf,
                                      penalty_cap = penalty_cap, support_penalty = support_penalty))
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
  C <- as.matrix(layer1$C_rel); Conf <- rc_layer2_confidence_matrix(layer1$reaction_confidence, C)
  common <- intersect(rownames(C), rownames(Conf)); C <- C[common,,drop=FALSE]; Conf <- Conf[common,,drop=FALSE]
  if (unit == "pool") return(list(C_rel = C, reaction_confidence = Conf, unit_meta = layer1$pool_meta, summary = "pool-level Layer1 matrices"))
  if (is.null(layer1$pool_meta)) stop("sample_celltype units require `layer1$pool_meta`.", call. = FALSE)
  pm <- layer1$pool_meta; required <- c("pool_id", sample_col, celltype_col); missing <- setdiff(required, colnames(pm)); if (length(missing)) stop("pool_meta missing: ", paste(missing, collapse=", "), call.=FALSE)
  pm <- pm[pm$pool_id %in% colnames(C), , drop=FALSE]
  group_cols <- c(sample_col, if (condition_col %in% colnames(pm)) condition_col else NULL, celltype_col)
  gid <- interaction(pm[, group_cols, drop=FALSE], sep="|", drop=TRUE)
  agg <- function(M) do.call(cbind, lapply(split(pm$pool_id, gid), function(cols) rowMedians_safe(M[, cols, drop=FALSE])))
  outC <- agg(C); outF <- agg(Conf)
  unit_meta <- unique(data.frame(unit_id = as.character(gid), pm[, group_cols, drop=FALSE], stringsAsFactors=FALSE))
  unit_meta <- unit_meta[match(colnames(outC), unit_meta$unit_id), , drop = FALSE]
  rownames(unit_meta) <- NULL; list(C_rel = outC, reaction_confidence = outF, unit_meta = unit_meta, summary = "Layer1 aggregated to sample × celltype medians")
}

rc_align_layer2_evidence <- function(M, rxns, fill = 1) {
  M <- as.matrix(M)
  out <- matrix(fill, nrow = length(rxns), ncol = ncol(M), dimnames = list(rxns, colnames(M)))
  common <- intersect(rxns, rownames(M))
  out[common, ] <- M[common, , drop = FALSE]
  out
}

rc_layer2_support_reactions <- function(gem, rxns) {
  meta <- gem$reaction_meta
  if (is.null(meta) || !is.data.frame(meta) || !"reaction_id" %in% colnames(meta)) return(character())
  text_cols <- intersect(c("type", "reaction_type", "category", "subsystem", "name", "reaction_name"), colnames(meta))
  if (!length(text_cols)) return(character())
  txt <- apply(meta[, text_cols, drop = FALSE], 1, paste, collapse = " ")
  support <- grepl("exchange|demand|sink|transport|support", txt, ignore.case = TRUE)
  intersect(as.character(meta$reaction_id[support]), rxns)
}

rc_layer2_confidence_matrix <- function(confidence, C_rel) {
  C_rel <- as.matrix(C_rel)
  if (is.data.frame(confidence) && all(c("reaction_id", "pool_id", "reaction_confidence") %in% colnames(confidence))) {
    out <- matrix(NA_real_, nrow = nrow(C_rel), ncol = ncol(C_rel), dimnames = dimnames(C_rel))
    rid <- as.character(confidence$reaction_id)
    pid <- as.character(confidence$pool_id)
    ok <- rid %in% rownames(out) & pid %in% colnames(out)
    if (any(ok)) {
      idx <- cbind(match(rid[ok], rownames(out)), match(pid[ok], colnames(out)))
      out[idx] <- as.numeric(confidence$reaction_confidence[ok])
    }
    return(out)
  }
  M <- as.matrix(confidence)
  if (is.null(rownames(M)) || is.null(colnames(M))) stop("`reaction_confidence` must be a matrix-like object with dimnames or a long data frame.", call. = FALSE)
  M
}
