#' Run relaxed mass-balance LP with metabolite slack diagnostics
#' @export
rc_run_relaxed_balance_lp <- function(microgem, penalties, target_reaction, target_direction = "forward", omega = 0.95, slack_penalty = NULL, slack_penalty_mode = c("relative", "absolute"), slack_penalty_multiplier = 5, metabolite_slack_weights = NULL, solver = "highs", time_limit = 60) {
  slack_penalty_mode <- match.arg(slack_penalty_mode)
  gv <- rc_validate_gem(microgem); n <- ncol(gv$S); m <- nrow(gv$S)
  vmax <- rc_compass_vmax_directional(gv$S, gv$lb, gv$ub, target_reaction, target_direction, solver, time_limit)
  p <- as.numeric(penalties[colnames(gv$S)]); p[!is.finite(p)] <- 20
  if (is.null(slack_penalty)) {
    slack_penalty <- if (slack_penalty_mode == "relative") slack_penalty_multiplier * max(p, na.rm = TRUE) else 100
  }
  w <- rep(1, m); names(w) <- rownames(gv$S); if (!is.null(metabolite_slack_weights)) w[names(metabolite_slack_weights)] <- metabolite_slack_weights
  Aeq <- cbind(gv$S, -gv$S, Matrix::Diagonal(m), -Matrix::Diagonal(m))
  lhs <- rhs <- rep(0, m)
  if (isTRUE(vmax$feasible)) {
    j <- match(target_reaction, colnames(gv$S)); at <- rep(0, 2*n + 2*m)
    if (target_direction == "reverse") { at[j] <- -1; at[n+j] <- 1 } else { at[j] <- 1; at[n+j] <- -1 }
    Aeq <- rbind(Aeq, at); lhs <- c(lhs, omega * vmax$vmax); rhs <- c(rhs, Inf)
  }
  ans <- rc_solve_lp(c(p, p, slack_penalty*w, slack_penalty*w), Aeq, lhs, rhs, c(rep(0,2*n), rep(0,2*m)), c(pmax(gv$ub,0), pmax(-gv$lb,0), rep(Inf,2*m)), solver, time_limit)
  sol <- ans$solution; sp <- sm <- rep(NA_real_, m)
  if (identical(ans$status, "optimal")) { sp <- sol[(2*n+1):(2*n+m)]; sm <- sol[(2*n+m+1):(2*n+2*m)] }
  total <- sum(sp + sm, na.rm = TRUE)
  role <- if (!is.null(microgem$metabolite_meta) && "is_currency" %in% colnames(microgem$metabolite_meta)) microgem$metabolite_meta$is_currency[match(rownames(gv$S), microgem$metabolite_meta$metabolite_id)] else FALSE
  data.frame(target_reaction=target_reaction, target_direction=target_direction, strict_feasible=isTRUE(vmax$feasible), relaxed_feasible=identical(ans$status,"optimal"), target_min_used=if (isTRUE(vmax$feasible)) omega * vmax$vmax else NA_real_, slack_penalty_used=slack_penalty, slack_penalty_multiplier=slack_penalty_multiplier, total_slack=total, top_slack_metabolites=paste(utils::head(rownames(gv$S)[order(sp+sm, decreasing=TRUE)],5), collapse=";"), cofactor_slack_used=sum((sp+sm)[role %in% TRUE], na.rm=TRUE), boundary_slack_used=NA_real_, medium_slack_used=NA_real_, solver_status=ans$status, stringsAsFactors=FALSE)
}
