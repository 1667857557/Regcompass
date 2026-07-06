#' Run selected flux variability analysis for a micro-GEM
#' @export
rc_run_selected_fva <- function(microgem, target_reaction, reactions = NULL, target_direction = "forward", fraction_of_vmax = 0.95, penalty_bound = NULL, solver = "highs", time_limit = 60) {
  gv <- rc_validate_gem(microgem); if (is.null(reactions)) reactions <- target_reaction; reactions <- intersect(reactions, gv$reactions)
  vmax <- rc_compass_vmax_directional(gv$S, gv$lb, gv$ub, target_reaction, target_direction, solver, time_limit)
  baseA <- gv$S; lhs <- rhs <- rep(0, nrow(gv$S))
  if (isTRUE(vmax$feasible)) { j <- match(target_reaction, gv$reactions); at <- rep(0, ncol(gv$S)); at[j] <- if (target_direction == "reverse") -1 else 1; baseA <- rbind(baseA, at); lhs <- c(lhs, fraction_of_vmax*vmax$vmax); rhs <- c(rhs, Inf) }
  out <- lapply(reactions, function(r) {
    j <- match(r, gv$reactions); obj <- rep(0, ncol(gv$S)); obj[j] <- 1
    mn <- rc_solve_lp(obj, baseA, lhs, rhs, gv$lb, gv$ub, solver, time_limit)
    mx <- rc_solve_lp(-obj, baseA, lhs, rhs, gv$lb, gv$ub, solver, time_limit)
    fmin <- if (identical(mn$status,"optimal")) mn$objective else NA_real_; fmax <- if (identical(mx$status,"optimal")) -mx$objective else NA_real_
    data.frame(reaction_id=r, unit_id=NA_character_, fva_min=fmin, fva_max=fmax, range_width=fmax-fmin, blocked_flag=!isTRUE(vmax$feasible), alternative_route_flag=is.finite(fmax-fmin) && (fmax-fmin) > 1e-6, stringsAsFactors=FALSE)
  })
  do.call(rbind, out)
}
