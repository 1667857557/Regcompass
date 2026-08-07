rc_make_gem <- function(S, lb = NULL, ub = NULL,
                        reaction_meta = NULL, metabolite_meta = NULL,
                        model_info = NULL, medium_policy = NULL,
                        gpr_table = NULL) {
  S <- RegCompassR:::.rc_as_dgCMatrix(S)
  if (is.null(colnames(S))) {
    colnames(S) <- paste0("R", seq_len(ncol(S)))
  }
  if (is.null(rownames(S))) {
    rownames(S) <- paste0("M", seq_len(nrow(S)))
  }
  rxns <- colnames(S)
  lb <- if (is.null(lb)) stats::setNames(rep(-1000, length(rxns)), rxns) else lb
  ub <- if (is.null(ub)) stats::setNames(rep(1000, length(rxns)), rxns) else ub
  gem <- list(
    S = S,
    lb = lb,
    ub = ub,
    reaction_meta = reaction_meta,
    metabolite_meta = metabolite_meta,
    model_info = model_info,
    medium_policy = medium_policy,
    gpr_table = gpr_table
  )
  RegCompassR:::rc_validate_gem(gem)
  gem
}
