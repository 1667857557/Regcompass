#' Construct a minimal GEM object from stoichiometry and bounds
#' @export
rc_make_gem <- function(S, lb = NULL, ub = NULL, reaction_meta = NULL, metabolite_meta = NULL, medium_policy = "base_bounds") {
  S <- as.matrix(S)
  if (is.null(colnames(S))) stop("`S` must have reaction IDs in colnames.", call. = FALSE)
  if (is.null(rownames(S))) rownames(S) <- paste0("met_", seq_len(nrow(S)))
  gem <- list(S = S, lb = lb, ub = ub, reaction_meta = reaction_meta, metabolite_meta = metabolite_meta, medium_policy = medium_policy)
  gv <- rc_validate_gem(gem)
  gem$lb <- gv$lb; gem$ub <- gv$ub; gem
}

#' Read a GEM object stored as an RDS file
#' @export
rc_read_gem <- function(file) {
  gem <- readRDS(file)
  rc_validate_gem(gem)
  gem
}
