#' Apply a simple medium policy to GEM bounds
#' @export
rc_layer2_apply_bounds <- function(gem, medium_policy = NULL) {
  gv <- rc_validate_gem(gem)
  out <- gem; out$lb <- gv$lb; out$ub <- gv$ub; out$medium_policy <- if (is.null(medium_policy)) "base_bounds" else medium_policy; out
}
