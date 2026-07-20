# Internal parallel map used by canonical matrix calculations.
rc_internal_lapply <- function(X, FUN, BPPARAM = NULL) {
  rc_parallel_lapply(X, FUN, BPPARAM = BPPARAM)
}
