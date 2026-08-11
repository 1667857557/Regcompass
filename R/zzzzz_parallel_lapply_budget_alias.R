# Activate the bounded dispatcher defined in zzzz_condition_parallel_budget.R.
# Kept as a late-loaded adapter so the canonical parallel helpers remain the
# source of backend construction and worker-cap validation.

.rc_parallel_lapply_budget_core <- .rc_parallel_lapply

rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (!length(X)) return(list())
  .rc_parallel_lapply_budget_core(X, FUN, BPPARAM = BPPARAM, ...)
}
