# Exact-source audit contract.
# No mathematical correction is applied to the pinned Python CORDA2 behavior.

.rc_corda_build_three_stage_exact_base <- .rc_corda_build_three_stage

.rc_corda_build_three_stage <- function(
    split, classes, options, solver, time_limit) {
  answer <- .rc_corda_build_three_stage_exact_base(
    split = split,
    classes = classes,
    options = options,
    solver = solver,
    time_limit = time_limit
  )
  answer$source_fidelity <- "exact_for_met_prod_NULL"
  answer$intentional_corrections <- character()
  answer
}
