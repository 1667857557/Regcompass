# The canonical condition-GRN implementation invokes the complete Pando runtime
# guard directly so that its explicit formals remain available to callers and
# the native self-test is executed only once. Keep the historical internal name
# as a compatibility alias without replacing the canonical implementation.
.rc_fit_condition_grns_by_cell_type_unchecked <-
  .rc_fit_condition_grns_by_cell_type
