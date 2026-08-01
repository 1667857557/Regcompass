# The canonical condition-GRN implementation is defined in
# condition_grn_contract.R. Preserve it and wrap the entry point only after the
# complete Pando runtime guard has been collated.
.rc_fit_condition_grns_by_cell_type_unchecked <-
  .rc_fit_condition_grns_by_cell_type

.rc_fit_condition_grns_by_cell_type <- function(..., BPPARAM = NULL) {
  .rc_require_pando_hybrid_runtime(BPPARAM = BPPARAM)
  .rc_fit_condition_grns_by_cell_type_unchecked(
    ...,
    BPPARAM = BPPARAM
  )
}
