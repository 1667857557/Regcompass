# Bootstrap-quality argument contract loaded after multitask_grn.R.

.rc_multitask_grn_defaults_core <- .rc_multitask_grn_defaults
.rc_validate_multitask_grn_args_core <- .rc_validate_multitask_grn_args

.rc_multitask_grn_defaults <- function() {
  out <- .rc_multitask_grn_defaults_core()
  out$min_bootstrap_success_fraction <- 0.8
  out
}

.rc_validate_multitask_grn_args <- function(args = list()) {
  out <- .rc_validate_multitask_grn_args_core(args)
  if (!is.numeric(out$alpha) || length(out$alpha) != 1L ||
      !is.finite(out$alpha) || out$alpha <= 0 || out$alpha >= 1) {
    stop(
      paste(
        "`multitask_args$alpha` must be in (0, 1): a positive lasso",
        "component is required for sparse bootstrap selection, and a positive",
        "ridge component is required for a unique symmetric deviation solution."
      ),
      call. = FALSE
    )
  }
  value <- out$min_bootstrap_success_fraction
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0 || value > 1) {
    stop(
      "`multitask_args$min_bootstrap_success_fraction` must be in (0, 1].",
      call. = FALSE
    )
  }
  out$min_bootstrap_success_fraction <- as.numeric(value)
  if (!identical(as.numeric(out$candidate_screen_threshold), 0)) {
    stop(
      paste(
        "`multitask_args$candidate_screen_threshold` must be 0 in the",
        "canonical model. Outcome-based full-data correlation screening would",
        "make cross-validated target reliability optimistic; use Pando's",
        "structural detection filters instead."
      ),
      call. = FALSE
    )
  }
  out$candidate_screen_threshold <- 0
  out
}
