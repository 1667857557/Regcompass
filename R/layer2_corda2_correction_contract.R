# Record intentional corrections to the pinned Python CORDA2 implementation.

.rc_corda_build_three_stage_before_correction_contract <-
  .rc_corda_build_three_stage

.rc_corda_build_three_stage <- function(
    split, classes, options, solver, time_limit) {
  answer <- .rc_corda_build_three_stage_before_correction_contract(
    split = split,
    classes = classes,
    options = options,
    solver = solver,
    time_limit = time_limit
  )
  answer$corrected_python_defects <- unique(c(
    answer$corrected_python_defects %||% character(),
    "remaining medium confidence flux is maximized",
    "opposite reversible target direction is blocked",
    "penalties use each directional variable's own current confidence"
  ))
  answer$directional_penalty_contract <- paste(
    "forward and reverse variables are assigned costs from their own",
    "current directional confidence; confidence promotion remains directional"
  )
  answer
}
