# Dispatch between the preserved serial CORDA2 core and stage-parallel core.

.rc_corda_build_three_stage_dispatch <- function(...) {
  if (.rc_corda_stage_parallel_requested()) {
    return(.rc_corda_build_three_stage(...))
  }
  progress_state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  answer <- do.call(
    .rc_corda_build_three_stage_serial_core,
    list(...)
  )
  if (!is.environment(progress_state)) return(answer)
  task <- progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    required <- list(
      c("corda2_step1", "corda2_step1_HC_dependencies", 4L),
      c("corda2_step2_1", "corda2_step2_1_MC_NC_dependencies", 5L),
      c("corda2_step2_2", "corda2_step2_2_MC_feasibility", 6L),
      c("corda2_step3", "corda2_step3_HC_OT_dependencies", 7L)
    )
    for (item in required) {
      if (!exists(item[[1L]], envir = progress_state$algorithm_flags,
                  inherits = FALSE)) {
        .rc_layer2_algorithm_once(
          item[[1L]], item[[2L]], as.integer(item[[3L]]),
          "skipped because this confidence class had no candidate directions"
        )
      }
    }
  }
  answer
}
