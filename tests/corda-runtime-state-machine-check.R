# The pinned-Python oracle suite defines deterministic synthetic models and the
# direct production CORDA2 functions.
source("tests/corda-synthetic-check.R", local = FALSE)

run_one_shot <- function(...) {
  original <- .rc_corda_highs_api_available
  assign(
    ".rc_corda_highs_api_available",
    function() FALSE,
    envir = .GlobalEnv
  )
  on.exit(assign(
    ".rc_corda_highs_api_available",
    original,
    envir = .GlobalEnv
  ), add = TRUE)
  run_build(...)$result
}

persistent <- list(
  stage3 = run_build(stage3_gem, c(SRC = 0L, H = 3L))$result,
  support = run_build(
    support_gem,
    c(N = -1L, M1 = 2L, M2 = 2L),
    n = 1L,
    support = 2L
  )$result,
  forced = run_build(
    forced_gem,
    c(SRC = 0L, M = 2L),
    n = 1L
  )$result,
  free = run_build(
    free_gem,
    c(SRC = 0L, M = 2L),
    n = 1L
  )$result
)

one_shot <- list(
  stage3 = run_one_shot(stage3_gem, c(SRC = 0L, H = 3L)),
  support = run_one_shot(
    support_gem,
    c(N = -1L, M1 = 2L, M2 = 2L),
    n = 1L,
    support = 2L
  ),
  forced = run_one_shot(
    forced_gem,
    c(SRC = 0L, M = 2L),
    n = 1L
  ),
  free = run_one_shot(
    free_gem,
    c(SRC = 0L, M = 2L),
    n = 1L
  )
)

mathematical_fields <- c(
  "included",
  "included_directional_variables",
  "initial_reaction_confidence",
  "final_reaction_confidence",
  "final_reaction_status",
  "final_confidence",
  "final_directional_confidence",
  "initial_directional_confidence",
  "confidence_after_stage1",
  "confidence_after_stage2_association",
  "confidence_after_stage2_support",
  "confidence_after_stage2_medium",
  "inclusion_stage",
  "inclusion_stage_direction",
  "stage1_associated",
  "stage2_nc_support_pairs",
  "stage2_nc_support_count",
  "stage2_promoted_nc",
  "stage2_promoted_mc",
  "stage3_associated_ot",
  "blocked_after_stage2",
  "blocked_before_stage3",
  "impossible_directional_targets",
  "redundancies",
  "algorithm",
  "python_reference_commit",
  "stage_update_policy",
  "source_semantics"
)

compare_state <- function(label, reference, candidate) {
  for (field in mathematical_fields) {
    if (!identical(reference[[field]], candidate[[field]])) {
      stop(
        "One-shot CORDA2 changed mathematical field `", field,
        "` in case `", label, "`.",
        call. = FALSE
      )
    }
  }

  reference_tasks <- reference$task_diagnostics
  candidate_tasks <- candidate$task_diagnostics
  drop_runtime <- intersect(
    c("backend", "solver_message"),
    union(colnames(reference_tasks), colnames(candidate_tasks))
  )
  reference_tasks <- reference_tasks[
    , setdiff(colnames(reference_tasks), drop_runtime), drop = FALSE
  ]
  candidate_tasks <- candidate_tasks[
    , setdiff(colnames(candidate_tasks), drop_runtime), drop = FALSE
  ]
  if (!isTRUE(all.equal(
    reference_tasks, candidate_tasks,
    tolerance = 1e-9,
    check.attributes = TRUE
  ))) {
    stop(
      "One-shot CORDA2 changed task-level mathematical diagnostics in case `",
      label, "`.",
      call. = FALSE
    )
  }

  persistent_performance <- reference$solver_performance
  one_shot_performance <- candidate$solver_performance
  stopifnot(
    isTRUE(persistent_performance$persistent_solver),
    persistent_performance$n_sparse_update_calls > 0L,
    persistent_performance$n_transmitted_numeric_values <
      persistent_performance$n_full_vector_numeric_values,
    identical(one_shot_performance$persistent_solver, FALSE),
    one_shot_performance$n_sparse_update_calls == 0L
  )
  invisible(TRUE)
}

for (label in names(persistent)) {
  compare_state(label, persistent[[label]], one_shot[[label]])
}

cat("Direct persistent and one-shot CORDA2 state machines are identical\n")
