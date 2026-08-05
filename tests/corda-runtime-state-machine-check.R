# Run the existing pinned-Python oracle suite first with the baseline
# full-vector persistent engine. This leaves its synthetic models, helpers and
# baseline results in the current global environment.
source("tests/corda-synthetic-check.R", local = FALSE)

baseline <- list(
  stage3 = stage3,
  support = support,
  forced = forced,
  free = free
)

# The runtime file captures these production functions at source time. They are
# not exercised by this isolated state-machine test, but must exist so the same
# package runtime override is loaded rather than a test-specific copy.
rc_regcompass_step_layer2 <- function(...) NULL
.rc_build_celltype_medium_union_gem_cache <- function(...) NULL
source("R/layer2_corda_runtime.R", local = FALSE)

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
  "source_semantics",
  "source_fidelity",
  "intentional_corrections"
)

compare_state <- function(label, reference, candidate) {
  for (field in mathematical_fields) {
    if (!identical(reference[[field]], candidate[[field]])) {
      stop(
        "Sparse persistent CORDA2 changed mathematical field `",
        field, "` in case `", label, "`.",
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
      "Sparse persistent CORDA2 changed task-level mathematical diagnostics ",
      "in case `", label, "`.",
      call. = FALSE
    )
  }

  performance <- candidate$solver_performance
  stopifnot(
    is.list(performance),
    performance$n_solves > 0L,
    performance$n_sparse_update_calls > 0L,
    performance$n_transmitted_numeric_values <
      performance$n_full_vector_numeric_values,
    performance$n_full_vector_numeric_values_avoided > 0,
    performance$transmitted_fraction_of_full < 1
  )
  invisible(TRUE)
}

sparse <- list(
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

for (label in names(baseline)) {
  compare_state(label, baseline[[label]], sparse[[label]])
  compare_oracle_confidence(label_map <- switch(
    label,
    stage3 = "stage3_unknown",
    support = "absent_support",
    forced = "positive_min_forced",
    free = "positive_min_free"
  ), sparse[[label]])
}

cat(
  "Full CORDA2 state machine is identical under sparse persistent updates\n"
)
