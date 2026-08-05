test_that("CORDA2 task storage keeps directional association IDs", {
  task <- data.frame(
    variable_id = "R::forward",
    reaction_id = "R",
    direction = "forward",
    stage = "corda2_stage1_high_associations",
    replicate = 1L,
    kind = "dependency",
    status = "optimal",
    target_flux = 1,
    objective = 100,
    backend = "highs_persistent_cpp_basis_reuse",
    solver_message = "optimal",
    noise_namespace = "not_applicable_to_python_corda2",
    opposite_direction_blocked = "",
    n_associated = 1L,
    associated = "A::reverse",
    corda2_redundancies = 1L,
    corda2_n_solves = 1L,
    stringsAsFactors = FALSE
  )
  associations <- RegCompassR:::.rc_corda_normalize_associations(task)
  tasks <- task
  tasks$task_key <- RegCompassR:::.rc_corda_task_keys(tasks)
  tasks$associated <- NULL
  summary <- RegCompassR:::.rc_corda_task_summary(tasks)

  expect_equal(nrow(tasks), 1L)
  expect_equal(nrow(summary), 1L)
  expect_identical(associations$associated_variable_id, "A::reverse")
  expect_identical(associations$associated_reaction_id, "A")
  expect_identical(associations$task_key, tasks$task_key)
})

test_that("empty association output has a stable three-column schema", {
  empty <- data.frame(
    task_key = character(),
    associated_reaction_id = character(),
    stringsAsFactors = FALSE
  )
  normalized <- RegCompassR:::.rc_corda_normalize_associations(empty)
  expect_identical(
    colnames(normalized),
    c("task_key", "associated_variable_id", "associated_reaction_id")
  )
  expect_equal(nrow(normalized), 0L)
})
