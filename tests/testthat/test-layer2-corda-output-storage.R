test_that("CORDA2 task storage keeps original directional IDs", {
  task <- data.frame(
    variable_id = "R",
    reaction_id = "R",
    direction = "forward",
    stage = "corda2_step1_HC_dependencies",
    replicate = 1L,
    kind = "dependency",
    status = "optimal",
    target_flux = 1,
    vmax = 10,
    objective = 100,
    backend = "highs_persistent_cpp",
    solver_message = "optimal",
    opposite_direction_blocked = "R_CORDA_rev_rxn",
    n_associated = 1L,
    associated = "A_CORDA_rev_rxn",
    corda2_n_solves = 2L,
    stringsAsFactors = FALSE
  )
  associations <- RegCompassR:::.rc_corda_normalize_associations(task)
  tasks <- task
  tasks$task_key <- RegCompassR:::.rc_corda_task_keys(tasks)
  tasks$associated <- NULL
  summary <- RegCompassR:::.rc_corda_task_summary(tasks)

  expect_equal(nrow(tasks), 1L)
  expect_equal(nrow(summary), 1L)
  expect_identical(
    associations$associated_variable_id, "A_CORDA_rev_rxn"
  )
  expect_identical(associations$associated_reaction_id, "A")
  expect_identical(associations$task_key, tasks$task_key)
})

test_that("empty association output has a stable schema", {
  normalized <- RegCompassR:::.rc_corda_normalize_associations(data.frame())
  expect_identical(
    colnames(normalized),
    c("task_key", "associated_variable_id", "associated_reaction_id")
  )
  expect_equal(nrow(normalized), 0L)
})
