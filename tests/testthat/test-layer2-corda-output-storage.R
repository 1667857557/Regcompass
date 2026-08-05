test_that("CORDA2 associations preserve variables and reaction IDs", {
  tab <- data.frame(
    variable_id = c("R1::forward", "R2::reverse"),
    reaction_id = c("R1", "R2"),
    direction = c("forward", "reverse"),
    stage = c("corda2_stage1_high_associations",
              "corda2_stage3_free_completion"),
    replicate = c(1L, 1L),
    kind = c("dependency", "dependency"),
    status = c("optimal", "optimal"),
    target_flux = c(1, 1),
    objective = c(2, 3),
    backend = c("highs", "highs"),
    solver_message = c("Optimal", "Optimal"),
    noise_namespace = c("not_applicable_to_corda2",
                        "not_applicable_to_corda2"),
    opposite_direction_blocked = c("", "R2::forward"),
    n_associated = c(2L, 1L),
    associated = c("S1::forward;S2::reverse", "S3::forward"),
    stringsAsFactors = FALSE
  )
  keys <- RegCompassR:::.rc_corda_task_keys(tab)
  edges <- RegCompassR:::.rc_corda_normalize_associations(tab)
  expect_length(keys, 2L)
  expect_false(anyDuplicated(keys))
  expect_equal(nrow(edges), 3L)
  expect_setequal(edges$task_key, keys)
  expect_setequal(
    edges$associated_variable_id,
    c("S1::forward", "S2::reverse", "S3::forward")
  )
  expect_setequal(edges$associated_reaction_id, c("S1", "S2", "S3"))
  expect_equal(nrow(edges), sum(tab$n_associated))
})

test_that("empty CORDA2 task storage remains schema-valid", {
  empty <- RegCompassR:::.rc_corda_empty_task_table()
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0L)
  expect_true(all(c(
    "task_key", "n_associated", "opposite_direction_blocked",
    "noise_namespace"
  ) %in% colnames(empty)))
  expect_identical(RegCompassR:::.rc_corda_task_keys(empty), character())
  edges <- RegCompassR:::.rc_corda_normalize_associations(empty)
  expect_equal(nrow(edges), 0L)
  expect_identical(
    colnames(edges),
    c("task_key", "associated_variable_id", "associated_reaction_id")
  )
})
