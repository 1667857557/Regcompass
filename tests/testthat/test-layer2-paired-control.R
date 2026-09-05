test_that("Layer 2 primary and RNA control share one target dispatch", {
  stage <- paste(readLines(file.path("R", "step_layer2.R"), warn = FALSE),
                 collapse = "\n")
  full <- paste(readLines(file.path("R", "microcompass_engine.R"), warn = FALSE),
                collapse = "\n")
  cell <- paste(readLines(
    file.path("R", "celltype_microcompass_reaction_parallel.R"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(stage, "control_layer1 = control_layer1", fixed = TRUE)
  expect_match(stage, "paired_step2_dispatch = TRUE", fixed = TRUE)
  expect_false(grepl("rna_only <- run_control(", stage, fixed = TRUE))
  expect_false(grepl("rna_only <- answer", stage, fixed = TRUE))
  expect_match(stage, "answer$comparison_paths <- NULL", fixed = TRUE)
  expect_match(full, "control_penalties = control_penalties", fixed = TRUE)
  expect_match(cell, "control_penalties = control_penalties", fixed = TRUE)
  expect_match(full, ".rc_compass_step2_route_solve", fixed = TRUE)
  expect_match(cell, ".rc_compass_step2_route_solve", fixed = TRUE)
})
