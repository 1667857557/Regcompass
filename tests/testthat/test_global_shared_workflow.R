test_that("global union collapses stratum modules into one structural identity", {
  membership <- data.frame(
    sample_id = c("S1", "S2", "S2"), module_id = c("M1", "M2", "M3"),
    reaction_id = c("R1", "R1", "R2"), stringsAsFactors = FALSE
  )
  core <- data.frame(
    sample_id = c("S1", "S2"), module_id = c("M1", "M3"),
    reaction_id = c("R1", "R2"), is_core = c(TRUE, TRUE), stringsAsFactors = FALSE
  )
  out <- .rc_global_union_tables(membership, core)
  expect_setequal(out$membership$reaction_id, c("R1", "R2"))
  expect_true(all(out$membership$sample_id == "GLOBAL"))
  expect_true(all(out$membership$module_id == "GLOBAL_UNION"))
  expect_setequal(out$core$reaction_id, c("R1", "R2"))
})

test_that("shared GEM rejects condition-specific medium bounds", {
  medium <- data.frame(
    medium_scenario_id = "m", exchange_reaction_id = "EX",
    lb = c(-1, -10), ub = 0, available = TRUE,
    condition = c("A", "B"), stringsAsFactors = FALSE
  )
  expect_error(.rc_validate_shared_medium(medium), "condition-independent")
  medium$condition <- "all"
  expect_invisible(.rc_validate_shared_medium(medium))
})

test_that("integrated workflow has an upstream barrier before global calibration and LP", {
  text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  upstream <- regexpr("rc_parallel_lapply", text, fixed = TRUE)[[1L]]
  barrier <- regexpr("Global barrier not reached", text, fixed = TRUE)[[1L]]
  calibration <- regexpr(".rc_merge_completed_layer1", text, fixed = TRUE)[[1L]]
  layer2 <- regexpr("rc_run_microcompass", text, fixed = TRUE)[[1L]]
  expect_true(all(c(upstream, barrier, calibration, layer2) > 0L))
  expect_lt(upstream, barrier)
  expect_lt(barrier, calibration)
  expect_lt(calibration, layer2)
})

test_that("Layer 2 defaults to metacells and evaluates the full shared task grid", {
  expect_identical(eval(formals(rc_run_microcompass)$unit)[[1L]], "metacell")
  text <- paste(deparse(body(rc_run_microcompass)), collapse = "\n")
  expect_match(text, "expand.grid\\(row_id = row_ids, unit_id = units", perl = TRUE)
  expect_false(grepl("unit_sample ==", text, fixed = TRUE))
})

test_that("parallel stages stop their worker backend", {
  text <- paste(deparse(body(rc_parallel_lapply)), collapse = "\n")
  expect_match(text, ".rc_stop_bpparam", fixed = TRUE)
})
