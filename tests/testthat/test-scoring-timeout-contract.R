test_that("scoring APIs expose no time-limit or sample parameter", {
  expect_false("time_limit" %in% names(formals(rc_compass_two_step_lp_directional)))
  expect_false("time_limit" %in% names(formals(rc_run_microcompass)))
  expect_false("time_limit" %in% names(formals(.rc_run_microcompass_engine)))
  expect_false("time_limit" %in% names(formals(.rc_score_existing_union_cache)))
  expect_false("sample_col" %in% names(formals(rc_run_microcompass)))
  expect_false("sample_col" %in% names(formals(.rc_score_existing_union_cache)))
})

test_that("only union-GEM construction receives the completion limit", {
  engine <- paste(deparse(body(.rc_run_microcompass_engine)), collapse = "\n")
  scoring <- paste(
    deparse(body(rc_compass_two_step_lp_directional)), collapse = "\n"
  )
  target_scoring <- paste(
    deparse(body(.rc_score_existing_union_cache_sample_core)), collapse = "\n"
  )

  expect_match(engine, "model_params$completion_time_limit", fixed = TRUE)
  expect_match(engine, ".rc_build_medium_specific_union_gem_cache", fixed = TRUE)
  expect_false(grepl("time_limit =", scoring, fixed = TRUE))
  expect_false(grepl("time_limit =", target_scoring, fixed = TRUE))
})

test_that("Stage 5 rejects retired timeout and sample arguments", {
  body_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  expect_match(body_text, "Scoring `time_limit`, sample aggregation", fixed = TRUE)
  expect_match(body_text, "completion_time_limit", fixed = TRUE)
  expect_match(body_text, "allowed_model_params", fixed = TRUE)
  expect_match(body_text, 'unit = "metacell"', fixed = TRUE)
  expect_false(grepl("params$sample_col", body_text, fixed = TRUE))
})

test_that("target-union results record unlimited metacell scoring", {
  target_core <- paste(
    deparse(body(.rc_score_existing_union_cache_sample_core)), collapse = "\n"
  )
  active <- paste(
    deparse(body(.rc_score_existing_union_cache)), collapse = "\n"
  )
  expect_match(target_core, 'scoring_time_limit = "none"', fixed = TRUE)
  expect_match(active, 'answer$params$unit <- "metacell"', fixed = TRUE)
  expect_match(active, 'answer$params$aggregation <- "none"', fixed = TRUE)
})

test_that("user examples contain no standalone scoring time_limit", {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[vapply(
    roots,
    function(path) file.exists(file.path(path, "DESCRIPTION")),
    logical(1)
  )]
  if (!length(roots)) skip("Source documentation is unavailable.")
  root <- normalizePath(roots[[1L]], mustWork = TRUE)
  paths <- c(
    file.path(root, "README.md"),
    list.files(file.path(root, "docs"), pattern = "\\.md$", full.names = TRUE),
    file.path(root, "vignettes", "regcompass-workflow.Rmd")
  )
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl("(?m)^\\s*time_limit\\s*=", text, perl = TRUE))
  expect_match(text, "completion_time_limit =", fixed = TRUE)
})
