test_that("full GEM cache reuses complete model per medium scenario", {
  S <- matrix(c(
    -1, 0, 0,
     1,-1, 0,
     0, 1,-1
  ), nrow = 3, dimnames = list(c("A_c", "B_c", "C_c"), c("R1", "R2", "R3")))
  reaction_meta <- data.frame(reaction_id = c("R1", "R2", "R3"), metabolic_module = c("a", "a", "b"), stringsAsFactors = FALSE)
  gem <- rc_make_gem(S, lb = c(R1 = 0, R2 = 0, R3 = 0), ub = c(R1 = 1000, R2 = 1000, R3 = 1000), reaction_meta = reaction_meta)
  dirs <- data.frame(reaction_id = c("R1", "R2"), target_direction = c("forward", "forward"), stringsAsFactors = FALSE)
  medium <- data.frame(medium_scenario_id = "base", exchange_reaction_id = character(), lb = numeric(), ub = numeric(), available = logical(), stringsAsFactors = FALSE)

  cache <- rc_build_full_gem_cache(gem = gem, dirs = dirs, medium_scenarios = medium)
  summary <- attr(cache, "summary")
  expect_equal(nrow(summary), 1)
  expect_equal(summary$build_strategy, "full_gem")
  expect_true(file.exists(summary$file[[1]]))
  expect_setequal(colnames(readRDS(summary$file[[1]])$S), colnames(S))
})
