test_that("module meso-GEM cache is built once per module and medium", {
  S <- matrix(c(
    -1, 0,
     1,-1,
     0, 1
  ), nrow = 3, dimnames = list(c("A_c", "B_c", "C_c"), c("R1", "R2")))
  reaction_meta <- data.frame(reaction_id = c("R1", "R2"), metabolic_module = c("module_a", "module_a"), stringsAsFactors = FALSE)
  met_meta <- data.frame(metabolite_id = rownames(S), compartment = c("c", "c", "c"), stringsAsFactors = FALSE)
  gem <- rc_make_gem(S, lb = c(R1 = 0, R2 = 0), ub = c(R1 = 1000, R2 = 1000), reaction_meta = reaction_meta, metabolite_meta = met_meta)
  dirs <- data.frame(reaction_id = c("R1", "R2"), target_direction = c("forward", "forward"), stringsAsFactors = FALSE)
  medium <- data.frame(medium_scenario_id = "base", exchange_reaction_id = character(), lb = numeric(), ub = numeric(), available = logical(), stringsAsFactors = FALSE)
  cache <- rc_build_module_gem_cache(gem = gem, dirs = dirs, medium_scenarios = medium, module_col = "metabolic_module")
  summary <- attr(cache, "summary")
  expect_equal(nrow(summary), 1)
  expect_true(file.exists(summary$file[[1]]))
})

test_that("module meso-GEM max_reactions guard never truncates reactions", {
  S <- diag(3)
  dimnames(S) <- list(paste0("M", 1:3), paste0("R", 1:3))
  reaction_meta <- data.frame(reaction_id = paste0("R", 1:3), metabolic_module = "large_module", stringsAsFactors = FALSE)
  gem <- rc_make_gem(S, lb = rep(0, 3), ub = rep(1000, 3), reaction_meta = reaction_meta)
  mg <- rc_build_module_meso_gem(gem, "large_module", max_reactions = 1)
  expect_setequal(colnames(mg$S), paste0("R", 1:3))
  expect_true(mg$build_params$max_reactions_exceeded)
})
