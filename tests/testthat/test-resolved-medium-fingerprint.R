make_resolved_medium_test_gem <- function() {
  reactions <- c("EX_glucose", "R1")
  list(
    S = Matrix::Matrix(
      matrix(c(-1, 1), nrow = 1,
             dimnames = list("m_e", reactions)),
      sparse = TRUE
    ),
    lb = stats::setNames(c(-1000, 0), reactions),
    ub = stats::setNames(c(1000, 1000), reactions),
    model_info = list(
      species = "human", source = "unit_test", model_version = "1"
    ),
    reaction_meta = data.frame(
      reaction_id = reactions,
      role = c("exchange", "internal"),
      metabolite_name = c("D-glucose", NA_character_),
      stringsAsFactors = FALSE
    )
  )
}

test_that("resolved medium fingerprint hashes final GEM bounds, not scenario labels", {
  gem <- make_resolved_medium_test_gem()
  medium <- data.frame(
    medium_scenario_id = c("scenario_A", "scenario_B"),
    exchange_reaction_id = c("EX_glucose", "EX_glucose"),
    condition = c("all", "all"),
    lb = c(-0.5, -0.5),
    ub = c(1000, 1000),
    available = c(TRUE, TRUE),
    bound_scope = c("both", "both"),
    stringsAsFactors = FALSE
  )
  cache <- RegCompassR:::.rc_build_full_gem_cache_core(
    gem = gem,
    dirs = data.frame(
      reaction_id = "R1", target_direction = "forward",
      stringsAsFactors = FALSE
    ),
    medium_scenarios = medium,
    cache_dir = tempfile("resolved-medium-cache-")
  )
  summary <- attr(cache, "summary")
  expect_equal(nrow(summary), 2L)
  expect_equal(length(unique(summary$resolved_medium_fingerprint)), 1L)
  expect_equal(summary$medium_fingerprint, summary$resolved_medium_fingerprint)
  expect_true(all(summary$n_changed_bounds_vs_reference >= 1L))
  expect_true(all(!summary$resolved_bounds_identical_to_reference))
  expect_identical(
    attr(cache, "medium_fingerprint_semantics"),
    "base_GEM_identity_plus_canonical_final_reaction_bounds"
  )
})

test_that("resolved fingerprint changes when final bounds change", {
  gem <- make_resolved_medium_test_gem()
  first <- RegCompassR:::rc_build_full_gem(
    gem,
    data.frame(
      medium_scenario_id = "A", exchange_reaction_id = "EX_glucose",
      condition = "all", lb = -0.5, ub = 1000, available = TRUE,
      bound_scope = "both", stringsAsFactors = FALSE
    )
  )
  second <- RegCompassR:::rc_build_full_gem(
    gem,
    data.frame(
      medium_scenario_id = "B", exchange_reaction_id = "EX_glucose",
      condition = "all", lb = -0.2, ub = 1000, available = TRUE,
      bound_scope = "both", stringsAsFactors = FALSE
    )
  )
  expect_false(identical(
    RegCompassR:::.rc_resolved_medium_bounds_fingerprint(gem, first),
    RegCompassR:::.rc_resolved_medium_bounds_fingerprint(gem, second)
  ))
})
