make_human_medium_test_gem <- function() {
  reactions <- c(
    "EX_glucose", "EX_lactate", "EX_glutamine", "EX_arginine",
    "EX_oxygen", "EX_urate", "EX_customfuel", "EX_unknown", "R1"
  )
  S <- Matrix::Matrix(
    matrix(
      c(-1, -1, -1, -1, -1, -1, -1, -1, 1),
      nrow = 1,
      dimnames = list("m_e", reactions)
    ),
    sparse = TRUE
  )
  list(
    S = S,
    lb = stats::setNames(c(rep(-1000, 8), 0), reactions),
    ub = stats::setNames(rep(1000, length(reactions)), reactions),
    model_info = list(species = "human"),
    reaction_meta = data.frame(
      reaction_id = reactions,
      role = c(rep("exchange", 8), "internal"),
      metabolite_name = c(
        "D-glucose", "L-lactate", "L-glutamine", "L-arginine",
        "oxygen", "urate", "customfuel", "unlisted nutrient", NA
      ),
      stringsAsFactors = FALSE
    )
  )
}

medium_row <- function(medium, reaction_id) {
  medium[as.character(medium$exchange_reaction_id) == reaction_id, , drop = FALSE]
}

test_that("Cantor 2017 HPLM is the publication-bound built-in scenario", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "cantor2017_hplm",
    species = "human",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  expect_true(nrow(medium) > 0L)
  expect_true(all(medium$medium_scenario_id == "cantor2017_hplm"))
  expect_true(all(medium$species == "Homo sapiens"))
  expect_true(all(medium$reference_doi == "10.1016/j.cell.2017.03.023"))
  expect_true(all(medium$component_reference_doi ==
                    "10.1016/j.cell.2017.03.023"))
  expect_true(all(medium$concentration_basis ==
                    "Cantor_2017_HPLM_formulation"))
  expect_true(all(medium$uptake_fraction == 1))
  expect_false(any(medium$target_exchange_flag))
  expect_false(any(medium$concentration_used_for_rate_bound))
  expect_identical(
    attr(medium, "medium_policy"),
    "published_paper_bound_presets_only"
  )

  expect_equal(medium_row(medium, "EX_glucose")$concentration_mM, 5)
  expect_equal(medium_row(medium, "EX_lactate")$concentration_mM, 1.6)
  expect_equal(
    medium_row(medium, "EX_glutamine")$concentration_mM,
    0.55000347
  )
  expect_false("EX_oxygen" %in% medium$exchange_reaction_id)
  expect_false("EX_unknown" %in% medium$exchange_reaction_id)
})

test_that("published concentrations do not become uptake-rate scalers", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "cantor2017_hplm",
    species = "human",
    exchange_limit = 2,
    strict_preset_matching = FALSE
  )
  expect_equal(medium_row(medium, "EX_glucose")$lb, -2)
  expect_equal(medium_row(medium, "EX_lactate")$lb, -2)
  expect_equal(medium_row(medium, "EX_glutamine")$lb, -2)
  expect_error(
    rc_make_medium_scenarios(
      gem,
      scenario = "cantor2017_hplm",
      species = "human",
      uptake_scale = 0.5,
      strict_preset_matching = FALSE
    ),
    "uptake_scale = 1"
  )
})

test_that("insufficiently documented legacy scenario identifiers are rejected", {
  gem <- make_human_medium_test_gem()
  retired <- c(
    "physiologic", "normal_human_plasma", "mouse_plasma", "rpmi1640",
    "dmem_high_glucose", "high_glucose", "low_glucose", "high_lactate",
    "low_lactate", "low_glutamine", "minimal", "compass_model_bounds",
    "permissive_all_exchange", "custom"
  )
  for (scenario in retired) {
    expect_error(
      rc_make_medium_scenarios(
        gem,
        scenario = scenario,
        species = "human",
        strict_preset_matching = FALSE
      ),
      "Unsupported or insufficiently documented medium scenario",
      info = scenario
    )
  }
})

test_that("custom media require published-paper provenance", {
  gem <- make_human_medium_test_gem()
  missing_reference <- data.frame(
    medium_scenario_id = "measured_custom",
    exchange_reaction_id = "EX_customfuel",
    lb = -0.4,
    ub = 1,
    available = TRUE,
    stringsAsFactors = FALSE
  )
  expect_error(
    rc_make_medium_scenarios(
      gem,
      scenario = NULL,
      species = "human",
      custom_medium = missing_reference
    ),
    "publication provenance columns"
  )

  invalid_doi <- transform(
    missing_reference,
    reference_label = "Unpublished example",
    reference_doi = "not-a-doi"
  )
  expect_error(
    rc_make_medium_scenarios(
      gem,
      scenario = NULL,
      species = "human",
      custom_medium = invalid_doi
    ),
    "valid published-paper"
  )
})

test_that("DOI-cited custom media are accepted without a custom scenario alias", {
  gem <- make_human_medium_test_gem()
  custom <- data.frame(
    medium_scenario_id = "published_custom_2024",
    exchange_reaction_id = "EX_customfuel",
    lb = -0.4,
    ub = 1,
    available = TRUE,
    reference_label = "Example et al., Example Journal 2024",
    reference_doi = "10.1234/example.2024.1",
    stringsAsFactors = FALSE
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = NULL,
    species = "human",
    custom_medium = custom
  )
  row <- medium_row(medium, "EX_customfuel")
  expect_equal(row$medium_scenario_id, "published_custom_2024")
  expect_equal(row$reference_doi, "10.1234/example.2024.1")
  expect_equal(row$lb, -0.4)
  expect_identical(
    attr(medium, "medium_policy"),
    "published_paper_bound_presets_only"
  )
})

test_that("built-in and DOI-cited custom environments can be returned together", {
  gem <- make_human_medium_test_gem()
  custom <- data.frame(
    medium_scenario_id = "published_custom_2024",
    exchange_reaction_id = "EX_customfuel",
    lb = -0.4,
    ub = 1,
    available = TRUE,
    reference_label = "Example et al., Example Journal 2024",
    reference_doi = "10.1234/example.2024.1",
    stringsAsFactors = FALSE
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "cantor2017_hplm",
    species = "human",
    custom_medium = custom,
    strict_preset_matching = FALSE
  )
  expect_setequal(
    unique(medium$medium_scenario_id),
    c("cantor2017_hplm", "published_custom_2024")
  )
})

test_that("mouse analyses require a DOI-cited custom environment", {
  gem <- make_human_medium_test_gem()
  gem$model_info$species <- "mouse"
  expect_error(
    rc_make_medium_scenarios(
      gem,
      scenario = "cantor2017_hplm",
      species = "mouse",
      strict_preset_matching = FALSE
    ),
    "requires a Human-GEM"
  )
  expect_error(
    rc_make_medium_scenarios(
      gem,
      scenario = NULL,
      species = "mouse"
    ),
    "Provide `scenario"
  )
})