make_human_medium_test_gem <- function() {
  reactions <- c(
    "EX_glucose", "EX_lactate", "EX_glutamine", "EX_arginine",
    "EX_leucine", "EX_oxygen", "EX_urate", "EX_customfuel",
    "EX_unknown", "R1"
  )
  S <- Matrix::Matrix(
    matrix(
      c(rep(-1, 9), 1),
      nrow = 1,
      dimnames = list("m_e", reactions)
    ),
    sparse = TRUE
  )
  list(
    S = S,
    lb = stats::setNames(c(rep(-1000, 9), 0), reactions),
    ub = stats::setNames(rep(1000, length(reactions)), reactions),
    model_info = list(species = "human"),
    reaction_meta = data.frame(
      reaction_id = reactions,
      role = c(rep("exchange", 9), "internal"),
      metabolite_name = c(
        "D-glucose", "L-lactate", "L-glutamine", "L-arginine",
        "L-leucine", "oxygen", "urate", "customfuel",
        "unlisted nutrient", NA
      ),
      stringsAsFactors = FALSE
    )
  )
}

medium_row <- function(medium, scenario_id, reaction_id) {
  medium[
    as.character(medium$medium_scenario_id) == scenario_id &
      as.character(medium$exchange_reaction_id) == reaction_id,
    , drop = FALSE
  ]
}

test_that("normal human plasma integrates cited HPLM and plasma evidence", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "normal_human_plasma",
    species = "human",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  expect_true(nrow(medium) > 0L)
  expect_true(all(medium$medium_scenario_id == "normal_human_plasma"))
  expect_true(all(medium$species == "Homo sapiens"))
  expect_true(all(grepl("10.1016/j.cell.2017.03.023", medium$reference_doi,
                        fixed = TRUE)))
  expect_true(all(grepl("10.1371/journal.pone.0016957", medium$reference_doi,
                        fixed = TRUE)))
  expect_equal(
    medium_row(medium, "normal_human_plasma", "EX_glucose")$concentration_mM,
    5
  )
  expect_equal(
    medium_row(medium, "normal_human_plasma", "EX_lactate")$concentration_mM,
    1.6
  )
  expect_equal(
    medium_row(medium, "normal_human_plasma", "EX_glutamine")$concentration_mM,
    0.55000347
  )
  expect_false("EX_unknown" %in% medium$exchange_reaction_id)
})

test_that("glucose challenge presets retain culture nutrients and override glucose", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = c("high_glucose", "low_glucose"),
    species = "human",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  expect_equal(
    medium_row(medium, "high_glucose", "EX_glucose")$concentration_mM,
    25
  )
  expect_equal(
    medium_row(medium, "low_glucose", "EX_glucose")$concentration_mM,
    1
  )
  for (scenario_id in c("high_glucose", "low_glucose")) {
    expect_true(all(c(
      "EX_glucose", "EX_glutamine", "EX_arginine", "EX_leucine",
      "EX_oxygen"
    ) %in% medium[medium$medium_scenario_id == scenario_id,
                   "exchange_reaction_id", drop = TRUE]))
    rows <- medium[medium$medium_scenario_id == scenario_id, , drop = FALSE]
    expect_true(all(rows$medium_background_id ==
                      "published_RPMI_DMEM_nutrient_union"))
    expect_true(all(rows$scenario_construction ==
                      "published_background_plus_named_nutrient_override"))
    expect_true(all(grepl("10.1016/j.ygyno.2015.06.036",
                          rows$challenge_reference_doi, fixed = TRUE)))
  }
})

test_that("lactate challenges retain their published basal environments", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = c("high_lactate", "low_lactate"),
    species = "human",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  high <- medium_row(medium, "high_lactate", "EX_lactate")
  low <- medium_row(medium, "low_lactate", "EX_lactate")
  expect_equal(high$concentration_mM, 20)
  expect_equal(low$concentration_mM, 0.5)
  expect_equal(high$medium_background_id,
               "published_DMEM_nutrient_background")
  expect_equal(low$medium_background_id,
               "published_plasma_like_nutrient_background")
  expect_equal(high$challenge_reference_doi,
               "10.3389/fonc.2019.01536")
  expect_equal(low$challenge_reference_doi,
               "10.14814/phy2.70450")
  expect_true("EX_glutamine" %in%
                medium[medium$medium_scenario_id == "high_lactate",
                       "exchange_reaction_id", drop = TRUE])
  expect_true("EX_arginine" %in%
                medium[medium$medium_scenario_id == "low_lactate",
                       "exchange_reaction_id", drop = TRUE])
})

test_that("low glutamine uses the Methods-defined 0.5 mM condition", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "low_glutamine",
    species = "human",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )
  row <- medium_row(medium, "low_glutamine", "EX_glutamine")
  expect_equal(row$concentration_mM, 0.5)
  expect_equal(row$challenge_reference_doi,
               "10.1186/s13578-015-0030-1")
  expect_true("EX_glucose" %in% medium$exchange_reaction_id)
  expect_true("EX_leucine" %in% medium$exchange_reaction_id)
})

test_that("target concentration caps remain explicit sensitivity assumptions", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = c(
      "high_glucose", "low_glucose", "high_lactate", "low_lactate",
      "low_glutamine"
    ),
    species = "human",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )
  row_for <- function(scenario_id, reaction_id) {
    medium_row(medium, scenario_id, reaction_id)
  }
  expect_equal(row_for("high_glucose", "EX_glucose")$uptake_fraction, 1)
  expect_equal(row_for("low_glucose", "EX_glucose")$uptake_fraction, 1 / 25)
  expect_equal(row_for("high_lactate", "EX_lactate")$uptake_fraction, 1)
  expect_equal(row_for("low_lactate", "EX_lactate")$uptake_fraction, 0.5 / 20)
  expect_equal(row_for("low_glutamine", "EX_glutamine")$uptake_fraction,
               0.5 / 4)
  expect_true(all(
    medium$concentration_used_for_rate_bound[
      medium$target_exchange_flag %in% TRUE
    ]
  ))
})

test_that("technical and hidden aliases are rejected", {
  gem <- make_human_medium_test_gem()
  retired <- c(
    "physiologic", "cantor2017_hplm", "rpmi1640", "dmem_high_glucose",
    "minimal", "compass_model_bounds", "permissive_all_exchange"
  )
  for (scenario in retired) {
    expect_error(
      rc_make_medium_scenarios(
        gem,
        scenario = scenario,
        species = "human",
        strict_preset_matching = FALSE
      ),
      "Unsupported biological medium scenario",
      info = scenario
    )
  }
})

test_that("user-defined media do not require publication metadata", {
  gem <- make_human_medium_test_gem()
  custom <- data.frame(
    medium_scenario_id = "measured_custom",
    exchange_reaction_id = "EX_customfuel",
    lb = -0.4,
    ub = 1,
    available = TRUE,
    stringsAsFactors = FALSE
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "custom",
    species = "human",
    custom_medium = custom
  )
  row <- medium_row(medium, "measured_custom", "EX_customfuel")
  expect_equal(row$lb, -0.4)
  expect_equal(row$evidence_source, "user_supplied_custom_medium")
  expect_identical(
    attr(medium, "medium_policy"),
    "published_plasma_or_culture_background_with_explicit_overrides"
  )
})

test_that("scenario NULL supports custom-only metabolite composition", {
  gem <- make_human_medium_test_gem()
  custom <- data.frame(
    metabolite_name = "customfuel",
    available = TRUE,
    concentration_mM = 3,
    uptake_fraction = 0.3,
    target_exchange_flag = TRUE,
    required_match = TRUE,
    stringsAsFactors = FALSE
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = NULL,
    species = "human",
    custom_metabolites = custom,
    exchange_limit = 2
  )
  row <- medium_row(medium, "custom", "EX_customfuel")
  expect_equal(row$lb, -0.6)
  expect_equal(row$concentration_mM, 3)
})

test_that("built-in and custom environments can be returned together", {
  gem <- make_human_medium_test_gem()
  custom <- data.frame(
    medium_scenario_id = "measured_custom",
    exchange_reaction_id = "EX_customfuel",
    lb = -0.4,
    ub = 1,
    available = TRUE,
    stringsAsFactors = FALSE
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "normal_human_plasma",
    species = "human",
    custom_medium = custom,
    strict_preset_matching = FALSE
  )
  expect_setequal(
    unique(medium$medium_scenario_id),
    c("normal_human_plasma", "measured_custom")
  )
})