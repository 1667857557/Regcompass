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

test_that("normal human plasma uses authoritative physiological-medium sources", {
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
  for (doi in c(
    "10.1016/j.cell.2017.03.023",
    "10.1016/j.cmet.2021.02.005",
    "10.1126/sciadv.aau7314"
  )) {
    expect_true(all(grepl(doi, medium$reference_doi, fixed = TRUE)))
  }
  expect_false(any(grepl(
    "10.1371/journal.pone.0016957",
    medium$reference_doi,
    fixed = TRUE
  )))
  expect_true(all(
    medium$medium_background_id == "authoritative_HPLM_2017_2021"
  ))
  expect_true(all(
    medium$scenario_construction ==
      "authoritative_HPLM_composition_without_cross_study_averaging"
  ))
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
  expect_false("EX_oxygen" %in% medium$exchange_reaction_id)
  expect_false("EX_unknown" %in% medium$exchange_reaction_id)
  expect_false(any(medium$concentration_used_for_rate_bound))
})

test_that("all culture challenges use the authoritative HPLM background", {
  gem <- make_human_medium_test_gem()
  scenarios <- c(
    "high_glucose", "low_glucose", "high_lactate", "low_lactate",
    "low_glutamine"
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = scenarios,
    species = "human",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  for (scenario_id in scenarios) {
    rows <- medium[medium$medium_scenario_id == scenario_id, , drop = FALSE]
    expect_true(nrow(rows) > 0L)
    expect_true(all(
      rows$medium_background_id == "authoritative_HPLM_2017_2021"
    ))
    expect_true(all(
      rows$scenario_construction ==
        "authoritative_HPLM_background_plus_named_concentration_provenance"
    ))
    expect_true(all(grepl(
      "10.1016/j.cell.2017.03.023",
      rows$background_reference_doi,
      fixed = TRUE
    )))
    expect_true(all(grepl(
      "10.1016/j.cmet.2021.02.005",
      rows$background_reference_doi,
      fixed = TRUE
    )))
    expect_true(all(
      rows$background_validation_reference_doi == "10.1126/sciadv.aau7314"
    ))
    expect_false(any(grepl(
      "10.1001/jama.1967.03120080053007|10.1016/0042-6822",
      rows$background_reference_doi
    )))
    expect_true(all(c(
      "EX_glucose", "EX_glutamine", "EX_arginine", "EX_leucine"
    ) %in% rows$exchange_reaction_id))
    expect_false(any(rows$concentration_used_for_rate_bound))
  }
})

test_that("challenge concentrations override only the named nutrient metadata", {
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

  expect_equal(
    medium_row(medium, "high_glucose", "EX_glucose")$concentration_mM,
    25
  )
  expect_equal(
    medium_row(medium, "low_glucose", "EX_glucose")$concentration_mM,
    1
  )
  expect_equal(
    medium_row(medium, "high_lactate", "EX_lactate")$concentration_mM,
    20
  )
  expect_equal(
    medium_row(medium, "low_lactate", "EX_lactate")$concentration_mM,
    0.5
  )
  expect_equal(
    medium_row(medium, "low_glutamine", "EX_glutamine")$concentration_mM,
    0.5
  )
  expect_equal(
    medium_row(medium, "high_glucose", "EX_glucose")$challenge_reference_doi,
    "10.1016/j.ygyno.2015.06.036"
  )
  expect_equal(
    medium_row(medium, "high_lactate", "EX_lactate")$challenge_reference_doi,
    "10.3389/fonc.2019.01536"
  )
  expect_equal(
    medium_row(medium, "low_lactate", "EX_lactate")$challenge_reference_doi,
    "10.14814/phy2.70450"
  )
  expect_equal(
    medium_row(medium, "low_glutamine", "EX_glutamine")$challenge_reference_doi,
    "10.1186/s13578-015-0030-1"
  )
})

test_that("published concentrations do not create implicit uptake ratios", {
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
  for (scenario_id in c(
    "high_glucose", "low_glucose", "high_lactate", "low_lactate",
    "low_glutamine"
  )) {
    target <- if (grepl("glucose", scenario_id)) {
      "EX_glucose"
    } else if (grepl("lactate", scenario_id)) {
      "EX_lactate"
    } else {
      "EX_glutamine"
    }
    row <- medium_row(medium, scenario_id, target)
    expect_equal(row$uptake_fraction, 1)
    expect_false(row$concentration_used_for_rate_bound)
  }
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

test_that("user-defined media remain unrestricted by publication metadata", {
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
    "authoritative_journal_composition_with_explicit_overrides"
  )
})

test_that("scenario NULL supports explicit custom metabolite uptake assumptions", {
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
