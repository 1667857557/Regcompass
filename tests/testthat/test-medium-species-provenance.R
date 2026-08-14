make_provenance_test_gem <- function(species = "human") {
  reactions <- c(
    "EX_glucose", "EX_lactate", "EX_glutamine", "EX_arginine",
    "EX_leucine", "EX_oxygen", "R1"
  )
  S <- Matrix::Matrix(
    matrix(c(rep(-1, 6), 1), nrow = 1,
           dimnames = list("m_e", reactions)),
    sparse = TRUE
  )
  list(
    S = S,
    lb = stats::setNames(c(rep(-1000, 6), 0), reactions),
    ub = stats::setNames(rep(1000, 7), reactions),
    model_info = list(
      species = species,
      source = if (identical(species, "human")) {
        "SysBioChalmers/Human-GEM"
      } else {
        "SysBioChalmers/Mouse-GEM"
      }
    ),
    reaction_meta = data.frame(
      reaction_id = reactions,
      role = c(rep("exchange", 6), "internal"),
      metabolite_name = c(
        "D-glucose", "L-lactate", "L-glutamine", "L-arginine",
        "L-leucine", "oxygen", NA
      ),
      stringsAsFactors = FALSE
    )
  )
}

test_that("human nutrient challenges are rejected for Mouse-GEM", {
  mouse_gem <- make_provenance_test_gem("mouse")
  human_only <- c(
    "normal_human_plasma", "high_glucose", "low_glucose",
    "high_lactate", "low_lactate", "low_glutamine"
  )
  for (scenario in human_only) {
    expect_error(
      rc_make_medium_scenarios(mouse_gem, scenario = scenario),
      "Human-derived medium scenarios",
      info = scenario
    )
  }
})

test_that("challenge outputs record authoritative background and intervention provenance", {
  gem <- make_provenance_test_gem("human")
  expected <- c(
    high_glucose = "10.1016/j.ygyno.2015.06.036",
    low_glucose = "10.1016/j.ygyno.2015.06.036",
    high_lactate = "10.3389/fonc.2019.01536",
    low_lactate = "10.14814/phy2.70450",
    low_glutamine = "10.1186/s13578-015-0030-1"
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = names(expected),
    species = "human",
    strict_preset_matching = FALSE
  )
  for (scenario in names(expected)) {
    rows <- medium[medium$medium_scenario_id == scenario, , drop = FALSE]
    expect_true(nrow(rows) > 0L)
    expect_true(all(rows$challenge_reference_doi == expected[[scenario]]))
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
    expect_true(all(grepl(
      "authoritative_HPLM_background_plus_named_nutrient_concentration_metadata",
      rows$scenario_construction,
      fixed = TRUE
    )))
    expect_true(all(grepl(
      "no_automatic_concentration_to_flux_mapping",
      rows$scenario_construction,
      fixed = TRUE
    )))
  }
})

test_that("challenge target rows retain concentration provenance without flux conversion", {
  gem <- make_provenance_test_gem("human")
  cases <- list(
    high_glucose = c("glucose", "25", "10.1016/j.ygyno.2015.06.036"),
    low_glucose = c("glucose", "1", "10.1016/j.ygyno.2015.06.036"),
    high_lactate = c("lactate", "20", "10.3389/fonc.2019.01536"),
    low_lactate = c("lactate", "0.5", "10.14814/phy2.70450"),
    low_glutamine = c("glutamine", "0.5", "10.1186/s13578-015-0030-1")
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = names(cases),
    species = "human",
    strict_preset_matching = FALSE
  )
  for (scenario in names(cases)) {
    expected <- cases[[scenario]]
    row <- medium[
      medium$medium_scenario_id == scenario &
        medium$preset_metabolite == expected[[1]],
      , drop = FALSE
    ]
    expect_equal(nrow(row), 1L)
    expect_equal(row$concentration_mM, as.numeric(expected[[2]]))
    expect_equal(row$component_reference_doi, expected[[3]])
    expect_equal(row$challenge_reference_doi, expected[[3]])
    expect_false(row$concentration_used_for_rate_bound)
    expect_equal(
      row$rate_bound_source,
      "background_model_cap_concentration_metadata_only"
    )
    expect_equal(
      row$assumption_level,
      "literature_concentration_metadata_without_flux_conversion"
    )
  }
})
