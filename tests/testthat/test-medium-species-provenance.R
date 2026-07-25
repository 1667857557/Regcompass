test_that("human nutrient challenges are rejected for Mouse-GEM", {
  mouse_gem <- list(
    model_info = list(
      species = "mouse",
      source = "SysBioChalmers/Mouse-GEM"
    )
  )
  human_only <- c(
    "normal_human_plasma", "high_glucose", "low_glucose",
    "high_lactate", "low_lactate", "low_glutamine"
  )
  for (scenario in human_only) {
    expect_error(
      rc_make_medium_scenarios(mouse_gem, scenario = scenario),
      "Human-derived medium presets"
    )
  }
})

test_that("reference catalog records explicit species and study DOIs", {
  refs <- .rc_medium_reference_catalog()
  human_challenges <- c(
    high_glucose = "10.1016/j.ygyno.2015.06.036",
    low_glucose = "10.1016/j.ygyno.2015.06.036",
    high_lactate = "10.3389/fonc.2019.01536",
    low_lactate = "10.14814/phy2.70450",
    low_glutamine = "10.1186/s13578-015-0030-1"
  )
  rows <- refs[match(names(human_challenges), refs$preset_id), , drop = FALSE]
  expect_true(all(rows$species == "Homo sapiens"))
  expect_equal(rows$reference_doi, unname(human_challenges))
  expect_equal(
    refs$species[refs$preset_id %in% c("rpmi1640", "dmem_high_glucose")],
    rep("not species-specific", 2)
  )
})

test_that("human challenge target rows retain their own provenance", {
  cases <- list(
    high_glucose = c("glucose", "25", "10.1016/j.ygyno.2015.06.036"),
    low_glucose = c("glucose", "1", "10.1016/j.ygyno.2015.06.036"),
    high_lactate = c("lactate", "20", "10.3389/fonc.2019.01536"),
    low_lactate = c("lactate", "0.5", "10.14814/phy2.70450"),
    low_glutamine = c("glutamine", "0.05", "10.1186/s13578-015-0030-1")
  )
  for (scenario in names(cases)) {
    expected <- cases[[scenario]]
    catalog <- .rc_medium_catalog(scenario, "human")
    row <- catalog[catalog$metabolite_name == expected[[1]], , drop = FALSE]
    expect_equal(row$concentration_mM, as.numeric(expected[[2]]))
    expect_equal(row$component_reference_doi, expected[[3]])
  }
})
