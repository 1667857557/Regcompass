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
    lb = stats::setNames(
      c(rep(-1000, 8), 0),
      reactions
    ),
    ub = stats::setNames(rep(1000, length(reactions)), reactions),
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
  medium[
    as.character(medium$exchange_reaction_id) == reaction_id,
    ,
    drop = FALSE
  ]
}

test_that("normal human plasma is a cited human uptake-availability preset", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "normal_human_plasma",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  expect_true(all(medium$species == "Homo sapiens"))
  expect_true(all(nzchar(medium$reference_doi)))
  expect_true(all(c(
    "EX_glucose", "EX_lactate", "EX_glutamine",
    "EX_arginine", "EX_oxygen", "EX_urate"
  ) %in% medium$exchange_reaction_id))
  expect_false("EX_unknown" %in% medium$exchange_reaction_id)
  expect_equal(medium_row(medium, "EX_glucose")$concentration_mM, 5)
  expect_equal(medium_row(medium, "EX_lactate")$concentration_mM, 1.5)

  constrained <- rc_apply_medium_constraints(gem, medium)$gem
  expect_lt(unname(constrained$lb["EX_glucose"]), 0)
  expect_lt(unname(constrained$lb["EX_lactate"]), 0)
  expect_equal(unname(constrained$lb["EX_unknown"]), 0)
  expect_equal(unname(constrained$lb["EX_customfuel"]), 0)
})

test_that("glucose presets change glucose uptake without changing lactate", {
  gem <- make_human_medium_test_gem()
  low <- rc_make_medium_scenarios(
    gem, scenario = "low_glucose", exchange_limit = 1,
    strict_preset_matching = FALSE
  )
  normal <- rc_make_medium_scenarios(
    gem, scenario = "normal_human_plasma", exchange_limit = 1,
    strict_preset_matching = FALSE
  )
  high <- rc_make_medium_scenarios(
    gem, scenario = "high_glucose", exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  glucose_lb <- c(
    low = medium_row(low, "EX_glucose")$lb,
    normal = medium_row(normal, "EX_glucose")$lb,
    high = medium_row(high, "EX_glucose")$lb
  )
  expect_true(abs(glucose_lb[["high"]]) > abs(glucose_lb[["normal"]]))
  expect_true(abs(glucose_lb[["normal"]]) > abs(glucose_lb[["low"]]))
  expect_equal(medium_row(low, "EX_glucose")$concentration_mM, 1)
  expect_equal(medium_row(high, "EX_glucose")$concentration_mM, 25)
  expect_equal(
    medium_row(low, "EX_lactate")$lb,
    medium_row(high, "EX_lactate")$lb
  )
})

test_that("lactate presets change lactate uptake without changing glucose", {
  gem <- make_human_medium_test_gem()
  low <- rc_make_medium_scenarios(
    gem, scenario = "low_lactate", exchange_limit = 1,
    strict_preset_matching = FALSE
  )
  normal <- rc_make_medium_scenarios(
    gem, scenario = "normal_human_plasma", exchange_limit = 1,
    strict_preset_matching = FALSE
  )
  high <- rc_make_medium_scenarios(
    gem, scenario = "high_lactate", exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  lactate_lb <- c(
    low = medium_row(low, "EX_lactate")$lb,
    normal = medium_row(normal, "EX_lactate")$lb,
    high = medium_row(high, "EX_lactate")$lb
  )
  expect_true(abs(lactate_lb[["high"]]) > abs(lactate_lb[["normal"]]))
  expect_true(abs(lactate_lb[["normal"]]) > abs(lactate_lb[["low"]]))
  expect_equal(medium_row(low, "EX_lactate")$concentration_mM, 0.5)
  expect_equal(medium_row(high, "EX_lactate")$concentration_mM, 20)
  expect_equal(
    medium_row(low, "EX_glucose")$lb,
    medium_row(high, "EX_glucose")$lb
  )
})

test_that("RPMI-1640 enables basal formulation nutrients but not plasma-only rows", {
  gem <- make_human_medium_test_gem()
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "rpmi1640",
    exchange_limit = 1,
    strict_preset_matching = FALSE
  )

  expect_true(all(c(
    "EX_glucose", "EX_glutamine", "EX_arginine", "EX_oxygen"
  ) %in% medium$exchange_reaction_id))
  expect_false("EX_lactate" %in% medium$exchange_reaction_id)
  expect_false("EX_urate" %in% medium$exchange_reaction_id)
  expect_equal(
    medium_row(medium, "EX_glucose")$concentration_mM,
    11.111111, tolerance = 1e-6
  )
  expect_match(
    medium_row(medium, "EX_glucose")$reference_label,
    "RPMI|Moore"
  )
})

test_that("custom metabolite presets map availability and relative uptake", {
  gem <- make_human_medium_test_gem()
  compounds <- data.frame(
    metabolite_name = "customfuel",
    metabolite_pattern = "customfuel",
    available = TRUE,
    concentration_mM = 3,
    uptake_fraction = 0.3,
    target_exchange_flag = TRUE,
    required_match = TRUE,
    reference_label = "Example human measurement",
    reference_doi = "10.0000/example",
    reference_pmid = "12345678",
    stringsAsFactors = FALSE
  )

  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "custom",
    custom_metabolites = compounds,
    exchange_limit = 2
  )
  row <- medium_row(medium, "EX_customfuel")

  expect_equal(row$lb, -0.6)
  expect_equal(row$ub, 2)
  expect_equal(row$concentration_mM, 3)
  expect_equal(row$reference_doi, "10.0000/example")

  constrained <- rc_apply_medium_constraints(gem, medium)$gem
  expect_equal(unname(constrained$lb["EX_customfuel"]), -0.6)
  expect_equal(unname(constrained$lb["EX_unknown"]), 0)
})

test_that("legacy human-like aliases map to the published presets", {
  gem <- make_human_medium_test_gem()
  expect_warning(
    blood <- rc_make_medium_scenarios(
      gem, scenario = "blood_like", exchange_limit = 1,
      strict_preset_matching = FALSE
    ),
    "blood_like -> normal_human_plasma"
  )
  expect_true(all(blood$medium_scenario_id == "normal_human_plasma"))

  expect_warning(
    culture <- rc_make_medium_scenarios(
      gem, scenario = "culture_like", exchange_limit = 1,
      strict_preset_matching = FALSE
    ),
    "culture_like -> rpmi1640"
  )
  expect_true(all(culture$medium_scenario_id == "rpmi1640"))
})
